#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <unordered_map>
#include <queue>
#include <chrono>
#include <cstdlib>

#define EARTH_RADIUS_KM 6371.0
#define DISTANCE_THRESHOLD_KM 100.0
#define CITY_LIMIT 100000

using namespace std;
using namespace std::chrono;

struct City {
    string name;
    string state;
    double lat;
    double lng;    
};

string removeQuotesAndTrim(const string& input) {
    string result = input;

    // Remove quotes
    result.erase(std::remove(result.begin(), result.end(), '\"'), result.end());

    // Trim spaces from the beginning
    size_t first = result.find_first_not_of(' ');
    if (first == string::npos) {
        return ""; // All spaces string or empty string
    }

    // Trim spaces from the end
    size_t last = result.find_last_not_of(' ');
    return result.substr(first, (last - first + 1));
}

double toRadians(double degree) {
    return degree * (M_PI / 180.0);
}

double calculateDistance(double lat1, double lng1, double lat2, double lng2) {
    // Haversine formula to calculate the distance between two points on the Earth
    double dLat = toRadians(lat2 - lat1);
    double dLng = toRadians(lng2 - lng1);
    lat1 = toRadians(lat1);
    lat2 = toRadians(lat2);

    double a = sin(dLat/2) * sin(dLat/2) + sin(dLng/2) * sin(dLng/2) * cos(lat1) * cos(lat2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));

    return EARTH_RADIUS_KM * c;
}

void printGraph(const unordered_map<string, vector<pair<string, double>>>& graph) {
    for (const auto& entry : graph) {
        std::cout << entry.first << ": ";
        for (const auto& neighbor : entry.second) {
            std::cout << "(" << neighbor.first << ", " << neighbor.second << " km) ";
        }
        std::cout << endl;
    }
}

void outputGraphToFile(const unordered_map<string, vector<pair<string, double>>>& graph, const string& filename) {
    ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        cerr << "Error opening output file" << endl;
        return;
    }

    // Output header line
    outputFile << "From,To,Distance" << endl;

    for (const auto& cityPair : graph) {
        for (const auto& connection : cityPair.second) {
            
            // Use concatenated city names
            outputFile << cityPair.first << "," 
                       << connection.first << ","  
                       << connection.second << endl;
        }
    }

    outputFile.close();
    std::cout << "Graph data written to " << filename << endl;
}

vector<City> readCitiesFromFile(const string& filename) {
    vector<City> cities;
    ifstream file(filename);
    string line;

    if (!file.is_open()) {
        cerr << "Error opening file" << endl;
        throw runtime_error("File opening failed");
    }

    getline(file, line); // Skip the first line (header)
    std::cout << "Creating cities vector..." << endl;

    while (getline(file, line)) {
        stringstream ss(line);
        string city, state, lat, lng, token;
        int column = 0;
        while (getline(ss, token, ',')) {
            column++;
            if (column == 1) {
                city = removeQuotesAndTrim(token);
            } else if (column == 3) {
                state = removeQuotesAndTrim(token);
            } else if (column == 7) {
                lat = removeQuotesAndTrim(token);
            } else if (column == 8) {
                lng = removeQuotesAndTrim(token);
            }
        }

        if (state == "AK" || state == "HI") { // Skip cities in Alaska and Hawaii
            continue;
        }

        try {
            City newCity = {city, state, stod(lat), stod(lng)};
            cities.push_back(newCity);
        } catch (const std::invalid_argument& e) {
           // cerr << "Invalid data for city: " << city << ". Skipping." << endl;
        }

        // Use only first CITY_LIMIT cities in CSV for testing
        if (cities.size() >= CITY_LIMIT) {
            break;
        }
    }

    file.close();
    return cities;
}

string getCityState(const City& city) {
    return city.name + " " + city.state;
}

unordered_map<string, vector<pair<string, double>>> createCityGraph(const vector<City>& cities) {
    unordered_map<string, vector<pair<string, double>>> graph;

    for (const auto& city1 : cities) {
        vector<pair<string, double>> cityDistances;
        for (const auto& city2 : cities) {
            if (city1.name != city2.name) {
                double distance = calculateDistance(city1.lat, city1.lng, city2.lat, city2.lng);
                if (distance <= DISTANCE_THRESHOLD_KM) {  
                    string city1Name = getCityState(city1); 
                    string city2Name = getCityState(city2);
                    
                    cityDistances.push_back(make_pair(city2Name, distance));
                }
            }
        }
        // Use concatenated name for graph key
        graph[getCityState(city1)] = cityDistances; 
    }

    return graph;
}

unordered_map<string, string> dijkstra(unordered_map<string, vector<pair<string, double>>>& graph, string start_city) {
    unordered_map<string, string> predecessor_list;
    priority_queue<pair<double, string>, vector<pair<double, string>>, greater<pair<double, string>>> minHeap;

    // Keep track of distances from the start city to each city in the graph
    unordered_map<string, double> dist_vec;

    // Initialize distances to infinity for all cities except the start city
    for (const auto& entry : graph) {
        dist_vec[entry.first] = numeric_limits<double>::max();
    }

    // Distance from start_city to itself is 0
    dist_vec[start_city] = 0;
    minHeap.push({0, start_city});

    // Main Dijkstra's algorithm loop:
    // - Get the city with the minimum distance from the priority queue
    // - Iterate through neighbors (connected cities) of the current city
    // - If a shorter path to the neighbor through the current city is found:
    //     - Update the distance to the neighbor
    //     - Set the predecessor of the neighbor to be the current city
    //     - Push the neighbor into the priority queue with its updated distance
    while (!minHeap.empty()) {
        string curr_city = minHeap.top().second;
        minHeap.pop();

        for (const auto& neighbor : graph[curr_city]) {
            string neighbor_city = neighbor.first;
            double edge_weight = neighbor.second;

            if (dist_vec[curr_city] + edge_weight < dist_vec[neighbor_city]) {
                dist_vec[neighbor_city] = dist_vec[curr_city] + edge_weight;
                predecessor_list[neighbor_city] = curr_city;
                minHeap.push({dist_vec[neighbor_city], neighbor_city});
            }
        }
    }

    // Return the map of predecessors, represents the shortest paths from start_city to all other cities
    return predecessor_list;
}


void printPathVector(vector<string> path_vec) {
    for(unsigned int i = 0; i < path_vec.size(); i++) {
        if(i == path_vec.size() - 1) {
            std::cout << path_vec[i] << endl;
            break;
        } else
            std::cout << path_vec[i] << ", ";
    }
    return;
}

vector<string> dijkstraPathFinding(unordered_map<string, vector<pair<string, double>>> graph, string start_city, string end_city) {
    unordered_map<string, string> predecessor_list = dijkstra(graph, start_city);
    vector<string> path_vec; // 0th index is start city, final index is end_city
    path_vec.push_back(end_city);
    while(predecessor_list[end_city] != start_city) {
        path_vec.push_back(predecessor_list[end_city]);
        end_city = predecessor_list[end_city];
    }
    path_vec.push_back(start_city);
    reverse(path_vec.begin(), path_vec.end());
    std::cout << "Dijkstra's path: ";
    printPathVector(path_vec);
    return path_vec;
}

struct compare {
    bool operator()(const pair<string, double>& l, const pair<string, double>& r) {
        return l.second > r.second;
    }
};

unordered_map<string, string> AStar(
    unordered_map<string, vector<pair<string, double>>> graph,
    string start_city, string end_city, unordered_map<string, City> city_map) {

    priority_queue<pair<string, double>, vector<pair<string, double>>, compare> Q;
    unordered_map<string, bool> closed;
    unordered_map<string, double> dist_vec;
    unordered_map<string, string> predecessor_list;

    for (const auto& pair : graph) {
        dist_vec[pair.first] = INT_MAX;
        predecessor_list[pair.first] = "";
        closed[pair.first] = false;
    }

    dist_vec[start_city] = 0;
    predecessor_list[start_city] = start_city;
    Q.push({start_city, 0.0});

    while (!Q.empty()) {
        string current = Q.top().first;
        Q.pop();
        closed[current] = true;

        if (current == end_city) break;

        for (const auto& neighbor : graph[current]) {
            if (closed[neighbor.first]) continue;

            double newDist = dist_vec[current] + neighbor.second;
            if (newDist < dist_vec[neighbor.first]) {
                dist_vec[neighbor.first] = newDist;
                predecessor_list[neighbor.first] = current;
                double heuristic = calculateDistance(city_map[neighbor.first].lat, city_map[end_city].lat, city_map[neighbor.first].lng, city_map[end_city].lng);
                Q.push({neighbor.first, newDist + heuristic});
            }
        }
    }

    return predecessor_list;
}


unordered_map<string, City> map_cities(vector<City> cities){
    unordered_map<string, City> city_map;
    for(auto i = cities.begin(); i != cities.end(); i++){
        city_map[i->name] = *i;
    }
    return city_map;
}

vector<string> AStarPathFinding(unordered_map<string, vector<pair<string, double>>> graph, string start_city, string end_city, vector<City> cities) {
    unordered_map<string, City> city_map = map_cities(cities);
    unordered_map<string, string> predecessor_list = AStar(graph, start_city, end_city, city_map);
    vector<string> path_vec; // 0th index is start city, final index is end_city
    path_vec.push_back(end_city);
    while(predecessor_list[end_city] != start_city) {
        path_vec.push_back(predecessor_list[end_city]);
        end_city = predecessor_list[end_city];
    }
    path_vec.push_back(start_city);
    reverse(path_vec.begin(), path_vec.end());
    std::cout << "A* path: ";
    printPathVector(path_vec);
    return path_vec;
}

unordered_map<string, vector<pair<string, double>>> createGraphFromInputFile(const string &filename) {
    ifstream file(filename);
    string line;
    unordered_map<string, vector<pair<string, double>>> graph;
    if (!file.is_open()) {
        cerr << "Error opening file" << endl;
        throw runtime_error("File opening failed");
    }
    getline(file, line); // Skip the first line (header)
    std::cout << "Creating cities graph from input file..." << endl;
    int cnt = 0;
    int done = 0;
    while (getline(file, line)) {
        cnt++;
        if(cnt % 600000 == 0) {
            done+=10;
            std::cout << done << " percent done with creating graph from input" << endl;
        }
        stringstream ss(line);
        string from, to, dist, token;
        int column = 0;
        while (getline(ss, token, ',')) {
            column++;
            if (column == 1) {
                from = removeQuotesAndTrim(token);
            } else if (column == 2) {
                to = removeQuotesAndTrim(token);
            } else if (column == 3) {
                dist = removeQuotesAndTrim(token);
            }
        }
        try {
            graph[from].push_back(make_pair(to, stod(dist)));
        } catch (const std::invalid_argument& e) {
            cerr << "Invalid data for cities: " << from << ", " << to << ". Skipping." << endl;
        }
    }
    std::cout << "Finished creating graph from input!" << endl << endl;
    file.close();
    return graph;
}

int main() {
    vector<City> cities = readCitiesFromFile("input/cities.csv");
    // auto graph = createCityGraph(cities);
    // outputGraphToFile(graph, "full_adjacency_list.csv");
    auto graph = createGraphFromInputFile("input/full_adjacency_list.csv");
    string startCity, endCity;

    std::cout << "Pathfinding Application - Returns the path between two cities and the time it took to calculate the result for each algorithm..." << endl;
    while (true) {
        std::cout << "Enter start city (e.g. New York NY) (or QUIT to exit): ";
        getline(cin, startCity);
        if (startCity == "QUIT") break;

        std::cout << "Enter end city (e.g. Los Angeles CA): "; 
        getline(cin, endCity);

        std::cout << endl;

        try {
            struct output_information {
                string start_city;
                string end_city;
                string dijkstra_path;
                string astar_path;
                double dijkstra_duration;
                double astar_duration;
            };
            struct output_information outinfo;
            outinfo.start_city = startCity;
            outinfo.end_city = endCity;

            std::cout << "Calculating path using Dijkstra's algorithm..." << endl;
            auto start = high_resolution_clock::now();
            auto dijksta_res = dijkstraPathFinding(graph, startCity, endCity);
            auto stop = high_resolution_clock::now();
            auto algo_duration = duration<double>(stop - start).count();
            std::cout << "Time taken by Dijkstra's algorithm: " << algo_duration << " seconds" << endl << endl;
            outinfo.dijkstra_duration = algo_duration;

            std::cout << "Calculating path using A* algorithm..." << endl;
            start = high_resolution_clock::now();
            auto astar_res = AStarPathFinding(graph, startCity, endCity, cities);
            stop = high_resolution_clock::now();
            algo_duration = duration<double>(stop - start).count();
            std::cout << "Time taken by A* algorithm: " << algo_duration << " seconds" << endl;
            outinfo.astar_duration = algo_duration;

            // open new file and write output information to it, overwriting if it has been seen before
            outinfo.dijkstra_path = "";
            outinfo.astar_path = "";
            for(unsigned long i = 0; i < dijksta_res.size(); i++) {
                if(i != dijksta_res.size() -1) {
                    outinfo.dijkstra_path+=dijksta_res[i];
                    outinfo.dijkstra_path+=" ";
                } else {
                    outinfo.dijkstra_path+=dijksta_res[i];
                }
            }
            for(unsigned long i = 0; i < astar_res.size(); i++) {
                if(i != dijksta_res.size() -1) {
                    outinfo.astar_path+=astar_res[i];
                    outinfo.astar_path+=" ";
                } else {
                    outinfo.astar_path+=astar_res[i];
                }
            }
            // open file and write 
            ofstream outputFile("output/" + startCity + " " + endCity + ".txt");
            if (!outputFile.is_open()) {
                cerr << "Error opening output file" << endl;
                return 0;
            }

            // Output header line
            outputFile << "Results when running " << outinfo.start_city << " to " << outinfo.end_city << endl;
            outputFile << "Dijkstra: " << outinfo.dijkstra_duration << " seconds"<< endl;
            outputFile << "Dijkstra Path: " << outinfo.dijkstra_path << endl;
            outputFile << "A*: " << outinfo.astar_duration << " seconds"<< endl;
            outputFile << "A* Path: " << outinfo.astar_path << endl;
            if(outinfo.dijkstra_path == outinfo.astar_path) {
                outputFile << "Dijkstra and Astar gave the same resulting path" << endl;
            }
            else {
                outputFile << "Dijkstra and Astar did not give the same resulting path" << endl;
            }
            outputFile.close();
            std::cout << "Results data written to " << "output/" + startCity + " " + endCity + ".txt" << endl << endl;;

        } catch (exception& e) {
            std::cout << "Error finding path: " << e.what() << endl;
        }
    }
    
    return 0;
}