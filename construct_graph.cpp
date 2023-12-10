#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <unordered_map>

#define EARTH_RADIUS_KM 6371.0      /* Radius of the Earth in Kilometers */
#define DISTANCE_THRESHOLD_KM 500.0 /* Distance Threshold */
#define CITY_LIMIT 50

using namespace std;

struct City {
    string name;
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
        cout << entry.first << ": ";
        for (const auto& neighbor : entry.second) {
            cout << "(" << neighbor.first << ", " << neighbor.second << " km) ";
        }
        cout << endl;
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
            // Separate fields with commas
            outputFile << cityPair.first << "," << connection.first << "," << connection.second << endl;
        }
    }

    outputFile.close();
    cout << "Graph data written to " << filename << endl;
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
    cout << "Creating cities vector..." << endl;

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
            City newCity = {city, stod(lat), stod(lng)};
            cities.push_back(newCity);
        } catch (const std::invalid_argument& e) {
            cerr << "Invalid data for city: " << city << ". Skipping." << endl;
        }

        // Use only first CITY_LIMIT cities in CSV for testing
        if (cities.size() >= CITY_LIMIT) {
            break;
        }
    }

    file.close();
    return cities;
}

unordered_map<string, vector<pair<string, double>>> createCityGraph(const vector<City>& cities) {
    unordered_map<string, vector<pair<string, double>>> graph;

    cout << "Calculating distances..." << endl;
    for (const auto& city1 : cities) {
        vector<pair<string, double>> cityDistances;
        for (const auto& city2 : cities) {
            if (city1.name != city2.name) {
                double distance = calculateDistance(city1.lat, city1.lng, city2.lat, city2.lng);
                if (distance <= DISTANCE_THRESHOLD_KM) {
                    cityDistances.push_back(make_pair(city2.name, distance));
                }
            }
        }
        graph[city1.name] = cityDistances;
    }

    return graph;
}

unordered_map<string,string> dijkstra(unordered_map<string, vector<pair<string, double>>> graph, string start_city)
{
    unordered_map<string,int> dist_vec; /* Keep track of distances */
 
    unordered_map<string,bool> visited; /* Keep track of visited nodes */

    unordered_map<string,string> predecessor_list;

    /* INITIALIZE VALUES */
    for(auto i = graph.begin(); i != graph.end(); i++) {
        dist_vec[i->first] = INT_MAX;
        visited[i->first] = false;
        predecessor_list[i->first] = "";
    }

    dist_vec[start_city] = 0;
    predecessor_list[start_city] = start_city;

    for (int iter = 0; iter < graph.size() - 1; iter++) {
        /* USE INFO AT HAND TO SELECT MINIMUM WEIGHT ADDITION EDGE */
        int curr_min = INT_MAX;
        string min_idx = "";
        for (auto i = graph.begin(); i != graph.end(); i++) {
            if(!visited[i->first] && dist_vec[i->first] <= curr_min) {
                curr_min = dist_vec[i->first];
                min_idx = i->first;
            }
        }
        visited[min_idx] = true;
        /* ITERATE THROUGH CONNECTED NODES */
        for(auto i = graph.begin(); i != graph.end(); i++) {
            string city_name = i->first;
            auto city_pair_vec = graph[city_name];
            double dist_to_main = INT_MAX;
            bool is_connected = false;
            for(auto i = city_pair_vec.begin(); i < city_pair_vec.end(); i++) {
                if(min_idx == i->first) {
                    dist_to_main = i->second;
                    is_connected = true;
                }
            }
            if(!visited[city_name] && is_connected && dist_vec[min_idx] + dist_to_main < dist_vec[city_name]) {
                dist_vec[city_name] = dist_vec[min_idx] + dist_to_main;
                predecessor_list[city_name] = min_idx;
            }
        }
    }
    return predecessor_list;
}

void printPathVector(vector<string> path_vec) {
    for(int i = 0; i < path_vec.size(); i++) {
        if(i == path_vec.size() - 1) {
            cout << path_vec[i] << endl;
            break;
        } else
            cout << path_vec[i] << ", ";
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
    printPathVector(path_vec);
    return path_vec;
}

int main() {
    try {
        vector<City> cities = readCitiesFromFile("cities.csv");
        auto graph = createCityGraph(cities);
        // printGraph(graph);
        outputGraphToFile(graph, "small_adjacency_list.csv");
        dijkstraPathFinding(graph, "Los Angeles", "Phoenix");
    } catch (const std::exception& e) {
        cerr << "Exception caught in main: " << e.what() << endl;
        return 1;
    }

    return 0;
}
