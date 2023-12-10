#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <unordered_map>

using namespace std;

const double EARTH_RADIUS_KM = 6371.0; // Radius of the Earth in kilometers
const double DISTANCE_THRESHOLD_KM = 500.0; // Distance threshold

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

    int edgeCount = 0;
    for (const auto& cityPair : graph) {
        edgeCount += cityPair.second.size();
    }

    outputFile << edgeCount << endl;

    for (const auto& cityPair : graph) {
        for (const auto& connection : cityPair.second) {
            outputFile << cityPair.first << " " << connection.first << " " << connection.second << endl;
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

        // Use only first 50 cities in CSV for testing
        if (cities.size() >= 50) {
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

int main() {
    try {
        vector<City> cities = readCitiesFromFile("cities.csv");
        auto graph = createCityGraph(cities);
        // printGraph(graph);
        // outputGraphToFile(graph, "small_adjacency_list.txt");
    } catch (const std::exception& e) {
        cerr << "Exception caught in main: " << e.what() << endl;
        return 1;
    }

    return 0;
}
