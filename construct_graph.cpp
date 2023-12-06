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

struct City {
    string name;
    double lat;
    double lng;
};

string removeQuotesAndTrim(const string& input) {
    string result;
    // Remove quotes
    for (char c : input) {
        if (c != '\"') {
            result += c;
        }
    }
    // Trim spaces
    size_t first = result.find_first_not_of(' ');
    if (string::npos == first) {
        return result;
    }
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

int main() {
    ifstream file("cities.csv");
    string line;
    if (!file.is_open()) {
        cerr << "Error opening file" << endl;
        return 1;
    }

    vector<City> cities;
    getline(file, line); // Skip the first line (header)
    cout << "Creating cities vector..." << endl;

    while (getline(file, line)) {
        stringstream ss(line);
        string city, lat, lng, token;
        int column = 0;
        while (getline(ss, token, ',')) {
            column++;
            if (column == 1) {
                city = removeQuotesAndTrim(token);
            } else if (column == 7) {
                lat = removeQuotesAndTrim(token);
            } else if (column == 8) {
                lng = removeQuotesAndTrim(token);
            }
        }

        try {
            City newCity = {city, stod(lat), stod(lng)};
            cities.push_back(newCity);
        } catch (const std::invalid_argument& e) {
            // Skips any invalid entries from the CSV file
            cerr << "Invalid data for city: " << city << ". Skipping." << endl;
        }
    }

    file.close();

    // Create a graph (defined as an adjacency list)
    unordered_map<string, vector<pair<string, double>>> graph; // Maps city to a list of other cities and their distances

    cout << "Calculating distances..." << endl;
    for (const auto& city1 : cities) {
        vector<pair<string, double>> cityDistances;
        for (const auto& city2 : cities) {
            if (city1.name != city2.name) {
                double distance = calculateDistance(city1.lat, city1.lng, city2.lat, city2.lng);
                cityDistances.push_back(make_pair(city2.name, distance));
            }
        }
        graph[city1.name] = cityDistances;
    }

    return 0;
}