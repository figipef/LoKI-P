/*
Configuration File Parser
*/

#include "configparser.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

void ConfigParser::load(const std::string& filename) {

    std::ifstream file("../" + filename);

    if (!file) throw std::runtime_error("Failed to open config file");

    std::string line;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            data[trim(key)] = trim(value);
        }
    }
}

std::string ConfigParser::get(const std::string& key) const {
    auto it = data.find(key);
    if (it == data.end())
        throw std::runtime_error("Missing config key: " + key);
    return it->second;
}

std::string ConfigParser::trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\r\n");
    size_t end = str.find_last_not_of(" \t\r\n");
    return (start == std::string::npos) ? "" : str.substr(start, end - start + 1);
}