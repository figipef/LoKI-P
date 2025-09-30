#include "inputparser.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <algorithm>

static std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t");
    size_t end = s.find_last_not_of(" \t\r\n");
    return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

static std::string strip_units(const std::string& value) {
    return value.substr(0, value.find_first_of(" \t"));
}

static std::vector<std::string> parseQuotedList(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ',')) {
        token = trim(token);
        if (!token.empty() && token.front() == '"' && token.back() == '"') {
            result.push_back(token.substr(1, token.size() - 2));
        } else {
            result.push_back(token);
        }
    }
    return result;
}

static std::vector<int> parseIntList(const std::string& line) {
    std::vector<int> result;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ',')) {
        result.push_back(std::stoi(trim(token)));
    }
    return result;
}

static std::vector<double> parseDoubleList(const std::string& line) {
    std::vector<double> result;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ',')) {
        result.push_back(std::stod(strip_units(trim(token))));
    }
    return result;
}

InputParser parseInputFile(const std::string& filename) {
    InputParser cfg;
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Cannot open input file: " + filename);

    std::string line;
    bool readingReactTable = false;
    bool readingWallReactTable = false;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::string trimmed = trim(line);

        // Check for block transitions
        if (readingReactTable && trimmed.find("WALL_REACT_TABLE") != std::string::npos) {
            readingReactTable = false;
            readingWallReactTable = true;
            continue;
        }

        if (readingReactTable) {
            if (trimmed.empty()) {
                readingReactTable = false;
                continue;
            }
            std::istringstream iss(trimmed);
            std::vector<int> row;
            int val;
            while (iss >> val) row.push_back(val);
            if (!row.empty()) cfg.reactTable.push_back(row);
            continue;
        }

        if (readingWallReactTable) {
            if (trimmed.empty()) {
                readingWallReactTable = false;
                continue;
            }
            std::istringstream iss(trimmed);
            std::vector<int> row;
            int val;
            while (iss >> val) row.push_back(val);
            if (!row.empty()) cfg.wallReactTable.push_back(row);
            continue;
        }

        auto pos = trimmed.find('=');
        if (pos == std::string::npos) continue;

        std::string key = trim(trimmed.substr(0, pos));
        std::string raw_value = trim(trimmed.substr(pos + 1));
        std::string value = strip_units(raw_value);

        try {
            if (key == "GAS_TEMPERATURE") cfg.gasTemperature = std::stod(value);
            else if (key == "GAS_PRESSURE") cfg.gasPressure = std::stod(value);
            else if (key == "ELECTRON_DENSITY") cfg.electronDensity = std::stod(value);
            else if (key == "ELECTRON_MEAN_ENERGY") cfg.electronMeanEnergy = std::stod(value);
            else if (key == "SECONDARY_ELECTRON_MEAN_ENERGY") cfg.secondaryElectronMeanEnergy = std::stod(value);
            else if (key == "GRID_SIZE") cfg.gridSize = std::stoi(value);
            else if (key == "GRID_INIT") cfg.gridInit = std::stoi(value);
            else if (key == "GRID_END") cfg.gridEnd = std::stoi(value);
            else if (key == "PERMITIVITY") cfg.permitivity = std::stod(value);
            else if (key == "LENGTH") cfg.length = std::stod(value);
            else if (key == "LEFT_POTENTIAL") cfg.leftPotential = std::stod(value);
            else if (key == "RIGHT_POTENTIAL") cfg.rightPotential = std::stod(value);
            else if (key == "POTENTIAL_TYPE") cfg.potentialType = raw_value;
            else if (key == "FREQ") cfg.freq = std::stod(value);
            else if (key == "LFA/LEA") cfg.lfaOrLea = raw_value;
            else if (key == "CHEMISTRY") cfg.chemistryOn = (value.find("ON") != std::string::npos);
            else if (key == "NUMBER_SPECIES") cfg.numberSpecies = std::stoi(value);
            else if (key == "NUMBER_REACTIONS") cfg.numberReactions = std::stoi(value);
            else if (key == "REACT_TABLE") readingReactTable = true;
            else if (key == "WALL_REACT_TABLE") readingWallReactTable = true;
            else if (key == "NAMES") cfg.speciesNames = parseQuotedList(raw_value);
            else if (key == "MASS") cfg.speciesMasses = parseDoubleList(raw_value);
            else if (key == "CHARGE") cfg.speciesCharges = parseIntList(raw_value);
            else if (key == "SPECIES") continue; // marker only
            else {
                std::cerr << "Warning: Unknown key \"" << key << "\" in config file.\n";
            }
        } catch (const std::exception& e) {
            std::cerr << "Error parsing key \"" << key << "\": " << e.what() << "\n";
        }
    }

    return cfg;
}
