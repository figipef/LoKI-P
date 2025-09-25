#include "inputparser.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

static void trim(std::string& s) { 
    size_t start = s.find_first_not_of(" \t");
    size_t end   = s.find_last_not_of(" \t\r\n");
    if (start == std::string::npos) {
        s.clear();
    } else {
        s = s.substr(start, end - start + 1);
    }
}

InputParser parseInputFile(const std::string& filename) {
    InputParser cfg;
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Cannot open input file: " + filename);

    std::string line;
    bool readingReactTable = false;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        if (readingReactTable) {
            if (line.find_first_not_of(" \t\r\n") == std::string::npos) {
                readingReactTable = false;
                continue;
            }
            std::istringstream iss(line);
            std::vector<int> row;
            int val;
            while (iss >> val) row.push_back(val);
            if (!row.empty()) cfg.reactTable.push_back(row);
            continue;
        }

        auto pos = line.find('=');
        if (pos == std::string::npos) continue;

        std::string key   = line.substr(0, pos);
        std::string value = line.substr(pos + 1);
        trim(key);
        trim(value);

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
        else if (key == "POTENTIAL_TYPE") cfg.potentialType = value;
        else if (key == "FREQ") cfg.freq = std::stod(value);
        else if (key == "LFA/LEA") cfg.lfaOrLea = value;
        else if (key == "CHEMISTRY") cfg.chemistryOn = (value.find("ON") != std::string::npos);
        else if (key == "NUMBER_SPECIES") cfg.numberSpecies = std::stoi(value);
        else if (key == "NUMBER_REACTIONS") cfg.numberReactions = std::stoi(value);
        else if (key == "REACT_TABLE") {
            readingReactTable = true;

        }
    }

    return cfg;
}
