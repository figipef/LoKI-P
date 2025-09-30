#ifndef INPUTPARSER_HPP
#define INPUTPARSER_HPP

#include <string>
#include <vector>

struct InputParser {
    // Physical parameters
    double gasTemperature = 0.0;
    double gasPressure = 0.0;
    double electronDensity = 0.0;
    double electronMeanEnergy = 0.0;
    double secondaryElectronMeanEnergy = 0.0;

    // Grid setup
    int gridSize = 0;
    int gridInit = 0;
    int gridEnd = 0;

    // Electrical parameters
    double permitivity = 0.0;
    double length = 0.0;
    double leftPotential = 0.0;
    double rightPotential = 0.0;
    std::string potentialType;
    double freq = 0.0;

    // Chemistry
    std::string lfaOrLea;
    bool chemistryOn = false;
    int numberSpecies = 0;
    int numberReactions = 0;

    // Species block
    std::vector<std::string> speciesNames;
    std::vector<double> speciesMasses;
    std::vector<int> speciesCharges;

    // Reaction tables
    std::vector<std::vector<int>> reactTable;
    std::vector<std::vector<int>> wallReactTable;
};

InputParser parseInputFile(const std::string& filename);

#endif // INPUTPARSER_HPP
