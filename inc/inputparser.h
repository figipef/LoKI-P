#ifndef SIMULATION_CONFIG_HPP
#define SIMULATION_CONFIG_HPP

#include <string>
#include <vector>

struct InputParser {
    double gasTemperature = 0.0;   // K
    double gasPressure = 0.0;      // Torr
    double electronDensity = 0.0;  // m^-3
    double electronMeanEnergy = 0.0; // eV
    double secondaryElectronMeanEnergy = 0.0; // eV
    
    int gridSize = 0;
    int gridInit = 0;
    int gridEnd = 0;

    double permitivity = 0.0;
    double length = 0.0; // m

    double leftPotential = 0.0;  // V
    double rightPotential = 0.0; // V
    std::string potentialType;
    double freq = 0.0; // s^-1

    std::string lfaOrLea;
    bool chemistryOn = false;
    int numberSpecies = 0;
    int numberReactions = 0;

    std::vector<std::vector<int>> reactTable;
};

/// Parse the simulation input file and return a SimulationConfig object.
InputParser parseInputFile(const std::string& filename);

#endif // SIMULATION_CONFIG_HPP