#ifndef INPUTPARSER_HPP
#define INPUTPARSER_HPP

#include <string>
#include <vector>

struct InputParser {
    double gasTemperature;
    double gasPressure;
    double electronDensity;
    double electronMeanEnergy;
    double secondaryElectronMeanEnergy;

    int gridSize;
    int plasmaInit;
    int plasmaEnd;

    double permitivity;
    std::vector<int> gridSizes;
    std::vector<double> lengths;
    std::vector<double> rperms;

    double leftPotential;
    double rightPotential;
    std::string potentialType;
    double freq;
    std::string lfaOrLea;

    bool chemistryOn;
    int numberSpecies;
    int numberReactions;

    std::vector<std::string> speciesNames;
    std::vector<double> speciesMasses;
    std::vector<int> speciesCharges;

    std::vector<std::vector<int>> reactTable;
    std::vector<std::vector<int>> wallReactTable;
};


InputParser parseInputFile(const std::string& filename);

#endif // INPUTPARSER_HPP
