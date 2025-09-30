#include "speciesloader.h"
#include <stdexcept>

std::vector<Specie> create_species_from_config(const InputParser& cfg) {
    std::vector<Specie> species;

    // Validate species data
    if (cfg.speciesNames.size() != cfg.numberSpecies ||
        cfg.speciesMasses.size() != cfg.numberSpecies ||
        cfg.speciesCharges.size() != cfg.numberSpecies ||
        //cfg.reactTable.size() != cfg.numberSpecies ||
        (!cfg.reactTable.empty() && cfg.reactTable[0].size() != cfg.numberReactions)) {
        throw std::runtime_error("Mismatch in species data dimensions.");
    }

    for (int i = 0; i < cfg.numberSpecies; ++i) {
        const std::string& name = cfg.speciesNames[i];
        int charge = cfg.speciesCharges[i];
        double mass = cfg.speciesMasses[i];
        const std::vector<int>& react_net = cfg.reactTable[i];

        Eigen::MatrixXd n(cfg.gridSize, 1);
        n.setZero(); // Initialize with zeros

        species.emplace_back(i, name, charge, mass, react_net, n);
    }

    // Add electron energy specie if LEA mode is active
    if (cfg.lfaOrLea == "LEA") {
        if (cfg.reactTable.empty()) {
            throw std::runtime_error("Missing reaction table for ELECTRON_ENERGY.");
        }

        Eigen::MatrixXd n(cfg.gridSize, 1);
        n.setZero(); // Initialize with zeros

        species.emplace_back(-1, "ELECTRON_ENERGY", -1, 1.0, cfg.reactTable[0], n);
    }

    return species;
}
