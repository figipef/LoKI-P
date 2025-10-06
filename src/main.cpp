#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <iomanip>

#include "configparser.h"
#include "inputparser.h"
#include "specie.h"
#include "speciesloader.h"
void inputconfirmer(InputParser cfg){
    try {

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "\n=== INPUT FILE SUMMARY ===\n";
        std::cout << "Gas Temperature: " << cfg.gasTemperature << " K\n";
        std::cout << "Gas Pressure:    " << cfg.gasPressure << " Torr\n";
        std::cout << "Electron Density: " << cfg.electronDensity << " m^-3\n";
        std::cout << "Electron Mean Energy: " << cfg.electronMeanEnergy << " eV\n";
        std::cout << "Secondary Electron Mean Energy: " << cfg.secondaryElectronMeanEnergy << " eV\n";
        std::cout << "Permitivitiy: " << cfg.permitivity << "\n";
        std::cout << "Plasma range: " << cfg.plasmaInit << " -> " << cfg.plasmaEnd << "\n";
        std::cout << "Potential: Left = " << cfg.leftPotential
                  << " V, Right = " << cfg.rightPotential << " V\n";
        std::cout << "Potential Type: " << cfg.potentialType << " (Freq = " << cfg.freq << " s^-1)\n";
        std::cout << "Chemistry: " << (cfg.chemistryOn ? "ON" : "OFF") << "\n";
        std::cout << "Number of species: " << cfg.numberSpecies << "\n";
        std::cout << "Number of reactions: " << cfg.numberReactions << "\n\n";

        std::cout << "--- Grid Configuration ---\n";
        for (size_t i = 0; i < cfg.gridSizes.size(); ++i) {
            std::cout << "  Segment " << (i+1)
                      << ": GRID=" << cfg.gridSizes[i]
                      << ", LENGTH=" << cfg.lengths[i] << " m\n";
        }
        std::cout << std::endl;

        std::cout << "--- Species ---\n";
        for (size_t i = 0; i < cfg.speciesNames.size(); ++i) {
            std::cout << "  " << cfg.speciesNames[i]
                      << " (m=" << cfg.speciesMasses[i]
                      << ", q=" << cfg.speciesCharges[i] << ")\n";
        }

        std::cout << "\n--- React Table ---\n";
        for (const auto& row : cfg.reactTable) {
            for (int val : row) std::cout << val << " ";
            std::cout << "\n";
        }

        std::cout << "\n--- Wall React Table ---\n";
        for (const auto& row : cfg.wallReactTable) {
            for (int val : row) std::cout << val << " ";
            std::cout << "\n";
        }

        std::cout << "\n=== END OF CONFIG SUMMARY ===\n\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
}

int main() {
    ConfigParser solver_cfg;
    InputParser input_cfg;

    try {
        solver_cfg = parseSolverConfig("../config/config.txt");
        std::cout << "Solver type: " << solver_cfg.solverType << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error parsing solver config: " << e.what() << "\n";
    }

    try {
        input_cfg = parseInputFile("../INPUT.txt");
        std::cout << "Parsed input file successfully.\n";
    } catch (const std::exception& e) {
        std::cerr << "Error parsing input file: " << e.what() << "\n";
        return 1;
    }

    inputconfirmer(input_cfg); // prints the  input file configuration details

    std::vector<Specie> particles = create_species_from_config(input_cfg);

    // Print particle info
    for (const auto& p : particles) {
        std::cout << "Particle ID: " << p.get_id() << "\n";
        std::cout << "Name: " << p.get_name() << "\n";
        std::cout << "Charge: " << p.get_charge() << "\n";
        std::cout << "Mass: " << p.get_mass() << "\n";
        std::cout << "q/m Ratio: " << p.get_qm_ratio() << "\n";
        std::cout << "Density Matrix:\n" << p.get_density() << "\n\n";
    }

    return 0;
}
