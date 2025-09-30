#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "configparser.h"
#include "inputparser.h"
#include "specie.h"
#include "speciesloader.h"

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
