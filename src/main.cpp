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

    std::vector<Specie> species = create_species_from_config(input_cfg);

    // Print particle info
    for (const auto& s : species) {
        std::cout << "Specie ID: " << s.get_id() << "\n";
        std::cout << "Name: " << s.get_name() << "\n";
        std::cout << "Charge: " << s.get_charge() << "\n";
        std::cout << "Mass: " << s.get_mass() << "\n";
        std::cout << "q/m Ratio: " << s.get_qm_ratio() << "\n";
        std::cout << "Density Matrix:\n" << s.get_density() << "\n\n";
    }

    return 0;
}
