#include <Eigen/Dense>
#include <iostream>

#include "configparser.h"
#include "inputparser.h"   

int main() { 

    ConfigParser solver_cfg;
    InputParser input_cfg;

    try {
        solver_cfg = parseSolverConfig("../config/config.txt");

    }
    catch (const std::exception& e) {
        std::cerr << "Error parsing config: " << e.what() << "\n";
    }
    
    std::cout << "compilation testing\n " << solver_cfg.solverType<<"\n"; 

    try {
        input_cfg = parseInputFile("../INPUT.txt");
    }
    catch (const std::exception& e) {
        std::cerr << "Error parsing file: " << e.what() << "\n";
    }

    return 0;
}