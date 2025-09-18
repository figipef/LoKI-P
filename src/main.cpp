#include <Eigen/Dense>
#include <iostream>

#include "configparser.h"

int main() { 

    ConfigParser config;
	
	std::string solverType;

	try {

        config.load("config/config.txt");
        solverType = config.get("Solver_Type");

    } catch (const std::exception& e) {
        std::cerr << "Configuration error: " << e.what() << '\n';
        return 1;
    }
    
    std::cout << "compilation testing\n " << solverType; 

    return 0;
}