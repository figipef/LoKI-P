#include "configparser.h"

#include <fstream>
#include <sstream>
#include <stdexcept>

// Helper trim
static std::string trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\r\n");
    size_t end   = str.find_last_not_of(" \t\r\n");
    return (start == std::string::npos) ? "" : str.substr(start, end - start + 1);
}

ConfigParser parseSolverConfig(const std::string& filename) {
    ConfigParser cfg;

    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Failed to open config file: " + filename);

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        auto pos = line.find('=');
        if (pos == std::string::npos) continue;

        std::string key   = trim(line.substr(0, pos));
        std::string value = trim(line.substr(pos + 1));

        if (key == "Solver") cfg.solver = value;
        else if (key == "Solver_Type") cfg.solverType = value;
        else if (key == "Convection_Scheme") cfg.convectionScheme = value;
        else if (key == "Convection_Limiter") cfg.convectionLimiter = value;
        else if (key == "Chemistry") {
            if (value == "Yes" || value == "ON" || value == "True")
                cfg.chemistry = true;
            else
                cfg.chemistry = false;
        } else if (key == "Grid_Type") {
            if (value == "Axial" || value == "Ax" || value == "ax" || value == "True")
                cfg.is_axial = true;
            else if (value == "Radial" || value == "Rad" || value == "False" || value == "rad")
                cfg.is_axial = false;
            else
                throw std::runtime_error("Failed to define geometry in config file. Please select either Radial or Axial geometries");
        } 
    }

    return cfg;
}
