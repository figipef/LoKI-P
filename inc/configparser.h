#pragma once
#include <string>

// Struct to hold solver configuration values
struct ConfigParser {
    std::string solver;           // "Explicit", "Implicit", etc.
    std::string solverType;       // "RK4", "Euler", etc.
    std::string convectionScheme; // "UNO2", "Upwind", etc.
    bool chemistry = false;       // true if "Yes"/"ON", false otherwise
    bool is_axial = true;         // true if geometry is axial, false if radial
};

// Parse a config file and fill a SolverConfig struct
ConfigParser parseSolverConfig(const std::string& filename);
