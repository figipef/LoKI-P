//#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <iomanip>
#include <mpi.h>

#include "configparser.h"
#include "inputparser.h"
#include "specie.h"
#include "speciesloader.h"

#include "grid.h"
#include "poisson.h"
#include "poisson_mpi.h"

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
        std::cout << "Grid Size: " << cfg.gridSize << "\n";
        std::cout << "Plasma range: " << cfg.plasmaInit << " -> " << cfg.plasmaEnd << "\n";
        std::cout << "Potential: Left = " << cfg.leftPotential
                  << " V, Right = " << cfg.rightPotential << " V\n";
        std::cout << "Potential Type: " << cfg.potentialType << " (Freq = " << cfg.freq << " s^-1)\n";
        std::cout << "Chemistry: " << (cfg.chemistryOn ? "ON" : "OFF") << "\n";
        std::cout << "Number of species: " << cfg.numberSpecies << "\n";
        std::cout << "Number of reactions: " << cfg.numberReactions << "\n\n";

        std::cout << "--- Grid Configuration ---\n";
        for (size_t i = 0; i < cfg.gridSizes.size(); ++i) {
            std::cout << "  Region " << (i+1)
                      << ": GRID=" << cfg.gridSizes[i]
                      << ", LENGTH=" << cfg.lengths[i] << " m"
                      << ", RPERM=" << ((i < cfg.rperms.size()) ? cfg.rperms[i] : cfg.permitivity)
                      << "\n";
        }

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

int main(int argc, char** argv){
    MPI_Init(&argc, &argv);  // Initialize MPI

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    ConfigParser solver_cfg;
    InputParser input_cfg;

    try {
        solver_cfg = parseSolverConfig("../config/config.txt");
        std::cout << "Solver type: " << solver_cfg.solverType << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error parsing solver config: " << e.what() << "\n";
        return 1;
    }

    try {
        input_cfg = parseInputFile("../INPUT.txt");
        std::cout << "Parsed input file successfully.\n";
    } catch (const std::exception& e) {
        std::cerr << "Error parsing input file: " << e.what() << "\n";
        return 1;
    }


    //inputconfirmer(input_cfg); // prints the  input file configuration details

    std::vector<Specie> species = create_species_from_config(input_cfg);


    // Print particle info
    //for (const auto& s : species) {
    //    std::cout << "Specie ID: " << s.get_id() << "\n";
    //    std::cout << "Name: " << s.get_name() << "\n";
    //    std::cout << "Charge: " << s.get_charge() << "\n";
    //    std::cout << "Mass: " << s.get_mass() << "\n";
    //    std::cout << "q/m Ratio: " << s.get_qm_ratio() << "\n";
    //    std::cout << "Density Matrix:\n" << s.get_density() << "\n\n";
    //}

    Grid grid(     
     input_cfg.gridSize,
     input_cfg.plasmaInit,
     input_cfg.plasmaEnd,
     input_cfg.gridSizes,
     input_cfg.lengths,
     input_cfg.rperms
    );

    //grid.print_summary();

    std::vector<double> eps(grid.size(), 1.0);
    PoissonMPI poisson(grid, eps, solver_cfg.is_axial);
    // -------------------------
    // Boundary conditions
    std::vector<int> bc = {0, 0};        // 0=Dirichlet, 1=Neumann
    std::vector<double> bv = {0.0, 1.0}; // Left=0V, Right=1V
    poisson.update_boundary_values(bc, bv);

    // -------------------------
    // Example charge density
    std::vector<double> rho(grid.size(), 0.0);

    // -------------------------
    // Solve Poisson equation
    std::vector<double> local_phi = poisson.solve(rho);
    std::vector<double> centers = grid.cell_centers();
    std::vector<double> boundaries = grid.boundaries();
    // Gather full solution on rank 0 inside solve(), or get local_phi
    std::vector<double> full_phi;
    if(rank==0){
        full_phi = local_phi; // solve() already returns full solution on rank 0
        for(int i=0;i<grid.size();i++){
            std::cout << "i="<<i<<" x="<<centers[i]<<" phi="<<full_phi[i]<<"\n";
        }
    }

    MPI_Finalize();
    poisson.print_summary();    

    return 0;
}
