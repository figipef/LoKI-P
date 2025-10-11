#define _USE_MATH_DEFINES

#include "poisson.h"
#include <iostream>
#include <stdexcept>
#include <cmath>

Poisson::Poisson(const Grid& grid,
                 const std::vector<double>& eps,
                 const bool& is_axial)
    : grid_(grid), eps_(eps), is_axial_(is_axial)
{
    int n = grid_.size();
    if (eps_.size() != static_cast<size_t>(n)) {
        throw std::runtime_error("Permittivity vector size must match grid size.");
    }

    phiw_.resize(n, 0.0);
    phie_.resize(n, 0.0);
    phic_.resize(n, 0.0);

    rhs_surface_charges.resize(n, 0.0);
    rhs_boundary_potential.resize(n, 0.0);

    boundary_conditions_.resize(2, 0.0); // default: both 0 V
    surface_charges_.resize(n + 1, 0.0); // one per boundary

    initialize_coefficients();// Define the vectors for Thomas Algorithm

    define_lhs_matrix(); 

    std::cout << "Initialized Poisson operator.\n";
}

void Poisson::initialize_coefficients()
{
    // Calculating the coeffecients to solve the poisson equation
    // Done a single time if permitivities don't change
    const std::vector<double>& centers = grid_.cell_centers();
    const std::vector<double>& boundaries = grid_.boundaries();
    int n = grid_.size();

    for (int i = 0; i < n; ++i) { 
        double eps_i = eps_[i];
        double eps_ip1 = (i + 1 < n) ? eps_[i + 1] : eps_i;

        double dx_i   = centers[i+1] - boundaries[i+1];
        double dx_ip1 = (i + 1 < n) ? (boundaries[i+1] - centers[i]) : 0;
        double east_coef = eps_i * eps_ip1 / (eps_ip1*dx_ip1 + eps_i*dx_i);

        if (!is_axial_){
            east_coef *= 2 * M_PI * boundaries[i+1];
        }

        phiw_[i] =  (i > 0) ? phie_[i-1] : 0;
        phie_[i] =  (i + 1 < n) ? east_coef: 0;
        phic_[i] = -(phiw_[i] + phie_[i]);
    }
}

void update_surface_charges(double left_charge, double right_charge)
{
    // Function to update the surface charges in 
    // the right hand side in the poisson equation

    const int lboundary_idx = grid_.plasma_start_index() //left boundary index for the dieletric surface charges
    const int rboundary_idx = grid_.plasma_end_index() //left boundary index for the dieletric surface charges

    const std::vector<double>& centers = grid_.cell_centers();
    const std::vector<double>& boundaries = grid_.boundaries();

    const std::vector<double> permitivity = grid_.permitivity();

    surface_charges_[lboundary_idx] = left_charge;
    surface_charges_[rboundary_idx+1] = right_charge;

    // Calculate the left and right geometry factors PI * delta_r^2
    double geom_factor_l = 1; 
    double geom_factor_r = 1;
    /*
    NEEDS TO BE ADDED FOR RADIAL CALCULATION!!!
    if (!is_axial_){
        geom_factor = M_PI * ();
    }
    */

    // STILL NEEDS TO BE MULTIPLIED BY THE GEOMETRY FACTOR IN CASE OF RADIAL CALCULATIONS
    // left case i + 1 = left boundary index
    double l_nume = centers[lboundary_idx] - boundaries[lboundary_idx] // numerator x(i+1) - x(i+1/2)
    // denominator eps(i) * (x(i+1) - x(i+1/2)) + eps(i+1) * (x(i+1/2) - x(i))
    double l_denom = permitivity[lboundary_idx - 1] * (centers[lboundary_idx] - boundaries[lboundary_idx]) + \\
                     permitivity[lboundary_idx] * (boundaries[lboundary_idx] - centers[lboundary_idx - 1]) 

    // right case i = right boundary index
    double r_nume = centers[rboundary_idx + 1] - boundaries[rboundary_idx + 1] // numerator x(i+1) - x(i+1/2)
    // denominator eps(i) * (x(i+1) - x(i+1/2)) + eps(i+1) * (x(i+1/2) - x(i))
    double r_denom = permitivity[rboundary_idx] * (centers[rboundary_idx + 1] - boundaries[rboundary_idx + 1]) + \\
                     permitivity[rboundary_idx + 1] * (boundaries[rboundary_idx + 1] - centers[rboundary_idx]) 


    rhs_surface_charges[lboundary_idx] = surface_charges_[lboundary_idx] * permitivity[lboundary_idx] * (-1) * l_nume / l_denom;
    rhs_surface_charges[lboundary_idx - 1] = surface_charges_[lboundary_idx] * permitivity[lboundary_idx-1] * l_nume / l_denom;

    rhs_surface_charges[rboundary_idx] = surface_charges_[rboundary_idx] * permitivity[rboundary_idx] * r_nume / r_denom;
    rhs_surface_charges[rboundary_idx + 1] = surface_charges_[rboundary_idx] * permitivity[rboundary_idx + 1] * (-1) *r_nume / r_denom;
}

void Poisson::set_boundary_condition(int side, double value)
{
    if (side < 0 || side > 1)
        throw std::invalid_argument("Boundary side must be 0 (left) or 1 (right).");
    boundary_conditions_[side] = value;
}

void Poisson::print_summary() const
{
    std::cout << "\n=== Poisson Operator Summary ===\n";
    std::cout << "Grid size: " << grid_.size() << "\n";
    std::cout << "Plasma start index: " << grid_.plasma_start_index() << "\n";
    std::cout << "Plasma end index: " << grid_.plasma_end_index() << "\n";
    std::cout << "Total length: " << grid_.total_length() << " m\n";

    std::cout << "\nBoundary conditions (V): "
              << "Left = " << boundary_conditions_[0]
              << ", Right = " << boundary_conditions_[1] << "\n";

    std::cout << "phic (central): ";
    for (auto v : phic_) std::cout << v << " ";
    std::cout << "\nphiw (west):   ";
    for (auto v : phiw_) std::cout << v << " ";
    std::cout << "\nphie (east):   ";
    for (auto v : phie_) std::cout << v << " ";
    std::cout << "\n=================================\n";
}
