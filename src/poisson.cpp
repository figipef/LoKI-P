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

    boundary_conditions_.resize(2, 0); // default: both Dirichlet
    boundary_values_.resize(2, 0); // default: both 0 V
    surface_charges_.resize(n + 1, 0.0); // one per boundary

    initialize_coefficients();// Define the vectors for Thomas Algorithm

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

    // Update the first and last values to account for ghost cells
    phiw_[0] = phiw_[1];
    phie_[n-1] = phie_[n-2];
    phic_[0] -= phiw_[0];
    phic_[n-1] -= phie_[n-1];
}

void Poisson::update_boundary_values(const std::vector<int>& boundary_conditions, const std::vector<double>& boundary_values)
{

    boundary_conditions_ = boundary_conditions;
    boundary_values_ = boundary_values;

    int n = grid_.size();

    if (boundary_conditions[0]){ // Neumann left side
        phic_[0] += phiw_[0];
    } else { // Dirichlet left side
        rhs_boundary_potential[0] = -boundary_values[0] * phiw_[0];
    }

    if (boundary_conditions[1]){ // Neumann right side
        phic_[n-1] += phie_[n-1];
    } else { // Dirichelet right side
        rhs_boundary_potential[n-1] = -boundary_values[1] * phie_[n-1];
    }
}

void Poisson::update_surface_charges(const std::vector<int>& boundary_idxs, const std::vector<double>& surface_charges)
{
    // Function to update the surface charges in 
    // the right hand side in the poisson equation

    const std::vector<double>& centers = grid_.cell_centers();
    const std::vector<double>& boundaries = grid_.boundaries();

    const std::vector<double> permitivity = grid_.permitivity();

    int N = boundary_idxs.size();

    for (int i = 0; i < N; i++ ){
        surface_charges_[boundary_idxs[i]] = surface_charges[i];
        rhs_surface_charges[boundary_idxs[i]] = 0.0;
    }

    for (int i = 0; i < N; i++ ){
        // Update the right hand side values
        calculate_rhs_charges(boundary_idxs[i], surface_charges_[i], centers, boundaries, permitivity);
    }
}

void Poisson::calculate_rhs_charges(int boundary_idx, double surface_charge, const std::vector<double>& centers, const std::vector<double>& boundaries, const std::vector<double>& permitivity)
{

    // STILL NEEDS TO BE MULTIPLIED BY THE GEOMETRY FACTOR IN CASE OF RADIAL CALCULATIONS
 
    // boundary at i, is between cell i and i - 1
    double nume = centers[boundary_idx] - boundaries[boundary_idx]; // numerator x(i+1) - x(i+1/2)
    // denominator eps(i) * (x(i+1) - x(i+1/2)) + eps(i+1) * (x(i+1/2) - x(i))
    double denom = permitivity[boundary_idx - 1] * (centers[boundary_idx] - boundaries[boundary_idx]) +
                   permitivity[boundary_idx] * (boundaries[boundary_idx] - centers[boundary_idx - 1]); 

    // STILL NEEDS TO BE MULTIPLIED BY THE GEOMETRY FACTOR IN CASE OF RADIAL CALCULATIONS    
    rhs_surface_charges[boundary_idx] += surface_charge * permitivity[boundary_idx] * (-1) * nume / denom;
    rhs_surface_charges[boundary_idx - 1] += surface_charge * permitivity[boundary_idx-1] * nume / denom;
}

std::vector<double> Poisson::solve_thomasalg(const std::vector<double>& charge_density)
{
    const int N = charge_density.size();
    const std::vector<double>& boundaries = grid_.boundaries();

    std::vector<double> rhs;
    rhs.resize(N,0.0);

    for (int i = 0; i < N; i++){

        double geom_factor; // Geometry factor for a axial

        if (!is_axial_){
            geom_factor = (boundaries[i+1] * boundaries[i+1] - boundaries[i]* boundaries[i]) * M_PI;
        } else {
            geom_factor = (boundaries[i+1] - boundaries[i]);
        }

        rhs[i] = geom_factor * charge_density[i] + rhs_surface_charges[i] + rhs_boundary_potential[i];
    }

    // temporary vectors for thomas algorithm
    std::vector<double> v1(N-1, 0.0);
    std::vector<double> v2(N, 0.0);
    std::vector<double> v3(N, 0.0);

    v1[0] = phie_[0]/phic_[0];
    v2[0] = rhs[0]/phic_[0];

    for (int i = 1; i < N-1; i++){
        v1[i] = phie_[i]/(phic_[i] - phiw_[i-1]*v1[i-1]);
    }
    for (int i = 1; i < N; i++){
        v2[i] = (rhs[i] - phiw_[i-1] *v2[i-1])/(phic_[i] - phiw_[i-1]*v1[i-1]);
    }

    v3[N-1] = v2[N-1];

    for (int i = N-1; i > 0; i--){
        v3[i-1] = v2[i-1] - v1[i-1]*v3[i];
    }

    return v3; // Potential Vector
}

void Poisson::print_summary() const
{
    std::cout << "\n=== Poisson Operator Summary ===\n";
    std::cout << "Grid size: " << grid_.size() << "\n";
    std::cout << "Plasma start index: " << grid_.plasma_start_index() << "\n";
    std::cout << "Plasma end index: " << grid_.plasma_end_index() << "\n";
    std::cout << "Total length: " << grid_.total_length() << " m\n";

    std::cout << "\nBoundary conditions (V): "
              << "Left = " << boundary_values_[0]
              << ", Right = " << boundary_values_[1] << "\n";

    std::cout << "phic (central): ";
    for (auto v : phic_) std::cout << v << " ";
    std::cout << "\nphiw (west):   ";
    for (auto v : phiw_) std::cout << v << " ";
    std::cout << "\nphie (east):   ";
    for (auto v : phie_) std::cout << v << " ";
    std::cout << "\n=================================\n";
}
