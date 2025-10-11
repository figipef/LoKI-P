#include "poisson.h"
#include <iostream>
#include <stdexcept>
#include <cmath>

Poisson::Poisson(const Grid& grid,
                 const std::vector<double>& eps)
    : grid_(grid), eps_(eps)
{
    int n = grid_.size();
    if (eps_.size() != static_cast<size_t>(n)) {
        throw std::runtime_error("Permittivity vector size must match grid size.");
    }

    phiw_.resize(n, 0.0);
    phie_.resize(n, 0.0);
    phic_.resize(n, 0.0);

    boundary_conditions_.resize(2, 0.0); // default: both 0 V
    surface_charges_.resize(n + 1, 0.0); // one per boundary

    initialize_coefficients();

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

        phiw_[i] =  (i > 0) ? phie_[i-1] : 0;
        phie_[i] =  (i + 1 < n) ? east_coef: 0;
        phic_[i] = -(phiw_[i] + phie_[i]);
    }
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
