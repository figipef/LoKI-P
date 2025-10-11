#pragma once

#include <vector>

#include "grid.h"

class Poisson {
public:
    Poisson(const Grid& grid,
            const std::vector<double>& eps,
            const bool& is_axial);

    // --- Accessors ---
    const std::vector<double>& phic() const { return phic_; }
    const std::vector<double>& phiw() const { return phiw_; }
    const std::vector<double>& phie() const { return phie_; }

    const std::vector<double>& boundary_conditions() const { return boundary_conditions_; }
    const std::vector<double>& surface_charges() const { return surface_charges_; }

    void update_boundary_values();
    void update_surface_charges(double left_charge, double right_charge);

    void solve();

    void print_summary() const;

private:
    const Grid& grid_;
    std::vector<double> eps_;

    // --- Coefficients ---
    std::vector<double> phic_;  // central
    std::vector<double> phiw_;  // west
    std::vector<double> phie_;  // east

    // --- Physical boundaries ---
    std::vector<double> boundary_conditions_; // [left, right] — in volts
    std::vector<double> surface_charges_;     // per boundary, default 0 (used later by solver)

    bool is_axial_; // Geometry definition for equation computations

    std::vector<double> rhs_surface_charges;
    std::vector<double> rhs_boundary_potential;

    // vectors to be used for thomas algorithm
    std::vector<double> a; // lower diagonal
    std::vector<double> b; // main diagonal
    std::vector<double> c; // upper diagonal

    void initialize_coefficients();

    void define_lhs_vectors();
};
