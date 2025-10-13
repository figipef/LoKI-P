#pragma once

#include <vector>

#include "grid.h"
#include "specie.h"

class Poisson {
public:
    Poisson(const Grid& grid,
            const std::vector<double>& eps,
            const bool& is_axial);

    // --- Accessors ---
    const std::vector<double>& phic() const { return phic_; }
    const std::vector<double>& phiw() const { return phiw_; }
    const std::vector<double>& phie() const { return phie_; }

    const std::vector<int>& boundary_conditions() const { return boundary_conditions_; }
    const std::vector<double>& surface_charges() const { return surface_charges_; }

    void update_boundary_values(const std::vector<int>& boundary_condition, const std::vector<double>& surface_charges);
    void update_surface_charges(const std::vector<int>& boundary_idxs, const std::vector<double>& surface_charges);

    std::vector<double> solve_thomasalg(const std::vector<double>& charge_density);

    void solve();

    void print_summary() const;

//private:
    const Grid& grid_;
    std::vector<double> eps_;

    // --- Coefficients --- | Thomas Algorithm vectors
    std::vector<double> phic_;  // central | central diagonal
    std::vector<double> phiw_;  // west    | lower diagonal
    std::vector<double> phie_;  // east    | upper diagonal

    // --- Physical boundaries ---
    std::vector<int> boundary_conditions_; // [left, right] — 0 Dirichlet | 1 Neumann
    std::vector<double> boundary_values_;     // [left, right] — Voltage values in V
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

    void calculate_rhs_charges(int boundary_idx, double surface_charge, const std::vector<double>& centers, const std::vector<double>& boundaries, const std::vector<double>& permitivity);
};
