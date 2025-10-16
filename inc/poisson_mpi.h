#pragma once

#include <mpi.h>
#include <vector>
#include "grid.h"  // Your Grid class

class PoissonMPI {
public:
    PoissonMPI(const Grid& grid,
               const std::vector<double>& eps,
               bool is_axial = false);

    void initialize_coefficients();
    void update_boundary_values(const std::vector<int>& boundary_conditions,
                                const std::vector<double>& boundary_values);
    void update_surface_charges(const std::vector<int>& boundary_idxs,
                                const std::vector<double>& surface_charges);
    std::vector<double> solve(const std::vector<double>& charge_density);
    void print_summary() const;

private:
    void calculate_rhs_charges(int boundary_idx, double surface_charge,
                               const std::vector<double>& centers,
                               const std::vector<double>& boundaries,
                               const std::vector<double>& permitivity);
    std::vector<double> solve_thomasalg(   
    const std::vector<double>& phiw_, // phiw
    const std::vector<double>& phic_, // phic
    const std::vector<double>& phie_, // phie
    const std::vector<double>& rhs);

    std::vector<double> solve_thomasalg_global(
    const std::vector<double>& a, // phiw
    const std::vector<double>& b, // phic
    const std::vector<double>& c, // phie
    const std::vector<double>& rhs);
    std::vector<double> solve_tridiagonal_thomas(
    const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    const std::vector<double>& d_in);

    const Grid& grid_;
    std::vector<double> eps_;
    bool is_axial_;

    std::vector<double> phiw_, phie_, phic_;
    std::vector<double> rhs_surface_charges, rhs_boundary_potential;
    std::vector<int> boundary_conditions_;
    std::vector<double> boundary_values_;
    std::vector<double> surface_charges_;

    int rank_, size_;
    int local_start_, local_end_;
    int local_N_;
};
