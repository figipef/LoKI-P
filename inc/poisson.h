#pragma once
#include <vector>
#include "grid.h"

class Poisson {
public:
    Poisson(const Grid& grid,
            const std::vector<double>& eps);

    // --- Accessors ---
    const std::vector<double>& phic() const { return phic_; }
    const std::vector<double>& phiw() const { return phiw_; }
    const std::vector<double>& phie() const { return phie_; }

    const std::vector<double>& boundary_conditions() const { return boundary_conditions_; }
    const std::vector<double>& surface_charges() const { return surface_charges_; }

    void set_boundary_condition(int side, double value);
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

    void initialize_coefficients();
};
