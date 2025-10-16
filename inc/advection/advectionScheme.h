#pragma once
#include <vector>
#include <memory>

class AdvectionScheme {
public:
    AdvectionScheme(size_t plasma_start, size_t plasma_end)
        : plasma_start_(plasma_start), plasma_end_(plasma_end) {}

    virtual ~AdvectionScheme() = default;

    virtual std::vector<double> compute_flux(
        const std::vector<double>& density,
        const std::vector<double>& velocity,
        const std::vector<double>& centers,
        const std::vector<double>& sizes,
        const double& dt) const = 0;

protected:
    size_t plasma_start_; // index where plasma starts and ends
    size_t plasma_end_;
};