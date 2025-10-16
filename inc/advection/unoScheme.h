// UNOScheme.h
#pragma once
#include "advection/AdvectionScheme.h"
#include "advection/UnoCalculator.h"
#include <memory>
#include <vector>

class UNOScheme : public AdvectionScheme {
public:
    UNOScheme(size_t plasma_start, size_t plasma_end, std::shared_ptr<UnoCalculator> gc)
        : AdvectionScheme(plasma_start, plasma_end), gc_calculator_(std::move(gc)) {}

    std::vector<double> compute_flux(
        const std::vector<double>& density,
        const std::vector<double>& velocity,
        const std::vector<double>& centers,
        const std::vector<double>& sizes,
        const double& dt) const override;

private:
    std::shared_ptr<UnoCalculator> gc_calculator_;
};
