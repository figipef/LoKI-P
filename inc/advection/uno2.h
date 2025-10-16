#pragma once
#include "advection/unoCalculator.h"
#include <algorithm>
#include <cmath>

class Uno2 : public UnoCalculator {
public:
    double compute(double G_dc, double G_cu, double x_d, double x_u, double dx_d, double u, double dt) const override {
        // Calculate UNO2 G_c
        return std::copysign(1.0, G_dc) * std::min(std::abs(G_dc), std::abs(G_cu));
    }
};