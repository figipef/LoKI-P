#pragma once
#include <vector>

class UnoCalculator {
public:
    virtual ~UnoCalculator() = default;
    virtual double compute(double G_dc, double G_cu, double x_d, double x_u, double dx_d, double u, double dt) const = 0;
};
