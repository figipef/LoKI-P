#include "advection/unoScheme.h"
#include <cmath> 

std::vector<double> UNOScheme::compute_flux(
    const std::vector<double>& density,
    const std::vector<double>& velocity,
    const std::vector<double>& centers,
    const std::vector<double>& sizes,
    const double& dt) const     
{
    std::vector<double> fluxes(density.size() + 1, 0.0); // The number of boundaries we have in our system

    for (size_t i = plasma_start_ + 1; i < plasma_end_; ++i) {

        double flux;

        double u;
        double c;
        double d;

        double G_dc;
        double G_cu;
        double x_d;
        double x_u;
        double dx_d;
        double v = velocity[i];

        double n_c;
        double dx_c;

        if (velocity[i] > 0) // velocity at the boundary (i - 1/2) (i is defined as the)
        {
            u = i-2;
            c = i-1;
            d = i;

        } else if (velocity[i] < 0)
        {
            u = i+1;
            c = i;
            d = i-1;

        } else 
        {
            fluxes[i] = 0;
            break;
        }

        // To calculate G_c
        G_dc = (density[d] - density[c])/(centers[d] - centers[c]);
        G_cu = (density[c] - density[u])/(centers[c] - centers[u]);

        x_d = centers[d];
        x_u = centers[u];

        dx_d = sizes[d];
            
        // To calculate Flux
        n_c = density[c];
        dx_c = sizes[c];

        double Gc = gc_calculator_->compute(G_dc, G_cu, x_d, x_u, dx_d, v, dt);

        // Generic UNO-type reconstruction (example form)
        double dens_half = n_c + 0.5 * std::copysign(1.0, v) * (sizes[c] - std::abs(v)*dt) * Gc;
        fluxes[i] = dens_half * v;
    }

    return fluxes;
}