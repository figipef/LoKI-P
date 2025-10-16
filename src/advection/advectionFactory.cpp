#include "advection/advectionFactory.h"
#include "advection/unoScheme.h"
#include "advection/unoCalculator.h"
#include "advection/uno2.h"
#include <stdexcept>

std::shared_ptr<AdvectionScheme> AdvectionFactory::create(
    const std::string& scheme_name,
    size_t plasma_start,
    size_t plasma_end,
    const std::string& limiter)
{
    if (scheme_name == "UNO") {
        std::shared_ptr<UnoCalculator> gc;
        if (limiter == "2") gc = std::make_shared<Uno2>();
        else throw std::runtime_error("Unknown UNO Gc variant: " + limiter);
        return std::make_shared<UNOScheme>(plasma_start, plasma_end, gc);
    }

    throw std::runtime_error("Unknown advection scheme: " + scheme_name);
}
