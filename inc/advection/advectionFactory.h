#pragma once
#include <memory>
#include <string>
#include "advection/advectionScheme.h"

class AdvectionFactory {
public:
    static std::shared_ptr<AdvectionScheme> create(
        const std::string& scheme_name,
        size_t plasma_start,
        size_t plasma_end,
        const std::string& limiter_name = "");
};
