#pragma once
#include "specie.h"
#include "inputparser.h"
#include <vector>

std::vector<Specie> create_species_from_config(const InputParser& cfg);
