#pragma once

#include <string>
#include <vector>
#include <utility>

#include "common_types.hpp"

std::vector<entry_t> load_data(std::string filename);
std::vector<std::vector<int>> load_seeds(std::string filename);
