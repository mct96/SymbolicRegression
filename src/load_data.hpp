#pragma once

#include <string>
#include <vector>
#include <utility>

namespace sr
{

using data_t = std::vector<double>;
using Xy = std::pair<data_t, double>;

std::vector<Xy> load_data(std::string filename);
std::vector<std::vector<int>> load_seeds(std::string filename);

}
