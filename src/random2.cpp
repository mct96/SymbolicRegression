#include "random2.hpp"

#include <algorithm>
#include <utility>

random_t::random_t(std::vector<int> seeds)
    :
    _seeds_v{seeds},
    _seeds(_seeds_v.begin(), _seeds_v.end())
{
    _gen.seed(_seeds);
}

int random_t::variable(int from, int to) const
{
    std::uniform_int_distribution<> dist{from, to};
    return dist(_gen);
}

int random_t::binary() const
{
    std::uniform_int_distribution<> dist{0, 1};
    return dist(_gen);
}

int random_t::function(int from, int to) const
{
    std::uniform_int_distribution<> dist{from, to};
    return dist(_gen);
}

bool random_t::choose() const
{
    return binary() == 1;
}

double random_t::real(double m1, double m2, double stddev) const
{
    std::normal_distribution<> dist{choose() ? m1 : m2, stddev};
    return dist(_gen);
}

int random_t::integer(int from, int to) const
{
    std::uniform_int_distribution<> dist{from, to};
    return dist(_gen);
}

random_t::operator std::string() const
{
    if (_seeds_v.size() == 0) return std::string{"seeds: Є."};
    
    std::string out{"seeds: "};
    for (int i = 0; i < _seeds_v.size()-1; ++i)
        out += std::to_string(_seeds_v[i]) + ", ";
    out += std::to_string(_seeds_v.back()) + ".";

    return out;
}

template <typename T>
std::vector<T> random_t::sample(const std::vector<T>& population,
                                std::size_t n_samples)
{
    auto b = population.begin(), e = population.end();
    std::vector<T> output{};

    std::sample(b, e, std::back_inserter(output), n_samples, _gen);

    return output;
}
