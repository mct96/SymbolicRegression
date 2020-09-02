#pragma once

#include <algorithm>
#include <random>
#include <utility>
#include <iterator>
#include <string>

constexpr int seeds_v[] = {7, 77, 777, 7777, 77777};
std::seed_seq seeds(seeds_v, seeds_v + 5);
std::mt19937 gen{seeds};

template <typename T>
std::vector<T> sample(const std::vector<T>& population,
                      std::size_t n_individuals)
{
    auto b = population.begin(), e = population.end();
    std::vector<T> output{};

    std::sample(b, e, std::back_inserter(output), n_individuals, gen);

    return output;
}

int rd_integer(int N)
{
    std::uniform_int_distribution<> dist{0, N-1};
    return dist(gen);
}

double rd_value() // multimodal uniform distribution. Two peaks: -1 and 1.
{
    std::uniform_int_distribution<> _switch{-1, 1};
    double c = _switch(gen);
    std::normal_distribution<> dist{c, 1};
    return dist(gen);
}

double rd_real()
{
    std::uniform_real_distribution<> dist{0, 1};
    return dist(gen);
}

inline int rd_function(int N)
{
    return rd_integer(N);
}

inline int rd_operator(int N)
{
    return rd_integer(N);
}

inline int rd_variable(int N)
{
    return rd_integer(N);
}

inline int rd_binary()
{
    return rd_integer(2);
}




