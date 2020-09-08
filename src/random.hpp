#pragma once

#include <random>
#include <string>
#include <vector>
#include <algorithm>

class random_t
{
public:
    random_t(std::vector<int> seeds);
        
    int var(int from, int to) const;
    int binary() const;
    int func(int from, int to) const;
    int oper(int from, int to) const;
    bool choose() const;
    double prob() const;
    
    double value(double m1, double m2, double stddev = 1.0) const;
    int integer(int from, int to) const;
    operator std::string() const;

    template <typename T>
    std::vector<T> sample(const std::vector<T>& population,
                          std::size_t n_samples);
private:
    std::vector<int> _seeds_v;
    std::seed_seq _seeds;
    std::mt19937& _gen;
};

template <typename T>
std::vector<T> random_t::sample(const std::vector<T>& population,
                                std::size_t n_samples)
{
    auto b = population.begin(), e = population.end();
    std::vector<T> output{};

    // NOTE This function ensures stability.
    std::sample(b, e, std::back_inserter(output), n_samples, _gen);

    return output;
}
