#pragma once

#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>

class random_t
{
public:
    random_t(std::vector<int> seeds);
    random_t(const random_t&) = default;
    void use_seed(std::vector<int> seed);
    
    int var(int from, int to) const;
    int binary() const;
    int func(int from, int to) const;
    int oper(int from, int to) const;
    bool choose() const;
    double prob() const;
    
    double value(double m1, double m2, double stddev = 1.0) const;
    int integer(int from, int to) const;
    //operator std::string() const;

    template <typename T>
    std::vector<T> sample(const std::vector<T>& population,
                          std::size_t n_samples);
private:
    std::shared_ptr<std::seed_seq> _seed_seq;
    static std::mt19937 _gen;
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
