#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <initializer_list>

class statistics_t
{
    friend class summarize_statistics_t;
public:
    statistics_t() = default;
    statistics_t(const statistics_t&) = default;
    statistics_t(statistics_t&&) = default;
    statistics_t(std::vector<double> samples);
    statistics_t(std::initializer_list<double> samples);
    ~statistics_t() = default;
    statistics_t& operator=(const statistics_t&) = default;
    statistics_t& operator=(statistics_t&&) = default;
    
    void update_values(std::vector<double> samples);

    double mean() const;
    double median() const;
    double max() const;
    double min() const;
    double var() const;
    double stddev() const;
    double total() const;

    operator std::string() const;

private:
    double _mean;
    double _median;
    double _max;
    double _min;
    double _var;
    double _stddev;
    double _total;
};

class summarize_statistics_t
{
    using Tp = std::pair<double, double>;
public:
    summarize_statistics_t() = default;
    summarize_statistics_t(const summarize_statistics_t&) = default;
    summarize_statistics_t(summarize_statistics_t&&) = default;
    summarize_statistics_t(const std::vector<statistics_t>& samples);
    summarize_statistics_t(std::initializer_list<statistics_t> samples);
    ~summarize_statistics_t() = default;
    summarize_statistics_t& operator=(const summarize_statistics_t&) = default;
    summarize_statistics_t& operator=(summarize_statistics_t&&) = default;
    
    void update_values(const std::vector<statistics_t>& samples);

    Tp min() const;
    Tp max() const;
    Tp mean() const;
    Tp median() const;
    Tp var() const;
    Tp stddev() const;

private:
    double mean(const std::vector<statistics_t>& samples,
                double(*)(const statistics_t&)) const;

    double stddev(const std::vector<statistics_t>& samples,
                  double mean,
                  double(*)(const statistics_t&)) const;

    Tp _min; // mean and standard deviation
    Tp _max;
    Tp _mean;
    Tp _median;
    Tp _var;
    Tp _stddev;
};
