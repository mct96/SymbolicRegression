#include "statistics.hpp"

statistics_t::statistics_t(std::initializer_list<double> samples)
    :
    statistics_t{std::vector<double>{samples.begin(), samples.end()}}
{
}

statistics_t::statistics_t(std::vector<double> samples)
{
    update_values(std::move(samples));
}

void statistics_t::update_values(std::vector<double> samples)
{
    std::sort(samples.begin(), samples.end());

    std::size_t n = samples.size();    
    _total = std::accumulate(samples.begin(), samples.end(), 0.0);
    _mean = _total / n;
    _min = samples.front();
    _max = samples.back();

    auto fn = [=](double acc, double x) {
        auto diff = (x-_mean); return acc + diff * diff; };
    double sum_sq = std::accumulate(samples.begin(), samples.end(), 0.0, fn);
    
    _var = sum_sq / n;
    _stddev = std::sqrt(_var);

    if (samples.size() % 2) {
        _median = samples[std::floor(n/2)];
    } else {
        auto mid = static_cast<std::size_t>(n/2);
        _median = (samples[mid-1] + samples[mid])/2.0;
    } 
}

double statistics_t::mean() const
{
    return _mean;
}

double statistics_t::median() const
{
    return _median;
}

double statistics_t::max() const
{
    return _max;
}

double statistics_t::min() const
{
    return _min;
}

double statistics_t::var() const
{
    return _var;
}

double statistics_t::stddev() const
{
    return _stddev;
}

double statistics_t::total() const
{
    return _total;
}

statistics_t::operator std::string() const
{
    std::string out{};
    out += (std::string{"mean: "} + std::to_string(_mean));
    out += std::string{"  median: "} + std::to_string(_median);
    out += std::string{"  stddev: "} + std::to_string(_stddev);
    out += std::string{"  var: "} + std::to_string(_var);
    out += std::string{"  max: "} + std::to_string(_max);
    out += std::string{"  min: "} + std::to_string(_min);
    return out;
}

summarize_statistics_t::summarize_statistics_t(
                             const std::vector<statistics_t>& samples)
{
    update_values(samples);
}

summarize_statistics_t::summarize_statistics_t(
                            std::initializer_list<statistics_t> samples)
    :
    summarize_statistics_t{std::vector<statistics_t>{ samples.begin(),
        samples.end()}}
{
    
}

void summarize_statistics_t::update_values(
                             const std::vector<statistics_t>& samples)
{
    using namespace std;
    auto f_min = [](const statistics_t& s) { return s._min; };
    auto f_max = [](const statistics_t& s) { return s._max; };
    auto f_mean = [](const statistics_t& s) { return s._mean; };
    auto f_median = [](const statistics_t& s) { return s._median; };
    auto f_stddev = [](const statistics_t& s) { return s._stddev; };
    auto f_var = [](const statistics_t& s) { return s._var; };

    _min.first = mean(samples, f_min);
    _min.second = stddev(samples, _min.first, f_min);

    _max.first = mean(samples, f_max);
    _max.second = stddev(samples, _max.first, f_max);
    
    _mean.first = mean(samples, f_mean);
    _mean.second = stddev(samples, _mean.first, f_mean);
    
    _median.first = mean(samples, f_median);
    _median.second = stddev(samples, _median.first, f_median);
    
    _stddev.first = mean(samples, f_stddev);
    _stddev.second = stddev(samples, _stddev.first, f_stddev);
    
    _var.first = mean(samples, f_var);
    _var.second = stddev(samples, _var.first, f_var);
}

std::pair<double, double> summarize_statistics_t::min() const
{
    return _min;
}

std::pair<double, double> summarize_statistics_t::max() const
{
    return _max;
}

std::pair<double, double> summarize_statistics_t::mean() const
{
    return _mean;
}

std::pair<double, double> summarize_statistics_t::median() const
{
    return _median;
}

std::pair<double, double> summarize_statistics_t::var() const
{
    return _var;
}

std::pair<double, double> summarize_statistics_t::stddev() const
{
    return _stddev;
}

double summarize_statistics_t::mean(const std::vector<statistics_t>& samples,
                                    double(*fn)(const statistics_t&)) const
{
    std::size_t n = samples.size();
    double total = 0.0;
    auto sum = [=](double x, const statistics_t& s) { return x + fn(s); };
    total += std::accumulate(samples.begin(), samples.end(), 0.0, sum);
    return total / n;
}

double summarize_statistics_t::stddev(const std::vector<statistics_t>& samples,
                                      double mean,
                                      double(*fn)(const statistics_t&)) const
{
    auto square_diff = [=](double acc, const statistics_t& s) -> double {
        double diff = fn(s)-mean; return acc + diff * diff;
    };

    auto b = samples.begin(), e = samples.end();
    double squared_diff = std::accumulate(b, e, 0.0, square_diff);
    return std::sqrt(squared_diff / samples.size());
}
                                                                       

// #include <iostream>

// using namespace std;

// int main()
// {
//     std::vector<double> dt1{1, 2.3, .3, -2, 5, 9};
//     std::vector<double> dt2{2, 2, -2, 3.4, -3, 4};
//     std::vector<double> dt3{2, 4, 2, 0, -3, 4};
//     std::vector<double> dt4{2, 3, -.2, -5, 3, .4};
//     std::vector<double> dt5{1, 0, -.2, 10, -13, 4};
//     std::vector<double> dt6{2, 3, -2, 3.4, -3, 4};
    
//     statistics_t stats1{dt1};
//     statistics_t stats2{dt2};
//     statistics_t stats3{dt3};
//     statistics_t stats4{dt4};
//     statistics_t stats5{dt5};
//     statistics_t stats6{dt6};
//     std::vector<statistics_t> vs{stats1, stats2, stats3, stats4, stats5, stats6};
//     summarize_statistics_t sst{vs};
//     cout << sst.min().first << " " << sst.min().second << endl;
//     return 0;
// }
