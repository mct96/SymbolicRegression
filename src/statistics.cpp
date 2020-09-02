#include "statistics.hpp"

#include <cmath>
#include <numeric>
#include <algorithm>

statistics_t::statistics_t(std::vector<double> values)
    :
    _values{values}
{
    update_values();
}

void statistics_t::update_values()
{
    std::size_t N = _values.size();
    std::sort(_values.begin(), _values.end());
    
    auto b = _values.begin(), e = _values.end();

    _total = std::accumulate(b, e, 0.0);
    _mean = _total / N;

    auto minmax = std::minmax_element(b, e);
    _min = *minmax.first;
    _max = *minmax.second;

    double sum_sq = 0.0;
    std::for_each(b, e, [&](double v) mutable {
        sum_sq += std::pow((v-_mean), 2); });
    
    _var = sum_sq / N;
    _stddev = std::sqrt(_var);

    if (_values.size() % 2) {
        _median = _values[std::floor(N/2)];
    } else {
        auto mid = static_cast<std::size_t>(N/2);
        _median = (_values[mid-1] + _values[mid])/2.0;
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
