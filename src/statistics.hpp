#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iterator>

class statistics_t
{
public:
    statistics_t() = default;
    statistics_t(const statistics_t&) = default;
    statistics_t& operator=(const statistics_t&) = default;
    
    template <typename Iter>
    statistics_t(Iter begin, Iter end);

    template <typename Iter>
    void update_values(Iter begin, Iter end);

    statistics_t merge(statistics_t stats, std::size_t n) const;
    
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
          

template <typename Iter>
statistics_t::statistics_t(Iter begin, Iter end)
    :
    statistics_t{}
{
    update_values(begin, end);
}

template <typename Iter>
void statistics_t::update_values(Iter begin, Iter end)
{
    using Tp = typename std::iterator_traits<Iter>::value_type;
    std::vector<Tp> buffer{};
    std::copy(begin, end, std::back_inserter(buffer));
    std::size_t n = buffer.size();
    
    std::sort(buffer.begin(), buffer.end());

    _total = std::accumulate(buffer.begin(), buffer.end(), 0.0);
    _mean = _total / n;

    auto minmax = std::minmax_element(buffer.begin(), buffer.end());
    _min = *minmax.first;
    _max = *minmax.second;

    double sum_sq = 0.0;
    std::for_each(buffer.begin(), buffer.end(), [&](double v) mutable {
        auto x = v - _mean; sum_sq += x * x; });
    
    _var = sum_sq / n;
    _stddev = std::sqrt(_var);

    if (buffer.size() % 2) {
        _median = buffer[std::floor(n/2)];
    } else {
        auto mid = static_cast<std::size_t>(n/2);
        _median = (buffer[mid-1] + buffer[mid])/2.0;
    } 
}

