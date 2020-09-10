#include "statistics.hpp"

// template <typename Iter>
// statistics_t::statistics_t(Iter begin, Iter end)
// {
//     update_values(begin, end);
// }


statistics_t statistics_t::merge(statistics_t prev, std::size_t n) const
{
    if (n < 1) throw std::invalid_argument{"statistics_t::merge"};
    
    auto fn = [=](double mean, double a) { return (mean * n + a)/(n + 1); };

    statistics_t merged{};
    merged._mean = fn(prev._mean, _mean);
    merged._median = fn(prev._median, _median);
    merged._max = fn(prev._max, _max);
    merged._min = fn(prev._min, _min);
    merged._var = fn(prev._var, _var);
    merged._stddev = fn(prev._stddev, _stddev);
    merged._total = fn(prev._total, _total);

    return merged;
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

// #include <iostream>
// #include <forward_list>
// using namespace std;

// int main()
// {
//     std::forward_list<double> dt1{{1, 2.3, .3, -2, 5, 9}};
//     std::forward_list<double> dt2{{2, 3, -2, 3.4, -3, 4}};
//     statistics_t stats1{dt1.begin(), dt1.end()};
//     statistics_t stats2{dt2.begin(), dt2.end()};
//     cout << static_cast<std::string>(stats1) << endl;
//     cout << static_cast<std::string>(stats2) << endl;
//     cout << static_cast<std::string>(stats2.merge(stats1, 1)) << endl;
    
//     return 0;
// }
