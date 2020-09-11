#include <numeric>
#include <algorithm>
#include <cmath>

constexpr double lo = -1000.0;
constexpr double hi = 1000.0;
constexpr double near_zero = 0.001;

double stan(double x)
{
    double c = std::cos(x);
    double result =  c != 0 ? std::tan(x) :
        std::sin(x)/std::copysign(near_zero, c);

    return std::clamp(result, lo, hi);
}

double ssqrt(double x)
{
    return std::clamp(std::sqrt(std::max(x, near_zero)), lo, hi);
}

double sinv(double x)
{
    double result = x != 0 ? 1/x : 1/std::copysign(near_zero, x);
    return std::clamp(result, lo, hi);
}

double slog(double x)
{
    double result = x > 0 ? std::log(x) : std::log(near_zero);
    return std::clamp(result, lo, hi);
}

double slog10(double x)
{
    double result = std::log10(std::max(x, near_zero));
    return std::clamp(result, lo, hi);
}

double scosh(double x)
{
    return std::clamp(std::cosh(x), lo, hi);
}

double ssinh(double x)
{
    return std::clamp(std::sinh(x), lo, hi);
}

double sdiv(double l, double r)
{
    double result = r != 0 ? l/r : l/std::copysign(near_zero, r);
    return std::clamp(result, lo, hi);
}

double spow(double l, double r)
{
    double base = std::max(l, 0.0);
    double power = std::clamp(r, -5.0, 5.0);
    return std::clamp(std::pow(base, power), lo, hi);
}
