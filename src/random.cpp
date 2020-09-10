#include "random.hpp"

#include <utility>

 // <-- It's a singleton $%$#$%@%!!!
std::mt19937 random_t::_gen{};

random_t::random_t(std::vector<int> seed)
    :
    _seed_seq{new std::seed_seq(seed.begin(), seed.end())}
{
    _gen.seed(*_seed_seq);
}


void random_t::use_seed(std::vector<int> seed)
{
    _seed_seq.reset(new std::seed_seq(seed.begin(), seed.end()));
    _gen.seed(*_seed_seq);
}

int random_t::var(int from, int to) const
{
    std::uniform_int_distribution<> dist{from, to-1};
    return dist(_gen);
}

int random_t::binary() const
{
    std::uniform_int_distribution<> dist{0, 1};
    return dist(_gen);
}

int random_t::func(int from, int to) const
{
    std::uniform_int_distribution<> dist{from, to - 1};
    return dist(_gen);
}

int random_t::oper(int from, int to) const
{
    std::uniform_int_distribution<> dist{from, to - 1};
    return dist(_gen);
}

bool random_t::choose() const
{
    return binary() == 1;
}

double random_t::prob() const
{
    std::uniform_real_distribution<> dist{0.0, 1.0};
    return dist(_gen);
}

double random_t::value(double m1, double m2, double stddev) const
{
    std::normal_distribution<> dist{choose() ? m1 : m2, stddev};
    return dist(_gen);
}

int random_t::integer(int from, int to) const
{
    std::uniform_int_distribution<> dist{from, to};
    return dist(_gen);
}

// random_t::operator std::string() const
// {
//     if (_seeds_v.size() == 0) return std::string{"seeds: Є."};
    
//     std::string out{"seeds: "};
//     for (int i = 0; i < _seeds_v.size()-1; ++i)
//         out += std::to_string(_seeds_v[i]) + ", ";
//     out += std::to_string(_seeds_v.back()) + ".";

//     return out;
// }


// #include <iostream>
// using namespace std;

// int main()
// {
//     random_t rd{{0, 1, 2}};
//     random_t rd2 = rd;
//     cout << rd2.integer(0, 100) << endl;
//     return 0;
// }
