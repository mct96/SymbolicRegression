#pragma once

#include <random>

constexpr int seeds_v[] = {7, 77, 777, 7777, 77777};
std::seed_seq seeds(seeds_v, seeds_v + 1);
std::mt19937 gen{seeds};

int rd_integer(int N) {
    std::uniform_int_distribution<> dist{0, N-1};
    return dist(gen);
}

// multimodal uniform distribution. Two peaks: -1 and 1.
double rd_value() {
    std::uniform_int_distribution<> _switch{-1, 1};
    double c = _switch(gen);
    std::normal_distribution<> dist{c, 1};
    return dist(gen);
}

inline int rd_function(int N) {
    return rd_integer(N);
}

inline int rd_operator(int N) {
    return rd_integer(N);
}

inline int rd_variable(int N) {
    return rd_integer(N);
}

inline int rd_binary() {
    return rd_integer(2);
}
