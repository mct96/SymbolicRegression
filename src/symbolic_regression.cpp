#include "symbolic_regression.hpp"

#include <exception>

namespace sr
{

std::size_t state_t::generation() const
{
    return _generation;
}
    
std::size_t state_t::population_sz() const
{
    return _population_sz;
}    

double state_t::max_fitness() const
{
    return _max_fitness;
}

double state_t::min_fitness() const
{
    return _min_fitness;
}

double state_t::avg_fitness() const
{
    return _avg_fitness;
}

void state_t::reset()
{
    _generation    = 0;
    _population_sz = 0;
    _max_fitness   = 0;
    _min_fitness   = 0;
    _avg_fitness   = 0;
}

parameters_t::parameters_t()
{
}

parameters_t::~parameters_t()
{
}
    
void parameters_t::population_sz(std::size_t value)
{
    _population_sz = value;
}

std::size_t parameters_t::population_sz() const
{
    return _population_sz;
}

void parameters_t::max_generation(std::size_t value)
{
    _max_generation = value;
}

std::size_t parameters_t::max_generation() const
{
    return _max_generation;
}

void parameters_t::max_depth(std::size_t value)
{
    _max_depth = value;
}

std::size_t parameters_t::max_depth() const
{
    return _max_depth;
}

void parameters_t::threshold(double value)
{
    _threshold = value;
}

double parameters_t::threshold() const
{
    return _threshold;
}

void parameters_t::prob_matation_one_point(double prob)
{
    _prob_m_one_point = prob;
}

double parameters_t::prob_mutation_one_point() const
{
    return _prob_m_one_point;
}

void parameters_t::prob_matation_expansion(double prob)
{
    
    _prob_m_expansion = prob;
}

double parameters_t::prob_mutation_expasion() const
{
    return _prob_m_expansion;
}

void parameters_t::prob_matation_reduction(double prob)
{
    _prob_m_reduction = prob;
}

double parameters_t::prob_mutation_reduction() const
{
    return _prob_m_reduction;
}

void parameters_t::prob_crossover(double prob)
{
    _prob_crossover = prob;
}

double parameters_t::prob_crossover() const
{
    return _prob_crossover;
}

void parameters_t::selection_method(selection_method_t method)
{
    _selection_method = method;
}

selection_method_t parameters_t::selection_method() const
{
    return _selection_method;
}

void parameters_t::eletism(bool enable)
{
    _eletism = enable;
}

bool parameters_t::eletism() const
{
    return _eletism;
}
    

void parameters_t::reset()
{
    _population_sz = 0;
    _max_depth = 0;
    _max_generation = 0;
    _threshold = 0.0;
    _prob_m_one_point = 0.0;
    _prob_m_expansion = 0.0;
    _prob_m_reduction = 0.0;
    _prob_crossover = 0.0;
    _selection_method = selection_method_t::roulette_wheel;
    _eletism = false;
}

symbolic_regression_t::symbolic_regression_t()
{
}
    
symbolic_regression_t::~symbolic_regression_t()
{
}
    
}
