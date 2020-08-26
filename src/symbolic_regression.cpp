#include "symbolic_regression.hpp"

#include <cmath>
#include <exception>

//namespace sr
//{

double eval_function(char codfunc, double x)
{
    double y = 0;

    switch (codfunc) {
    case 0: y = sin(x); break;
    case 1: y = cos(x); break;        
    case 2: y = tan(x); break;
    case 3: y = sinh(x); break;
    case 4: y = cosh(x); break;
    case 5: y = tanh(x); break;
    case 6: y = log10(x); break;
    case 7: y = log(x); break;
    default:
        throw std::exception{};
    }

    return y;
}

double eval_operator(char codoper, double lop, double rop)
{
    double result = 0;

    switch (codoper) {
    case 0: result = lop + rop; break;
    case 1: result = lop - rop; break;
    case 2: result = lop * rop; break;
    case 3: result = lop / rop; break;
    case 4: result = lop ^ rop; break;
    default:
        throw std::exception{};
    }

    return result;
}

double eval(individual_t& individual,
            std::size_t pos,
            data_t& variables,
            data_t& constants)
{
    gene_t gene = individual[pos];
    std::size_t child_offset = pow(2, pos); // in a heap tree, the children are
                                // 2^pos distant from parents. 

    if (gene._class == class_t::op) {                          // eval operator.
        double lop = eval(individual, child_offset, vars);     // extract lop.
        double rop = eval(individual, child_offset + 1, vars); // extract rop.
        return eval_operator(gene._cod, lop, rop);
    } else if (gene._class == class_t::func) {           // eval function.
        double x = eval(individual, child_offset, vars); // eval arguments.
        return eval_function(gene._cod, x);
    } else if (gene._class == class_t::var) { // eval variable.
        return variables[gene._cod];
    } else if (gene._class == class_t::cons) { // eval constant.
        return constants[gene._cod];
    } else {
        throw std::exception{};
    }
}

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
    _population_sz    = 0;
    _max_depth        = 0;
    _max_generation   = 0;
    _threshold        = 0.0;
    _prob_m_one_point = 0.0;
    _prob_m_expansion = 0.0;
    _prob_m_reduction = 0.0;
    _prob_crossover   = 0.0;
    _selection_method = selection_method_t::roulette_wheel;
    _eletism          = false;
}

symbolic_regression_t::symbolic_regression_t(parameters_t params)
{
}
    
symbolic_regression_t::~symbolic_regression_t()
{
}

    
    
//}
