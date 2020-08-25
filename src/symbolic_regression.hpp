#pragma once

#include <vector>

namespace sr
{
    
class state_t
{
public:
    std::size_t generation() const;
    std::size_t population_sz() const;
    double max_fitness() const;
    double min_fitness() const;
    double avg_fitness() const;

    void reset();
private:
    std::size_t _generation;
    std::size_t _population_sz;
    double _max_fitness;
    double _min_fitness;
    double _avg_fitness;
};

enum class selection_method_t{ roulette_wheel, tournament };

class parameters_t
{
public:
    parameters_t();
    ~parameters_t();

    void population_sz(std::size_t value);
    std::size_t population_sz() const;
    
    void max_generation(std::size_t value);
    std::size_t max_generation() const;

    void max_depth(std::size_t value);
    std::size_t max_depth() const;

    void threshold(double value);
    double threshold() const;

    void prob_matation_one_point(double prob);
    double prob_mutation_one_point() const;

    void prob_matation_expansion(double prob);
    double prob_mutation_expasion() const;

    void prob_matation_reduction(double prob);
    double prob_mutation_reduction() const;

    void prob_crossover(double prob);
    double prob_crossover() const;

    void selection_method(selection_method_t method);
    selection_method_t selection_method() const;

    void eletism(bool enable);
    bool eletism() const;
    
    void reset();
    
private:
    std::size_t _population_sz;
    std::size_t _max_depth;
    std::size_t _max_generation;

    double _threshold;
    
    double _prob_m_one_point;
    double _prob_m_expansion;
    double _prob_m_reduction;

    double _prob_crossover;

    selection_method_t _selection_method;

    bool _eletism;
};

class symbolic_regression_t
{
public:
    symbolic_regression_t();
    ~symbolic_regression_t();

private:

};

}
