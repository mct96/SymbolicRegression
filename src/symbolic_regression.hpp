#pragma once

#include <vector>
#include <cmath>
#include <iostream>

#include "random.hpp"
#include "node.hpp"
#include "common_types.hpp"
//namespace sr
//{

enum class gen_t {full_t, grow_t, ramped_t};
enum class sel_t {roul_t, tour_t};
enum class err_t {mae_t, mse_t, rmse_t};

class state_t;
class params_t;
class gp_operators_t;
class symbolic_regression_t;

// state_t class is used to show the state of convergence of the algorithm.
// its purpouse is to plot the progress in realtime.
class state_t
{
    friend class symbolic_regression_t;
public:
    state_t(std::vector<double> pop_fitness);
    
    std::size_t generation() const;
    std::size_t population_sz() const;
    std::size_t unique_individuals() const;
    std::size_t better_than_previous() const;
        
    double min() const;
    double max() const;
    double mean() const;
    double stddev() const;
    double var() const;
    double median() const;
    double total() const;

    operator std::string() const;
    
    void reset();
private:
    std::size_t _generation;
    std::size_t _population_sz;
    std::size_t _unique_individuals;
    std::size_t _better_than_previous;
    
    statistics_t _stats;
};

// params_t are the params of the algorithm.
class params_t
{
    friend class symbolic_regression_t;
public:
    params_t(std::vector<std::vector<int>> seeds,
                 std::size_t population_sz = 500,
                 std::size_t max_depth = 7,
                 std::size_t max_generation = 100,
                 double threshold = 0.001,
                 double prob_crossover = 0.8,
                 double prob_mutation = 0.15,
                 double prob_op_mutation = 0.15,
                 double prob_reproduction = 0.05,
                 sel_t selection_method = sel_t::tour_t,
                 std::size_t k = 5,
                 gen_t generation_method = gen_t::ramped_t,
                 err_t error_metric = err_t::mae_t,
                 bool eletism = false);
    
    ~params_t();

    operator std::string() const;
    
    // size of population.
    void population_sz(std::size_t value);
    std::size_t population_sz() const;

    // generation limit.
    void max_generation(std::size_t value);
    std::size_t max_generation() const;

    // max depth of an individual. (consider to use 2^n)
    void max_depth(std::size_t value);
    std::size_t max_depth() const;

    void seeds(std::vector<std::vector<int>> sds);
    std::vector<std::vector<int>> seeds() const;
    
    // tolerable error: [0, 1].
    void threshold(double value);
    double threshold() const;

    // probability of one point mutation.
    void prob_mutation(double prob);
    double prob_mutation() const;

    void prob_one_point_mutation(double prob);
    double prob_one_point_mutation() const;
    
    // probability of crossover.
    void prob_crossover(double prob);
    double prob_crossover() const;

    // prob = 1 - prob_mutation - prob_crossover.
    void prob_reproduction(double prob);
    double prob_reproduction() const;
    
    // selection method.
    void selection_method(sel_t method);
    sel_t selection_method() const;
    void tournament(std::size_t k);
    std::size_t tournament() const;
    
    // generation method.
    void generation_method(gen_t method);
    gen_t generation_method() const;

    // error metric.
    void error_metric(err_t metric);
    err_t error_metric() const;
    
    // enable eletism?
    void eletism(bool enable);
    bool eletism() const;
    
    void reset();
    
private:    
    const std::vector<const std::vector<int>> _seeds;
    
    std::size_t _population_sz;
    std::size_t _max_depth;
    std::size_t _max_generation;
    std::size_t _n_vars;

    double _threshold;

    double _prob_crossover;
    double _prob_mutation;
    double _prob_op_mutation;
    double _prob_reproduction;

    sel_t _sel;
    std::size_t _k;
    
    gen_t _gen;
    err_t _err;

    bool _eletism;

    std::string selec_met_str() const;
    std::string gen_met_str() const;
    std::string err_met_str() const;
    std::string seeds_str() const;
};

// gp_operators_t implementation of operators for gp.
class gp_operators_t
{
public:
    gp_operators_t(random_t& random,
                   individual_handler_t& handler);
    
    // eval fitness. (used to sort).
    double fitness(const individual_t& individual,
                   const std::vector<entry_t>& input_data,
                   err_t metric);

    void population_fitness(individuals_t& population,
                            const std::vector<entry_t>& input_data,
                            err_t metric);
    
    // 2 types of individual's generation (growth, full).
    individual_t full_gen(std::size_t max_depth);
    individual_t grow_gen(std::size_t max_depth);

    individual_t gen_individual(gen_t generation_method,
                                std::size_t max_depth,
                                std::size_t n_vars);

    // 2 types of selection (roullet wheel, tournement).
    individual_t selection(const individuals_t& population,
                           sel_t sm,
                           std::size_t k);

    individual_t selection_rw(const individuals_t& population);
    individual_t selection_t(const individuals_t& population, std::size_t k);

    
    individual_t crossover(const individual_t& p1,
                           const individual_t& p2,
                           std::size_t max_depth);

    
    individual_t mutation_op(const individual_t& individual,
                             std::size_t n_vars,
                             double prob);
    
    individual_t mutation_sb(const individual_t& individual,
                             std::size_t n_vars);

    individual_t reproduction(const individuals_t& population,
                              sel_t sm);

    std::size_t count_unique(const individuals_t& population);

private:   
    void grow_gen_recursive(individual_t& individual,
                            std::size_t n_vars,
                            std::size_t max_depth,
                            std::size_t cur_point);
    
    random_t& _rd;
    individual_handler_t& _hdl;
};
    
class symbolic_regression_t
{
public:
    symbolic_regression_t(params_t params,
                          random_t& random,
                          individual_handler_t& handler);
    
    ~symbolic_regression_t();

    const state_t& state();

    void params(params_t params);

    params_t params() const;
            
    std::string report() const;

    void train(std::vector<entry_t> training_set);
    
    double test(std::vector<entry_t> test_set) const;
    
private:
    void clear_line() const;
    
    void print_progress() const;
    
    void update_state();
    
    void initialize_population(const std::vector<entry_t>& t_set);

    void next_generation(const std::vector<entry_t>& t_set);

    void do_crossover(individuals_t& individuals,
                      std::size_t n,
                      const std::vector<entry_t>& t_set);

    void do_mutation(individuals_t& individuals,
                     std::size_t n,
                     const std::vector<entry_t>& t_set);

    void do_reproduction(individuals_t& n_individuals,
                         std::size_t n,
                         const std::vector<entry_t>& t_set);

    bool check_stop_condition() const;
    
    std::vector<state_t> _states;
    individuals_t _population;
    params_t  _params;
    gp_operators_t _gpo;
    individual_handler_t _hdl;
};

//}
