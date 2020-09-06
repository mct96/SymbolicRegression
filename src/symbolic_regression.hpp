#pragma once

#include <vector>
#include <cmath>
#include <iostream>

#include "random.hpp"

//namespace sr
//{

static std::size_t n_functions = 5, n_operators = 4; 

struct gene_t;

// individual_t represents an individual: each gene can be a function,
// an operator or a variable.
using individual_t = std::vector<gene_t>;

// individual and score (fitness).
using individuals_t = std::vector<std::pair<individual_t, double>>;

// variable_t are terminals. Its size should be defined when the dataset
// is loaded. For example, if dataset has N column + y output, the 
// variables_t's length should be N.
using data_t = std::vector<double>;
using training_set_t = std::vector<std::pair<std::vector<double>, double>>;

enum class selection_method_t{ roulette_wheel, tournament };
enum class generation_method_t{ full, grow, ramped_hh };
enum class error_metric_t { mse, rmse, mae };

// consteval for C++ 20.
inline std::size_t lchild(std::size_t pos) { return 2 * pos + 1; };
inline std::size_t rchild(std::size_t pos) { return 2 * pos + 2; };
inline std::size_t parent(std::size_t pos) {
    return static_cast<int>((pos + 1)/2 - 1); };



// the eval_* functions are responsable for evaluate a tree.
double eval_function(gene_t func, double x);
double eval_operator(gene_t oper, double lop, double rop);
double eval(const individual_t& individual, std::size_t pos,
            const data_t& variables);

void print     (std::ostream& out, const individual_t& u, std::size_t pos); 
void print_func(std::ostream& out, const individual_t& u, std::size_t pos);
void print_oper(std::ostream& out, const individual_t& u, std::size_t pos);
void print_var (std::ostream& out, const individual_t& u, std::size_t pos);
void print_cons(std::ostream& out, const individual_t& u, std::size_t pos);
std::ostream& operator<<(std::ostream& out, const individual_t& u); 

class state_t;
class parameters_t;
class gp_operators_t;
class symbolic_regression_t;

// a gene is represented by two portions: class and code. Class can be an op,
// func, var, nil. "op" represents an binary operator (+, -, /, *, ^); "func"
// represents a function (trigonometric functions, hiberbolic functions,
// logarithmics functions); "var" represents a variable and the domain is defi-
// ned by the dataset; "cons" represents numerics constants (their value is
// random generated with a distribution with two peaks in -1 and 1).
enum class class_t: unsigned char { operator_t, function_t, variable_t,
                                    constant_t, nill_t };
// a gene can have from 0-7 class and 0-32 different codes.
struct gene_t {
    class_t _class = class_t::nill_t;
    union { unsigned short _code = 0; float _value; };
};

std::size_t get_depth(std::size_t pos);
std::size_t get_max_depth(const individual_t& pos);

// state_t class is used to show the state of convergence of the algorithm.
// its purpouse is to plot the progress in realtime.
class state_t
{
    friend class symbolic_regression_t;
public:
    state_t();
    
    std::size_t generation() const;
    std::size_t population_sz() const;
    std::size_t unique_individuals() const;
    std::size_t better_than_med() const;
        
    double max_fitness() const;
    double min_fitness() const;
    double avg_fitness() const;
    double med_fitness() const;
    double std_fitness() const;

    operator std::string() const;
    
    void reset();
private:
    std::size_t _generation;
    std::size_t _population_sz;
    std::size_t _unique_individuals;
    std::size_t _better_than_med;
    
    double _max_fitness;
    double _min_fitness;
    double _avg_fitness;
    double _med_fitness;
    double _std_fitness;
};

// parameters_t are the parameters of the algorithm.
class parameters_t
{
    friend class symbolic_regression_t;
public:
    parameters_t(std::size_t n_vars,
                 std::vector<int> seeds,
                 std::size_t population_sz = 1000,
                 std::size_t max_depth = 5,
                 std::size_t max_generation = 50,
                 double threshold = 1.0,
                 double prob_crossover = 0.8,
                 double prob_mutation = 0.15,
                 double prob_op_mutation = 0.30,
                 double prob_reproduction = 0.05,
                 selection_method_t selection_method = selection_method_t::tournament,
                 std::size_t k = 5,
                 generation_method_t generation_method = generation_method_t::ramped_hh,
                 error_metric_t error_metric = error_metric_t::mae,
                 bool eletism = false);
    
    ~parameters_t();

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

    void n_vars(std::size_t n);
    std::size_t n_vars() const;

    void seeds(std::vector<int> sds);
    std::vector<int> seeds() const;
    
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
    void selection_method(selection_method_t method);
    selection_method_t selection_method() const;
    void tournament(std::size_t k);
    std::size_t tournament() const;
    
    // generation method.
    void generation_method(generation_method_t method);
    generation_method_t generation_method() const;

    // error metric.
    void error_metric(error_metric_t metric);
    error_metric_t error_metric() const;
    
    // enable eletism?
    void eletism(bool enable);
    bool eletism() const;
    
    void reset();
    
private:    
    std::vector<int> _seeds;
    
    std::size_t _population_sz;
    std::size_t _max_depth;
    std::size_t _max_generation;
    std::size_t _n_vars;

    double _threshold;

    double _prob_crossover;
    double _prob_mutation;
    double _prob_op_mutation;
    double _prob_reproduction;

    selection_method_t _selection_method;
    std::size_t _k;
    
    generation_method_t _generation_method;
    error_metric_t _error_metric;

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
    gp_operators_t(std::vector<int> seeds);
    
    // eval fitness. (used to sort).
    double fitness(const individual_t& individual,
                   const training_set_t& input_data,
                   error_metric_t metric);

    void population_fitness(individuals_t& population,
                            const training_set_t& input_data,
                            error_metric_t metric);
    
    // 2 types of individual's generation (growth, full).
    individual_t full_gen(std::size_t max_depth, std::size_t n_vars);
    individual_t grow_gen(std::size_t max_depth, std::size_t n_vars);

    individual_t gen_individual(generation_method_t generation_method,
                                std::size_t max_depth,
                                std::size_t n_vars);

    // 2 types of selection (roullet wheel, tournement).
    individual_t selection(const individuals_t& population,
                           selection_method_t sm,
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
                              selection_method_t sm);

    std::size_t count_unique(const individuals_t& population);
private:

    
    void grow_gen_recursive(individual_t& individual,
                            std::size_t n_vars,
                            std::size_t max_depth,
                            std::size_t cur_point);
    
    std::size_t rd_gene(const individual_t& individual);
    
    gene_t rd_terminal(std::size_t n_vars);

    void copy_subtree(const individual_t& src,
                      std::size_t src_point,
                      individual_t& dst,
                      std::size_t dst_point);

    void clear_subtree(individual_t& individual, std::size_t point);

    random_t _rd;
};
    
class symbolic_regression_t
{
public:
    symbolic_regression_t(parameters_t params,
                          const training_set_t& data);
    
    ~symbolic_regression_t();

    const state_t& state();

    void parameters(parameters_t params);

    parameters_t parameters() const;
            
    std::string report() const;

    void train(bool verbose = true);
    
private:
    void clear_line() const;
    
    void print_progress() const;
    
    void update_state();
    
    void initialize_population();

    void next_generation();

    void do_crossover(individuals_t& individuals, std::size_t n);

    void do_mutation(individuals_t& individuals, std::size_t n);

    void do_reproduction(individuals_t& n_individuals, std::size_t n);

    bool check_stop_condition() const;
    
    std::vector<state_t> _states;
    state_t _cur_state;
    individuals_t _population;
    parameters_t  _parameters;
    gp_operators_t _gpo;
    const training_set_t& _input_data;
};

//}
