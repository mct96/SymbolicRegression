#pragma once

#include <vector>
#include <iostream>

//namespace sr
//{
struct gene_t;

enum class selection_method_t{ roulette_wheel, tournament };
enum class generation_method_t{ full, growth };
enum class error_metric_t { mse, rmse, mae };

class state_t;
class parameters_t;
class gp_operators_t;
class symbolic_regression_t;

inline std::size_t lchild(std::size_t pos) { return 2 * pos + 1; };
inline std::size_t rchild(std::size_t pos) { return 2 * pos + 2; };
inline std::size_t parent(std::size_t pos) {
    return static_cast<int>((pos + 1)/2 - 1); };

// individual_t represents an individual: each gene can be a function,
// an operator or a variable.
using individual_t = std::vector<gene_t>;
using individuals_t = std::vector<individual_t>;

// variable_t are terminals. Its size should be defined when the dataset
// is loaded. For example, if dataset has N column + y output, the 
// variables_t's length should be N.
using data_t = std::vector<double>;

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

// state_t class is used to show the state of convergence of the algorithm.
// its purpouse is to plot the progress in realtime.
class state_t
{
    friend class symbolic_regression;
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

// parameters_t are the parameters of the algorithm.
class parameters_t
{
public:
    parameters_t();
    ~parameters_t();

    // size of population.
    void population_sz(std::size_t value);
    std::size_t population_sz() const;

    // generation limit.
    void max_generation(std::size_t value);
    std::size_t max_generation() const;

    // max depth of an individual. (consider to use 2^n)
    void max_depth(std::size_t value);
    std::size_t max_depth() const;

    // tolerable error.
    void threshold(double value);
    double threshold() const;

    // probability of one point mutation.
    void prob_matation_one_point(double prob);
    double prob_mutation_one_point() const;

    // probability of expansion mutation.
    void prob_matation_expansion(double prob);
    double prob_mutation_expasion() const;

    // probability of reduction mutation.
    void prob_matation_reduction(double prob);
    double prob_mutation_reduction() const;

    // probability of crossover.
    void prob_crossover(double prob);
    double prob_crossover() const;

    // selection method.
    void selection_method(selection_method_t method);
    selection_method_t selection_method() const;

    // enable eletism?
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



// gp_operators_t implementation of operators for gp.
class gp_operators_t
{
public:
    // eval fitness. (used to sort).
    double fitness(const individual_t& individual,
                   const data_t& variables,
                   double target);

    // MSE of population.
    double population_error(const individuals_t& population,
                            const data_t& variables,
                            double target,
                            error_metric_t error_metric = error_metric_t::mae);

    // 2 types of individual's generation (growth, full).
    individual_t gen_individual(generation_method_t generation_method);

    // 2 types of selection (roullet wheel, tournement).
    individuals_t selection_rw(const individuals_t& population);
    individuals_t selection_t(const individuals_t& population, std::size_t sz);

    // crossover: copy subtree.
    individuals_t crossover(const individual_t& p1, const individual_t& p2);

    // 3 types of mutation (one point, expansion, reduction). Note: one point
    // mutation change the code not the class, except for variable that can
    // exchange with constant.
    void mutation_op(individual_t& individual, double p);
    void mutation_ex(individual_t& individual, double p);
    void mutation_rd(individual_t& individual, double p);

private:
    void copy_subtree(const individual_t& src, individual_t& dst);
};
    
class symbolic_regression_t
{
public:
    symbolic_regression_t(parameters_t params);
    ~symbolic_regression_t();


private:

};

//}
