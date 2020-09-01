#include "symbolic_regression.hpp"
#include "random.hpp"

#include <functional>
#include <cmath>
#include <type_traits>
#include <stdexcept>
#include <iterator>


//namespace sr
//{

// safe operators and functions.
double stan(double x)
{
    if (cos(x) == 0)
        x += 0.01;

    return tan(x);
}

double slog10(double x)
{
    return log10(std::max(0.00001, x));
}

double slog(double x)
{
    return log(std::max(0.00001, x));
}

double div(double a, double b)
{
    if (b == 0) return 1000 * a;
    return a / b;
}

double eval_function(gene_t func, double x)
{
    auto code = func._code;
 
    double y = 0;
    switch (code) {
    case 0: y = sin(x); break;
    case 1: y = cos(x); break;        
    case 2: y = stan(x); break;
        //    case 3: y = sinh(x); break;
        //    case 4: y = cosh(x); break;
        //    case 5: y = tanh(x); break;
    case 3: y = slog10(x); break;
    case 4: y = slog(x); break;
    default:
        throw std::invalid_argument{"eval function"};
    }
    
    return y;
}

double eval_operator(gene_t oper, double lop, double rop)
{
    auto code = oper._code;
    
    double result = 0;
    switch (code) {
    case 0: result = lop + rop; break;
    case 1: result = lop - rop; break;
    case 2: result = lop * rop; break;
    case 3: result = div(lop, rop); break;
        //case 4: result = pow(lop, rop); break;
    default:
        throw std::invalid_argument{"eval operator"};
    }

    return result;
}

double eval(const individual_t& individual,
            std::size_t pos,
            const data_t& variables)
{
    gene_t gene = individual[pos];

    if (gene._class == class_t::operator_t) {                  // eval operator.
        double lop = eval(individual, lchild(pos), variables); // extract lop.
        double rop = eval(individual, rchild(pos), variables); // extract rop.
        return eval_operator(gene, lop, rop);
        
    } else if (gene._class == class_t::function_t) {         // eval function.
        double x = eval(individual, lchild(pos), variables); // eval arguments.
        return eval_function(gene, x);
        
    } else if (gene._class == class_t::variable_t) { // eval variable.
        return variables[gene._code];
        
    } else if (gene._class == class_t::constant_t) { // eval constant.
        return gene._value;
        
    } else {
        throw std::invalid_argument{"general eval"};
    }
}

void print_func(std::ostream& out, const individual_t& u, std::size_t pos)
{
    gene_t func = u[pos];

    switch (func._code) {
    case 0: out << "sin"  ; break;        
    case 1: out << "cos"  ; break;
    case 2: out << "tan"  ; break;
        //    case 3: out << "sinh" ; break;
        //    case 4: out << "cosh" ; break;
        //    case 5: out << "tanh" ; break;
    case 3: out << "log10"; break;
    case 4: out << "log"  ; break;
    default:
        throw std::invalid_argument{"print function"};
    }

    out << "("; print(out, u, lchild(pos)); out << ")";
}

void print_oper(std::ostream& out, const individual_t& u, std::size_t pos)
{
    gene_t oper = u[pos], loper = u[lchild(pos)], roper = u[rchild(pos)];
    
    // print left operand. surround lower precedence operator with parentheses.
    if (loper._class == class_t::operator_t && loper._code <= oper._code) {
        out << "("; print(out, u, lchild(pos)); out << ")";
    } else {
        print(out, u, lchild(pos));
    }

    // print operator.
    switch (oper._code) {
    case 0: out << " + "; break;
    case 1: out << " - "; break;
    case 2: out << " * "; break;
    case 3: out << " / "; break;
        //case 4: out << " ^ "; break;
    default:
        throw std::invalid_argument{"eval operator"};
    }

    // print right operand. surround lower precedence operator with parentheses.
    if (roper._class == class_t::operator_t && roper._code <= oper._code) {
        out << "("; print(out, u, rchild(pos)); out << ")";
    } else {
        print(out, u, rchild(pos));
    }
}

void print_var(std::ostream& out, const individual_t& u, std::size_t pos)
{
    gene_t var = u[pos];
    out << "x_" << static_cast<int>(var._code);
}

void print_cons(std::ostream& out, const individual_t& u, std::size_t pos)
{
    gene_t cons = u[pos];
    out << cons._value;
}

void print(std::ostream& out, const individual_t& u, std::size_t pos)
{
    // just boot.
    gene_t head = u[pos];

    switch (head._class) {
    case class_t::operator_t: print_oper(out, u, pos); break;
    case class_t::function_t: print_func(out, u, pos); break;
    case class_t::variable_t: print_var (out, u, pos); break;
    case class_t::constant_t: print_cons(out, u, pos); break;
    default:
        throw std::invalid_argument{"general print"};
    }
}

std::ostream& operator<<(std::ostream& out, const individual_t& u)
{
    print(out, u, 0);
    out << std::flush;
    return out;
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

double state_t::med_fitness() const
{
    return _med_fitness;
}

double state_t::std_fitness() const
{
    return _std_fitness;
}

void state_t::reset()
{
    _generation    = 0;
    _population_sz = 0;
    _max_fitness   = 0;
    _min_fitness   = 0;
    _avg_fitness   = 0;
}

parameters_t::parameters_t(std::size_t n_vars)
    :
    _n_vars{n_vars}
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

void parameters_t::n_vars(std::size_t n)
{
    _n_vars = n;
}

std::size_t parameters_t::n_vars() const
{
    return _n_vars;
}

void parameters_t::threshold(double value)
{
    _threshold = value;
}

double parameters_t::threshold() const
{
    return _threshold;
}

void parameters_t::prob_mutation(double prob)
{
    _prob_mutation = prob;
}

double parameters_t::prob_mutation() const
{
    return _prob_mutation;
}

void parameters_t::prob_one_point_mutation(double prob)
{
    _prob_op_mutation = prob;
}

double parameters_t::prob_one_point_mutation() const
{
    return _prob_op_mutation;
}

void parameters_t::prob_crossover(double prob)
{
    _prob_crossover = prob;
}

double parameters_t::prob_crossover() const
{
    return _prob_crossover;
}

double parameters_t::prob_reproduction() const
{
    return _prob_reproduction;
}

void parameters_t::selection_method(selection_method_t method)
{
    _selection_method = method;
}

selection_method_t parameters_t::selection_method() const
{
    return _selection_method;
}

void parameters_t::tournament(std::size_t k)
{
    _k = k;
}

std::size_t parameters_t::tournament() const
{
    return _k;
}

void parameters_t::generation_method(generation_method_t method)
{
    _generation_method = method;
}

generation_method_t parameters_t::generation_method() const
{
    return _generation_method;
}

void parameters_t::error_metric(error_metric_t metric)
{
    _error_metric = metric;
}

error_metric_t parameters_t::error_metric() const
{
    return _error_metric;
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
    _population_sz     = 0;
    _max_depth         = 0;
    _max_generation    = 0;
    _threshold         = 0.0;
    _prob_mutation     = 0.0;
    _prob_reproduction = 0.0;
    _prob_crossover    = 0.0;
    _selection_method  = selection_method_t::roulette_wheel;
    _eletism           = false;
}

double gp_operators_t::fitness(const individual_t& individual,
                               const training_set_t& input_data,
                               error_metric_t err_m)
{
    std::vector<double> error{};

    for (auto xy: input_data) {
        std::vector<double> x = xy.first;

        double y = xy.second;
        double y_pred = eval(individual, 0, x);
        double dy = y - y_pred;
        
        error.push_back(err_m == error_metric_t::mae ? fabs(dy) : dy * dy);
    }

    auto b = error.begin(), e = error.end();
    double me = std::accumulate(b, e, 0.0)/error.size(); // mean error.

    return (err_m == error_metric_t::rmse) ? sqrt(me) : me; 
}

void gp_operators_t::population_fitness(individuals_t& population,
                                        const training_set_t& input_data,
                                        error_metric_t err_m)
{
    for (auto& individual : population) {
        double error = fitness(individual.first, input_data, err_m);
        individual.second = error;
    }
}

individual_t gp_operators_t::full_gen(std::size_t max_depth,
                                      std::size_t n_vars)
{
    if (max_depth <= 1) return individual_t{rd_terminal(n_vars)};
    
    individual_t individual((int)pow(2, max_depth) - 1);

    for (std::size_t i = 0; i < (std::size_t)pow(2, max_depth - 2) - 1; ++i) {
        individual[i]._class = class_t::operator_t;
        individual[i]._code = rd_operator(n_operators);
    }
    
    for (std::size_t i = pow(2, max_depth - 2) - 1;
                     i < pow(2, max_depth - 1) - 1; ++i) {
        auto function_or_operator = rd_binary();

        if (function_or_operator) {
            individual[i]._class = class_t::function_t;
            individual[i]._code = rd_function(n_functions);
            individual[lchild(i)] = rd_terminal(n_vars);
        } else {
            individual[i]._class = class_t::operator_t;
            individual[i]._code = rd_operator(n_operators);
            individual[lchild(i)] = rd_terminal(n_vars);
            individual[rchild(i)] = rd_terminal(n_vars);
        }
    }

    return individual;
}

individual_t gp_operators_t::grow_gen(std::size_t max_depth,
                                      std::size_t n_vars)
{
    individual_t new_individual(pow(2, max_depth) - 1);
    clear_subtree(new_individual, 0);
    grow_gen_recursive(new_individual, n_vars, max_depth, 0);
    return new_individual;
}

individual_t gp_operators_t::gen_individual(generation_method_t gm,
                                            std::size_t max_depth,
                                            std::size_t n_vars)
{
    auto gm_f = generation_method_t::full, gm_g = generation_method_t::grow;
    if (gm == gm_f)
        return full_gen(max_depth, n_vars);
    else if (gm == gm_g)
        return grow_gen(max_depth, n_vars);
    else if (gm == generation_method_t::ramped_hh)
        return gen_individual((rd_value() < 0.5) ? gm_f : gm_g,
                              max_depth,
                              n_vars);
    return {};
}

individual_t gp_operators_t::selection(const individuals_t& population,
                                       selection_method_t sm,
                                       std::size_t k)
{
    auto sm_rw = selection_method_t::roulette_wheel,
         sm_t  = selection_method_t::tournament;

    if (sm == sm_rw)
        return selection_rw(population);
    else if (sm == sm_t)
        return selection_t(population, k);
    else
        throw std::invalid_argument{"gp_operators_t::selection"};
}
           
            

individual_t gp_operators_t::selection_rw(const individuals_t& population)
{
    auto b = population.begin(), e = population.end();
    auto sum_callable = [](auto acc, auto v) { return acc + v.second; };
    double total_fitness = std::accumulate(b, e, 0.0, sum_callable);

    double prob = rd_real();
    
    // HYPOTHESIS population should be sorted here.
    for (auto individual: population) {
        auto individual_prob = individual.second / total_fitness;
        if (prob < individual_prob)
            return individual.first;
        else
            prob -= individual_prob;
    }

    throw std::invalid_argument{"gp_operators_t::selection_rw"};
}

individual_t gp_operators_t::selection_t(const individuals_t& population,
                                          std::size_t k)
{
    auto selected = sample<typename individuals_t::value_type>(population, k);
    return selected[0].first;
}

individual_t gp_operators_t::crossover(const individual_t& p1,
                                       const individual_t& p2,
                                       std::size_t max_depth)
{
    std::size_t p1_pt = 0, p2_pt = 0;
    do {
        p1_pt = rd_gene(p1); p2_pt = rd_gene(p2);
    } while ((p1_pt + (max_depth - p2_pt - 1)) >= max_depth);
    
    std::size_t p1_parent_pt = parent(p1_pt);

    if (p1_parent_pt < 0) {
        individual_t new_individual(pow(2, max_depth) - 1);
        copy_subtree(p2, p2_pt, new_individual, 0);
        return new_individual;
    } else {
        individual_t new_individual = p1;
        clear_subtree(new_individual, p1_pt);
        copy_subtree(p2, p2_pt, new_individual, p1_pt);
        return new_individual;
    }
}

individual_t gp_operators_t::mutation_op(const individual_t& individual,
                                         std::size_t n_vars,
                                         double prob)
{
    individual_t new_individual = individual;
    for (std::size_t point = 0; point < individual.size() ; ++point) {

        if (rd_real() > prob) continue;
        
        gene_t gene = individual[point];
        switch (gene._class) {
        case class_t::operator_t:
            gene._code = rd_operator(n_operators); break;
        case class_t::function_t:
            gene._code = rd_function(n_functions); break;
        case class_t::variable_t:
            gene._code = rd_variable(n_vars); break;
        case class_t::constant_t:
            gene._value = rd_value(); break;
        }

        new_individual[point] = gene;
    }
    
    return new_individual;
}

individual_t gp_operators_t::mutation_sb(const individual_t& individual,
                                        std::size_t n_vars)
{
    std::size_t point = rd_gene(individual);
    std::size_t point_depth = get_depth(point);
    std::size_t tree_depth = get_depth(individual.size()-1);

    individual_t new_individual = individual;
    individual_t subtree = full_gen(tree_depth - point_depth, n_vars);

    clear_subtree(new_individual, point);
    copy_subtree(subtree, 0, new_individual, point);

    return new_individual;
}

void gp_operators_t::grow_gen_recursive(individual_t& individual,
                                        std::size_t n_vars,
                                        std::size_t max_depth,
                                        std::size_t point)
{
    if (max_depth - 1 == get_depth(point)) {
        if (rd_binary()) {
            individual[point]._class = class_t::variable_t;
            individual[point]._code = rd_variable(n_vars);
        } else {
            individual[point]._class = class_t::constant_t;
            individual[point]._value = rd_value();
        }
    } else {
        if (rd_binary()) {
            individual[point]._class = class_t::function_t;
            individual[point]._code = rd_operator(n_functions);
            if (rd_binary())
                individual[lchild(point)] = rd_terminal(n_vars);
            else
                grow_gen_recursive(individual, n_vars, max_depth,
                                   lchild(point));
        } else {
            individual[point]._class = class_t::operator_t;
            individual[point]._code = rd_operator(n_operators);
            if (rd_binary())
                individual[lchild(point)] = rd_terminal(n_vars);
            else
                grow_gen_recursive(individual, n_vars, max_depth,
                                   lchild(point));
            if (rd_binary())
                individual[rchild(point)] = rd_terminal(n_vars);
            else
                grow_gen_recursive(individual, n_vars, max_depth,
                                   rchild(point));
        }
    }
}

std::size_t gp_operators_t::rd_gene(const individual_t& individual)
{
    std::vector<std::size_t> points{};
    for (std::size_t pos = 0; pos < individual.size(); ++pos) {
        gene_t gene = individual[pos];
        if (gene._class != class_t::nill_t)
            points.push_back(pos);
    }

    return points[rd_integer(points.size())];
}

gene_t gp_operators_t::rd_terminal(std::size_t n_vars)
{
    gene_t gene;
    auto constant_or_variable = rd_binary();
    if (constant_or_variable) {
        gene._class = class_t::variable_t;
        gene._code = rd_variable(n_vars);
    } else {
        gene._class = class_t::constant_t;
        gene._value = rd_value();
    }

    return gene;
}

void gp_operators_t::copy_subtree(const individual_t& src,
                                 std::size_t src_point,
                                 individual_t& dst,
                                 std::size_t dst_point)
{
    dst[dst_point] = src[src_point];
    
    auto src_lc = lchild(src_point), dst_lc = lchild(dst_point);
    if (src_lc < src.size() && dst_lc < dst.size())
        copy_subtree(src, src_lc, dst, dst_lc);

    auto src_rc = rchild(src_point), dst_rc = rchild(dst_point);
    if (src_rc < src.size() && dst_rc < dst.size())
        copy_subtree(src, src_rc, dst, dst_rc);
}

void gp_operators_t::clear_subtree(individual_t& individual, std::size_t point)
{
    if (point >= individual.size())
        return;
 
    individual[point]._class = class_t::nill_t;
    clear_subtree(individual, lchild(point));
    clear_subtree(individual, rchild(point));
}

symbolic_regression_t::symbolic_regression_t(parameters_t params,
                                             const training_set_t& data)
    :
    _parameters{params},
    _input_data{data}
{

}
    
symbolic_regression_t::~symbolic_regression_t()
{
    //    report();
}

const state_t& symbolic_regression_t::state()
{
    return _state.back();
}

void symbolic_regression_t::parameters(parameters_t params)
{
    _parameters = params;
}

parameters_t symbolic_regression_t::parameters() const
{
    return _parameters;
}

void symbolic_regression_t::initialize_population()
{
    auto params = _parameters;
    auto em = params._error_metric;
    auto method = params._generation_method;    
    auto max_depth = params._max_depth;
    auto n_vars = params._n_vars;
    
    for (std::size_t i = 0; i < params._population_sz; ++i) {
        auto individual = _gpo.gen_individual(method, max_depth, n_vars);
        auto fitness = _gpo.fitness(individual, _input_data, em); 
        _population.emplace_back(individual, fitness);
    }

    auto b = _population.begin(), e = _population.end();
    std::sort(b, e, [](auto a, auto b) {return a.second < b.second;});

    update_state();
}

void symbolic_regression_t::next_generation()
{
    auto params = _parameters; 
    auto em = params._error_metric;
    auto pop_sz = params._population_sz;
    
    double pc = params._prob_crossover, pm = params._prob_mutation;
    std::size_t crossover_sz    = pc * pop_sz;
    std::size_t mutation_sz     = pm * pop_sz;
    std::size_t reproduction_sz = pop_sz - crossover_sz - mutation_sz;

    individuals_t new_pop{};
    do_crossover(new_pop, crossover_sz);
    do_mutation(new_pop, mutation_sz);
    do_reproduction(new_pop, reproduction_sz);    
    _gpo.population_fitness(new_pop, _input_data, em);

    using Tp = typename individuals_t::value_type;
    auto cmp_less = [](Tp a, Tp b) { return a.second < b.second; };
    std::sort(new_pop.begin(), new_pop.end(), cmp_less);

    if (params._eletism) {
        new_pop.insert(new_pop.end(), _population.begin(), _population.end());

        auto b = new_pop.begin(), e = new_pop.end();
        std::inplace_merge(b, b + pop_sz, e, cmp_less);

        b = new_pop.begin(), e = new_pop.end();
        new_pop.erase(b + pop_sz, e);
    }

    _population = new_pop;
    update_state();
}

void symbolic_regression_t::report()
{
    for (auto individual: _population)
        std::cout << individual.first << "  " << individual.second << std::endl;
}

void symbolic_regression_t::update_state()
{
    std::vector<double> all_fitness{};
    for (auto individual: _population)
        all_fitness.push_back(individual.second);

    state_t current_state{};
    // TODO calc statistics.
    _state.push_back(current_state);
}


void symbolic_regression_t::do_crossover(individuals_t& individuals,
                                         std::size_t n)
{
    auto params = _parameters;
    std::size_t max_depth = params._max_depth;    
    selection_method_t sm = params._selection_method;
    std::size_t k = params._k;
    error_metric_t em = params._error_metric;

    for (std::size_t i = 0; i < n; ++i) {
        auto p1 = _gpo.selection(_population, sm, k);
        auto p2 = _gpo.selection(_population, sm, k);
        individual_t offspring = _gpo.crossover(p1, p2, max_depth);
        individuals.emplace_back(offspring, 0.0);
    }
}

void symbolic_regression_t::do_mutation(individuals_t& individuals,
                                        std::size_t n)
{
    auto params = _parameters;
    auto max_depth = params._max_depth;
    auto n_vars = params._n_vars;
    auto sm = params._selection_method;
    auto k = params._k;
    auto em = params._error_metric;
    auto prob_op_mut = params._prob_op_mutation;
    std::size_t op_sz = (std::size_t)0.5 * n, sb_sz = n - op_sz;
    
    for (std::size_t i = 0; i < op_sz; ++i) {
        individual_t p = _gpo.selection(_population, sm, k);
        individual_t offspring = _gpo.mutation_op(p, n_vars, prob_op_mut);
        double fitness = _gpo.fitness(offspring, _input_data, em);
        individuals.emplace_back(offspring, fitness);
    }

    for (std::size_t i = 0; i < sb_sz; ++i) {
        individual_t p = _gpo.selection(_population, sm, k);
        individual_t offspring = _gpo.mutation_sb(p, n_vars);
        double fitness = _gpo.fitness(offspring, _input_data, em);
        individuals.emplace_back(offspring, fitness);
    }
}

void symbolic_regression_t::do_reproduction(individuals_t& individuals,
                                            std::size_t n)
{
    auto params = _parameters;
    auto sm = params._selection_method;
    auto k = params._k;
    auto em = params._error_metric;
    
    for (std::size_t i = 0; i < n; ++i) {
        individual_t offspring = _gpo.selection(_population, sm, k);
        double fitness = _gpo.fitness(offspring, _input_data, em);
        individuals.emplace_back(offspring, fitness);
    }
}

// #include <iostream>
// #include <iomanip>
// using namespace std;
// int main() {
    
//     training_set_t training{};
//     for (double x = -1; x < 1; x += 0.15)  {
//         for (double y = -1; y < 1; y += 0.15) {
//             training.emplace_back(std::vector<double>{x, y},
//                                   cos(x)*x + sin(y)/2 + x * y);
//         }
//     }
       
//     parameters_t params{2};
//     params.selection_method(selection_method_t::tournament);
//     params.error_metric(error_metric_t::rmse);
//     params.tournament(5);
//     params.max_depth(5);
//     params.eletism(true);
    
//     symbolic_regression_t sm{params, training};
//     sm.initialize_population();
//     for (auto g = 0; g < 30; ++g)
//         sm.next_generation();
//     sm.report();
//     return 0;
// }
