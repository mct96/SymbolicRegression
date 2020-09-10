#include <functional>
#include <cmath>
#include <type_traits>
#include <stdexcept>
#include <iterator>
#include <set>
#include <iomanip>

#include "symbolic_regression.hpp"

//namespace sr
//{

state_t::state_t(std::vector<double> pop_fitness)
    :
    _stats{pop_fitness.begin(), pop_fitness.end()}
{
    reset();
}

std::size_t state_t::generation() const
{
    return _generation;
}
    
std::size_t state_t::population_sz() const
{
    return _population_sz;
}    

std::size_t state_t::unique_individuals() const
{
    return _unique_individuals;
}

std::size_t state_t::better_than_previous() const
{
    return _better_than_previous;
}

state_t state_t::merge(state_t prev, std::size_t n) const
{
    if (n < 1) throw std::invalid_argument{"state_t::merge"};
    
    auto fn = [=](double mean, double a) { return (mean * n + a)/(n + 1); };
    state_t merged{};
    merged._generation = fn(prev._generation, _generation);

    merged._population_sz = fn(prev._population_sz, _population_sz);

    merged._unique_individuals = fn(prev._unique_individuals,
                                    _unique_individuals);

    merged._better_than_previous = fn(prev._better_than_previous,
                                      _better_than_previous);
    
    merged._stats = _stats.merge(prev._stats, n);

    return merged;    
}

double state_t::min() const
{
    return _stats.min();
}

double state_t::max() const
{
    return _stats.max();
}

double state_t::mean() const
{
    return _stats.mean();
}

double state_t::stddev() const
{
    return _stats.stddev();
}

double state_t::var() const
{
    return _stats.var();
}

double state_t::median() const
{
    return _stats.median();
}

double state_t::total() const
{
    return _stats.total();
}

state_t::operator std::string() const
{
    using namespace std;
    using namespace std::string_literals;
    
    string out{};
    out += to_string(_generation);
    out += ", "s + to_string(_population_sz);
    out += ", "s + to_string(_unique_individuals);
    out += ", "s + to_string(_better_than_previous);
    out += ", "s + to_string(_stats.min());
    out += ", "s + to_string(_stats.max());
    out += ", "s + to_string(_stats.mean());
    out += ", "s + to_string(_stats.stddev());
    out += ", "s + to_string(_stats.var());
    out += ", "s + to_string(_stats.median());
    return out;
}

void state_t::reset()
{
    _generation = 0;
    _population_sz = 0;
    _unique_individuals = 0;
    _better_than_previous = 0;
}

params_t::params_t(std::vector<std::vector<int>> seeds,
                   std::size_t population_sz,
                   std::size_t max_depth,
                   std::size_t max_generation,
                   double threshold,
                   double prob_crossover,
                   double prob_mutation,
                   double prob_op_mutation,
                   double prob_reproduction,
                   sel_t selection_method,
                   std::size_t k,
                   gen_t generation_method,
                   err_t error_metric,
                   bool eletism)
    :
    _seeds{seeds},
    _population_sz{population_sz},
    _max_depth{max_depth},
    _max_generation{max_generation},
    _threshold{threshold},
    _prob_crossover{prob_crossover},
    _prob_mutation{prob_mutation},
    _prob_op_mutation{prob_op_mutation},
    _prob_reproduction{prob_reproduction},
    _sel{selection_method},
    _k{k},
    _gen{generation_method},
    _err{error_metric},
    _eletism{eletism}
{
}

params_t::~params_t()
{
}

params_t::operator std::string() const
{
    using namespace std;
    using namespace std::string_literals;
    
    string out{};
    //    out += "; seeds: "s + seeds_str() + "\n";
    out += "; population size: "s + to_string(_population_sz) + "\n";
    out += "; max. depth: "s + to_string(_max_depth) + "\n";
    out += "; max. generation: "s + to_string(_max_generation) + "\n";
    out += "; # variables : "s + to_string(_n_vars) + "\n";
    out += "; threshold: "s + to_string(_threshold) + "\n";
    out += "; prob. mutation: "s + to_string(_prob_mutation) + "\n";
    out += "; prob. one point mut.: "s + to_string(_prob_op_mutation) + "\n";
    out += "; prob. crossover: "s + to_string(_prob_crossover) + "\n";
    out += "; selection method: "s + selec_met_str() + "\n";
    out += "; generation method: "s + gen_met_str() + "\n";
    out += "; error metric: "s + err_met_str() + "\n";
    out += "; eletism: "s + (_eletism ? "true" : "false") + "\n";
    return out;
}

void params_t::population_sz(std::size_t value)
{
    _population_sz = value;
}

std::size_t params_t::population_sz() const
{
    return _population_sz;
}

void params_t::max_generation(std::size_t value)
{
    _max_generation = value;
}

std::size_t params_t::max_generation() const
{
    return _max_generation;
}

void params_t::max_depth(std::size_t value)
{
    _max_depth = value;
}

std::size_t params_t::max_depth() const
{
    return _max_depth;
}

void params_t::seeds(std::vector<std::vector<int>> seeds)
{
    _seeds = seeds;
}

std::vector<std::vector<int>> params_t::seeds() const
{
    return _seeds;
}

void params_t::threshold(double value)
{
    _threshold = value;
}

double params_t::threshold() const
{
    return _threshold;
}

void params_t::prob_mutation(double prob)
{
    if (prob > 1)
        throw std::invalid_argument{"mutation probability > 1 (" +
                std::to_string(prob) + ")."};

    if (prob + _prob_crossover > 1) {
        _prob_mutation = prob;
        _prob_crossover = 1 - prob;
        _prob_reproduction = 0;
    } else {        
        _prob_mutation = prob;
        _prob_reproduction = 1 - prob - _prob_crossover;
    }
}

double params_t::prob_mutation() const
{
    return _prob_mutation;
}

void params_t::prob_one_point_mutation(double prob)
{
    _prob_op_mutation = prob;
}

double params_t::prob_one_point_mutation() const
{
    return _prob_op_mutation;
}

void params_t::prob_crossover(double prob)
{
    if (prob > 1) throw std::invalid_argument{"crossover probability > 1 (" +
                          std::to_string(prob) + ")."};

    if (prob + _prob_mutation > 1) {
        _prob_crossover = prob;
        _prob_mutation = 1 - prob;
        _prob_reproduction = 0;
    } else {
        _prob_crossover = prob;
        _prob_reproduction = 1 - prob - _prob_mutation;
    }
}

double params_t::prob_crossover() const
{
    return _prob_crossover;
}

void params_t::prob_reproduction(double prob)
{
    if (prob > 1) throw std::invalid_argument{"reproduction probability > 1 (" +
                          std::to_string(prob) + ")."};

    // HYPOTHESIS scale the values. WORKS!!1: 1) calculates the proportion of
    // each prob and 2) scale for the new possible probability (the complement
    // of prob).
    double total = _prob_crossover + _prob_mutation;
    _prob_mutation = (1 - prob) * (_prob_mutation / total);
    _prob_crossover = (1 - prob) * (_prob_crossover / total);
    _prob_reproduction = prob;
}

double params_t::prob_reproduction() const
{
    return _prob_reproduction;
}

void params_t::selection_method(sel_t method)
{
    _sel = method;
}

sel_t params_t::selection_method() const
{
    return _sel;
}

void params_t::tournament(std::size_t k)
{
    _k = k;
}

std::size_t params_t::tournament() const
{
    return _k;
}

void params_t::generation_method(gen_t method)
{
    _gen = method;
}

gen_t params_t::generation_method() const
{
    return _gen;
}

void params_t::error_metric(err_t metric)
{
    _err = metric;
}

err_t params_t::error_metric() const
{
    return _err;
}

void params_t::eletism(bool enable)
{
    _eletism = enable;
}

bool params_t::eletism() const
{
    return _eletism;
}
    

void params_t::reset()
{
    _population_sz     = 500;
    _max_depth         = 7;
    _max_generation    = 100;
    _threshold         = 0.001;
    _prob_crossover    = 0.8;
    _prob_mutation     = 0.15;
    _prob_op_mutation  = 0.05;
    _prob_reproduction = 0.05;
    _sel               = sel_t::tour_t;
    _k                 = 5;
    _gen               = gen_t::ramped_t;
    _err               = err_t::mae_t;
    _eletism           = false;
}

std::string params_t::selec_met_str() const
{
    using namespace std;

    switch (_sel) {
    case sel_t::roul_t:
        return "roulette wheel"; break;
    case sel_t::tour_t:
        return to_string(_k) + "-tournament"; break;
    }

    return "";
}

std::string params_t::gen_met_str() const
{
    switch (_gen) {
    case gen_t::full_t:
        return "full"; break;
    case gen_t::grow_t:
        return "grow"; break;
    case gen_t::ramped_t:
        return "ramped half-and-half"; break;
    default:
        throw std::invalid_argument{"params_t::gen_met_str"};
    }

}

std::string params_t::err_met_str() const
{
    switch (_err) {
    case err_t::mae_t:
        return "mae"; break;
    case err_t::mse_t:
        return "mse"; break;
    case err_t::rmse_t:
        return "rmse"; break;
    default:
        throw std::invalid_argument{"params_t::err_met_str"};
    }
}

// std::string params_t::seeds_str() const
// {
//     std::string out{};

    
//     for (int i = 0; i < _seeds.size(); ++i)
//         out += std::to_string(_seeds[i]) + ", ";

//     out += std::to_string(_seeds.back());

//     return out;
// }

gp_operators_t::gp_operators_t(random_t& random,
                               individual_handler_t& handler)
    :
    _rd{random},
    _hdl{handler}
{
}

double gp_operators_t::fitness(const individual_t& individual,
                               const std::vector<entry_t>& input_data,
                               err_t err_m)
{
    std::vector<double> error{};

    for (auto xy: input_data) {
        vars_t x = xy.first;
        double y = xy.second;
        double y_pred = _hdl.eval(individual, x);
        double dy = y - y_pred;

        //        if (std::isnan(y_pred)) throw std::logic_error{"NAN"};
        error.push_back(err_m == err_t::mae_t ? fabs(dy) : dy * dy);
    }

    auto b = error.begin(), e = error.end();
    double me = std::accumulate(b, e, 0.0)/error.size(); // mean error.
    
    return (err_m == err_t::rmse_t) ? sqrt(me) : me; 
}

void gp_operators_t::population_fitness(individuals_t& population,
                                        const std::vector<entry_t>& input_data,
                                        err_t err_m)
{
    for (auto& individual : population) {
        double error = fitness(individual.first, input_data, err_m);
        individual.second = error;
    }
}

individual_t gp_operators_t::full_gen(std::size_t max_depth)
{
    individual_t ind{max_depth};
    if (max_depth <= 1) {
        ind[0] = _hdl.rd_term();
        return ind;
    }    

    for (std::size_t i = 0; i < (std::size_t)pow(2, max_depth - 2) - 1; ++i)
        ind[i] = _hdl.rd_oper();
    
    for (std::size_t i = pow(2, max_depth - 2) - 1;
                     i < pow(2, max_depth - 1) - 1; ++i) {
        if (_rd.choose()) {
            ind[i] = _hdl.rd_func();
            ind[ind.lchild(i)] = _hdl.rd_term();
        } else {
            ind[i] = _hdl.rd_oper();
            ind[ind.lchild(i)] = _hdl.rd_term();
            ind[ind.rchild(i)] = _hdl.rd_term();
        }
    }

    return ind;
}

individual_t gp_operators_t::grow_gen(std::size_t max_depth)
{
    individual_t ind{max_depth};
    grow_gen_recursive(ind, max_depth, 0);
    return ind;
}

individual_t gp_operators_t::gen_individual(gen_t gen, std::size_t max_depth)
{
    if (gen == gen_t::full_t)
        return full_gen(max_depth);
    else if (gen == gen_t::grow_t)
        return grow_gen(max_depth);
    else if (gen == gen_t::ramped_t) 
        return gen_individual(_rd.choose() ?
                              gen_t::full_t : gen_t::grow_t, max_depth);
    else
        throw std::invalid_argument{"gp_operators_t::gen_individual"};
}

individual_t gp_operators_t::selection(const individuals_t& population,
                                       sel_t sel,
                                       std::size_t k)
{
    if (sel == sel_t::roul_t)
        return selection_rw(population);
    else if (sel == sel_t::tour_t)
        return selection_t(population, k);
    else
        throw std::invalid_argument{"gp_operators_t::selection"};
}
           
individual_t gp_operators_t::selection_rw(const individuals_t& population)
{
    auto b = population.begin(), e = population.end();
    auto sum_callable = [](auto acc, auto v) { return acc + v.second; };
    double total_fitness = std::accumulate(b, e, 0.0, sum_callable);

    double prob = _rd.prob();
    
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
    using Tp = typename individuals_t::value_type;
    auto selected = _rd.sample<Tp>(population, k);
    return selected[0].first;
}

individual_t gp_operators_t::crossover(const individual_t& p1,
                                       const individual_t& p2)
{
    if (p1.max_depth() != p2.max_depth())
        throw std::invalid_argument{"gp_operators::crossover"};

    std::size_t p1_pt = 0, p2_pt = 0, total_sz = 0;
    do {
        p1_pt = _hdl.rd_point(p1); p2_pt = _hdl.rd_point(p2);
        total_sz = p1_pt + (p2.max_depth() - p2_pt);
    } while (total_sz > p1.max_depth());

    if (p1.has_parent(p1_pt)) {
        individual_t n_ind = p1;
        n_ind.clear_subtree(p1_pt);
        n_ind.copy_subtree(p1_pt, p2, p2_pt);
        return n_ind;
    } else {
        individual_t n_ind{p1.max_depth()};
        n_ind.copy_subtree(0, p2, p2_pt);
        return n_ind;
    }
}

individual_t gp_operators_t::mutation_op(const individual_t& ind,
                                         double prob)
{
    individual_t n_ind = ind;
    for (std::size_t p = 0; p < ind.size() ; ++p) {

        if (_rd.prob() > prob) continue;
        
        switch (n_ind[p]._code) {
        case code_t::oper:
            n_ind[p] = _hdl.rd_oper(); break;
        case code_t::func:
            n_ind[p] = _hdl.rd_func(); break;
        case code_t::var:
            n_ind[p] = _hdl.rd_var(); break;
        case code_t::cons:
            n_ind[p] = _hdl.rd_cons(); break;
        }
    }
    
    return n_ind;
}

individual_t gp_operators_t::mutation_sb(const individual_t& ind)
{
    std::size_t pt = _hdl.rd_point(ind);
    std::size_t pt_depth = individual_t::depth(pt);
    std::size_t tree_depth = ind.max_depth();

    individual_t n_ind = ind;
    individual_t sub_tree = full_gen(tree_depth - pt_depth); // to replace.

    n_ind.clear_subtree(pt); // clear old data.
    n_ind.copy_subtree(pt, sub_tree, 0);  // copy new tree.

    return n_ind;
}

std::size_t gp_operators_t::count_unique(const individuals_t& population)
{
    std::set<double> counter{};
    for (auto individual: population)
        counter.insert(individual.second);

    return counter.size();
}

void gp_operators_t::grow_gen_recursive(individual_t& ind,
                                        std::size_t max_depth,
                                        std::size_t pt)
{
    if (max_depth - 1 == individual_t::depth(pt)) { // last depth, gen terminals.
        ind[pt] = _hdl.rd_term();
    } else {
        if (_rd.choose()) { // gen functions or operators.
            ind[pt] = _hdl.rd_func();
            if (_rd.choose()) // expand or stop ?
                ind[ind.lchild(pt)] = _hdl.rd_term();
            else
                grow_gen_recursive(ind, max_depth, ind.lchild(pt));
        } else {
            ind[pt] = _hdl.rd_oper();
            if (_rd.choose()) // expand or stop ?
                ind[ind.lchild(pt)] = _hdl.rd_term();
            else
                grow_gen_recursive(ind, max_depth, ind.lchild(pt));
            if (_rd.choose())
                ind[ind.rchild(pt)] = _hdl.rd_term();
            else
                grow_gen_recursive(ind, max_depth, ind.rchild(pt));
        }
    }
}

symbolic_regression_t::symbolic_regression_t(params_t params,
                                             random_t& rd,
                                             individual_handler_t& hdl)
    :
    _params{params},
    _rd{rd},
    _hdl{hdl},
    _gpo{rd, hdl}
{
}

symbolic_regression_t::~symbolic_regression_t()
{
}

const state_t& symbolic_regression_t::state()
{
    return _states.back().back();
}

void symbolic_regression_t::params(params_t params)
{
    _params = params;
}

params_t symbolic_regression_t::params() const
{
    return _params;
}

std::string symbolic_regression_t::report() const
{
    std::string report_text{};
    report_text = static_cast<std::string>(_params) +  "\n";

    //for (auto state: _states)
    //    report_text += static_cast<std::string>(state) + "\n";
        
    return report_text;
}

void symbolic_regression_t::train(const std::vector<entry_t>& t_set)
{
    for (auto seed: _params._seeds) {
        _rd.use_seed(seed);
        initialize_state();
        initialize_population(t_set);
        update_state();
    
        while (!check_stop_condition()) {
            cout << "new gen." << endl;
            next_generation(t_set);
            update_state();
        }
        cout << "finished" << endl;
        cout << _population[0].second << endl;
    }
}

double symbolic_regression_t::predict(const vars_t& vars) const
{
    individual_t first = _population[0].first;
    return _hdl.eval(first, vars);
}

void symbolic_regression_t::clear_line() const
{
    std::cout << "\r" << std::flush;
    for (int i = 0; i < 120; ++i)
        std::cout << " ";
    std::cout << "\r" << std::flush;
}

void symbolic_regression_t::print_progress() const
{
    if (_states.size() == 0) return;

    state_t cur_state = _states.back().back();
    std::cout << std::setprecision(3);
    std::cout << cur_state.generation() << " "
              << cur_state.population_sz() << " "
              << cur_state.unique_individuals() << " "
              << cur_state.better_than_previous() << " "
              << cur_state.min() << " "
              << cur_state.max() << " "
              << cur_state.mean() << " "
              << cur_state.stddev() << " "
              << cur_state.var() << " "
              << cur_state.median() << std::endl;
             
}

void symbolic_regression_t::initialize_state()
{
    _states.emplace_back(std::vector<state_t>{});
}

void symbolic_regression_t::update_state()
{
    std::vector<double> pop_fit{};
    for (auto individual: _population)
        pop_fit.push_back(individual.second);

    state_t state{pop_fit};    
    state._generation = _states.back().size();
    state._population_sz = _population.size();
    state._unique_individuals = _gpo.count_unique(_population);
    _states.back().emplace_back( state);
}

void symbolic_regression_t::initialize_population(
                               const std::vector<entry_t>& t_set)
{
    auto err = _params._err;
    auto gen = _params._gen;    
    auto max_depth = _params._max_depth;

    for (std::size_t i = 0; i < _params._population_sz; ++i) {
        auto individual = _gpo.gen_individual(gen, max_depth);
        auto fitness = _gpo.fitness(individual, t_set, err);
        _population.emplace_back(individual, fitness);
    }

    auto b = _population.begin(), e = _population.end();
    std::sort(b, e, [](auto a, auto b) {return a.second < b.second;});
}

void symbolic_regression_t::next_generation(const std::vector<entry_t>& t_set)
{
    auto pop_sz = _params._population_sz;
    
    double pc = _params._prob_crossover, pm = _params._prob_mutation,
           pr = _params._prob_reproduction;
    
    std::size_t c_sz = pc * pop_sz, m_sz = pm * pop_sz,
        r_sz = pop_sz - c_sz - m_sz;
    
    individuals_t new_pop{};
    do_crossover(new_pop, c_sz, t_set);
    do_mutation(new_pop, m_sz, t_set);
    do_reproduction(new_pop, r_sz, t_set);

    using Tp = typename individuals_t::value_type;
    auto cmp_less = [](Tp a, Tp b) { return a.second < b.second; };
    std::sort(new_pop.begin(), new_pop.end(), cmp_less);

    if (_params._eletism) {
        new_pop.insert(new_pop.end(), _population.begin(), _population.end());

        auto b = new_pop.begin(), e = new_pop.end();
        std::inplace_merge(b, b + pop_sz, e, cmp_less);

        b = new_pop.begin(), e = new_pop.end();
        new_pop.erase(b + pop_sz, e);
    }

    _population = new_pop;
}

void symbolic_regression_t::do_crossover(individuals_t& individuals,
                                         std::size_t n,
                                         const std::vector<entry_t>& t_set)
{
    std::size_t max_depth = _params._max_depth;    
    sel_t sel = _params._sel;
    std::size_t k = _params._k;
    err_t err = _params._err;

    for (std::size_t i = 0; i < n; ++i) {
        auto p1 = _gpo.selection(_population, sel, k);
        auto p2 = _gpo.selection(_population, sel, k);
        individual_t offspring = _gpo.crossover(p1, p2);
        double fitness = _gpo.fitness(offspring, t_set, err);
        individuals.emplace_back(offspring, fitness);
    }
}

void symbolic_regression_t::do_mutation(individuals_t& individuals,
                                        std::size_t n,
                                        const std::vector<entry_t>& t_set)
{
    auto max_depth = _params._max_depth;
    auto sel = _params._sel;
    auto k = _params._k;
    auto err = _params._err;
    auto prob_op_mut = _params._prob_op_mutation;
    
    std::size_t op_sz = std::floor(0.5 * n), sb_sz = n - op_sz;

    for (std::size_t i = 0; i < op_sz; ++i) {
        individual_t p = _gpo.selection(_population, sel, k);
        individual_t offspring = _gpo.mutation_op(p, prob_op_mut);
        double fitness = _gpo.fitness(offspring, t_set, err);
        individuals.emplace_back(offspring, fitness);
    }

    for (std::size_t i = 0; i < sb_sz; ++i) {
        individual_t p = _gpo.selection(_population, sel, k);
        individual_t offspring = _gpo.mutation_sb(p);
        double fitness = _gpo.fitness(offspring, t_set, err);
        individuals.emplace_back(offspring, fitness);
    }
}

void symbolic_regression_t::do_reproduction(individuals_t& individuals,
                                            std::size_t n,
                                            const std::vector<entry_t>& t_set)
{
    auto sel = _params._sel;
    auto k = _params._k;
    auto err = _params._err;
    
    for (std::size_t i = 0; i < n; ++i) {
        individual_t offspring = _gpo.selection(_population, sel, k);
        double fitness = _gpo.fitness(offspring, t_set, err);
        individuals.emplace_back(offspring, fitness);
    }
}

bool symbolic_regression_t::check_stop_condition() const
{
    return _params._max_generation <= _states.back().size();
}

#include <iostream>
#include <cmath>

using namespace std;

int main()
{
    std::vector<int> seeds{1, 2, 3};
    random_t rd{seeds};
    individual_handler_t hdl{rd};
    hdl.add_function("sin", [](double x) { return sin(x); });
    hdl.add_function("cos", [](double x) { return cos(x); });
    hdl.add_function("tan", [](double x) { return tan(x); });
    hdl.add_function("sinh", [](double x) { return sinh(x); });
    hdl.add_function("cosh", [](double x) { return cosh(x); });
    hdl.add_function("tanh", [](double x) { return tanh(x); });
    hdl.add_function("log", [](double x) { return log(abs(x)+1); });
    hdl.add_function("log10", [](double x) { return log10(abs(x)+1); });
    hdl.add_operator("+", 1, [](double l, double r) { return l + r; });
    hdl.add_operator("-", 2, [](double l, double r) { return l - r; });
    hdl.add_operator("*", 3, [](double l, double r) { return l * r; });
    auto safe_div = [](double l, double r) { return r == 0 ? l/0.001 : l / r; };
    hdl.add_operator("/", 4, safe_div);
    hdl.add_variable(1);

    std::vector<entry_t> dataset{};
    dataset.emplace_back(vars_t{0.0}, safe_div(1, 0.0));
    dataset.emplace_back(vars_t{0.1}, safe_div(1, 0.1));
    dataset.emplace_back(vars_t{0.2}, safe_div(1, 0.2));
    dataset.emplace_back(vars_t{0.3}, safe_div(1, 0.3));
    dataset.emplace_back(vars_t{0.4}, safe_div(1, 0.4));
    dataset.emplace_back(vars_t{0.5}, safe_div(1, 0.5));
    dataset.emplace_back(vars_t{0.6}, safe_div(1, 0.6));


    
    params_t params{{seeds}};
    params.max_generation(60);
    symbolic_regression_t sreg{params, rd, hdl};
    sreg.train(dataset);
    cout << hdl.str(sreg._population[0].first) << endl;

    return 0;
}
