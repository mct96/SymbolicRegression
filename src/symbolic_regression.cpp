#include "symbolic_regression.hpp"

#include <cmath>
#include <stdexcept>

//namespace sr
//{

double eval_function(gene_t func, double x)
{
    auto code = func._code;
 
    double y = 0;
    switch (code) {
    case 0: y = sin(x); break;
    case 1: y = cos(x); break;        
    case 2: y = tan(x); break;
    case 3: y = sinh(x); break;
    case 4: y = cosh(x); break;
    case 5: y = tanh(x); break;
    case 6: y = log10(x); break;
    case 7: y = log(x); break;
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
    case 3: result = lop / rop; break;
    case 4: result = pow(lop, rop); break;
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
    case 3: out << "sinh" ; break;
    case 4: out << "cosh" ; break;
    case 5: out << "tanh" ; break;
    case 6: out << "log10"; break;
    case 7: out << "log"  ; break;
    default:
        throw std::invalid_argument{"print function"};
    }

    out << "("; print(out, u, lchild(pos)); out << ")";
}

void print_oper(std::ostream& out, const individual_t& u, std::size_t pos)
{
    gene_t oper = u[pos], loper = u[lchild(pos)], roper = u[rchild(pos)];
    
    // print left operand. surround lower precedence operator with parentheses.
    if (loper._class == class_t::operator_t && loper._code < oper._code) {
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
    case 4: out << " ^ "; break;
    default:
        throw std::invalid_argument{"eval operator"};
    }

    // print right operand. surround lower precedence operator with parentheses.
    if (roper._class == class_t::operator_t && roper._code < oper._code) {
        out << "("; print(out, u, rchild(pos)); out << ")";
    } else {
        print(out, u, rchild(pos));
    }
}

void print_var(std::ostream& out, const individual_t& u, std::size_t pos)
{
    gene_t var = u[pos];
    out << "x_" << static_cast<int>(var._value);
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

double gp_operators_t::fitness(const individual_t& individual,
                               const data_t& variables,
                               double target)
{
    double predicted = eval(individual, 0, variables);
    double error = abs(target - predicted);
    double rel = error / target;
    return max(0, 1 - rel * rel);
}

double gp_operators_t::population_error(const individuals_t& population,
                                        const data_t& variables,
                                        double target,
                                        error_metric_t error_metric)
{
    std::vector<double> individual_error(population.size());
    for (const auto& individual: population) {
        auto predicted = eval(individual, 0, variables);
        individual_error.push_back(target - predicted);
    }

    auto b = individual_error.begin(), e = individual_error.end();
    switch (error_metric) {
    case error_metric_t::mae:
        std::for_each(b, e, abs);
        break;
    case error_metric_t::mse:
    case error_metric_t::rmse:
        std::for_each(b, e, [](double e) { return e * e; });
        break;
    }
    
    double avg_error = std::accumulate(b, e, 0) / population.size();
    auto is_rsme = error_metric == error_metric_t::rmse;
    return is_rmse ? sqrt(avg_error) : avg_error;
}

individual_t gp_operators_t::full_gen(std::size_t N)
{
}

individual_t gp_operators_t::grow_gen(std::size_t N)
{
}

individual_t gp_operators_t::gen_individual(generation_method_t gm,
                                            std::size_t N)
{
    if (gm == generation_method_t::full)
        return full_gen(N);
    else if (gm == generation_method_t::grow)
        return grow_gen(N);
}

individuals_t gp_operators_t::selection_rw(const individuals_t& population)
{
}

individuals_t gp_operator_t::selection_t(const individuals_t& population,
                                         std::size_t sz)
{
}

individuals_t gp_operator_t::crossover(const individual_t& p1,
                                       const individual_t& p2)
{
    
}

void gp_operator_t::mutation_op(individual_t& individual, double p)
{
}

void gp_operator_t::mutation_ex(individual_t& individual, double p)
{
}

void gp_operator_t::mutation_rd(individual_t& individual, double p)
{
}

void gp_operator_t::copy_subtree(const individual_t& src,
                                 individual_t& dst,
                                 std::size_t sp)
{
    auto maxsz = dst.size();    // longest child.
    auto buffer = std::vector<std::size_t>{sp}; // to do list.

    while (buffer.size() > 0) {
        auto point = buffer.pop_back();
        dst[point] = src[point]; // copy.

        auto lc = lchild(point), rc = rchild(point);
        if (lc < maxsz)         // if lchild exists.
            buffer.push_back(lc);
        
        if (rc < maxdz)         // if rchild exists.
            buffer.push_back(rc);
    }
}

symbolic_regression_t::symbolic_regression_t(parameters_t params)
{
}
    
symbolic_regression_t::~symbolic_regression_t()
{
}

    
    
//}

#include <iostream>

using namespace std;
int main() {
    individual_t individual(15);
    individual[ 0]._class = class_t::operator_t;
    individual[ 0]._code  = 0;
    individual[ 1]._class = class_t::operator_t;
    individual[ 1]._code  = 4;
    individual[ 2]._class = class_t::operator_t;
    individual[ 2]._code  = 4;
    individual[ 3]._class = class_t::function_t;
    individual[ 3]._code  = 0;
    individual[ 4]._class = class_t::constant_t;
    individual[ 4]._value = 2.0;
    individual[ 5]._class = class_t::function_t;
    individual[ 5]._code  = 1;
    individual[ 6]._class = class_t::constant_t;
    individual[ 6]._value = 2.0;
    individual[ 7]._class = class_t::variable_t;
    individual[ 7]._code  = 0;
    individual[11]._class = class_t::variable_t;
    individual[11]._code  = 0;

    double x = eval(individual, 0, {1.2});
    cout << individual << " = " << x << endl;

    individual[ 0]._class = class_t::operator_t;
    individual[ 0]._code  = 2;
    individual[ 1]._class = class_t::operator_t;
    individual[ 1]._code  = 0;
    individual[ 2]._class = class_t::operator_t;
    individual[ 2]._code  = 0;
    individual[ 3]._class = class_t::constant_t;
    individual[ 3]._value = 3.0;
    individual[ 4]._class = class_t::constant_t;
    individual[ 4]._value = 2.0;
    individual[ 5]._class = class_t::constant_t;
    individual[ 5]._value = 2.0;
    individual[ 6]._class = class_t::constant_t;
    individual[ 6]._value = 1.0;

    x = eval(individual, 0, data_t{});
    cout << individual << " = " << x << endl;
    return 0;
}
