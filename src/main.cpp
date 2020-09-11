#include <iostream>
#include <string>
#include <cctype>
#include <functional>
#include <fstream>

#include "common_types.hpp"
#include "symbolic_regression.hpp"
#include "load_data.hpp"
#include "argh.h"

extern bool help_message(argh::parser& cmdl);
extern std::string load_input(argh::parser& cmdl, bool verbose);
extern std::string output_file(argh::parser& cmdl, bool verbose);
extern std::size_t load_n_vars(argh::parser& cmdl, bool verbose);
extern std::size_t load_pop_sz(argh::parser& cmdl, bool verbose);
extern std::size_t load_max_pop(argh::parser& cmdl, bool verbose);
extern gen_t load_gen_method(argh::parser& cmdl, bool verbose);
extern std::size_t load_max_depth(argh::parser& cmdl, bool verbose);
extern sel_t load_sel_method(argh::parser& cmdl, bool verbose);
extern std::size_t load_k_tournament(argh::parser& cmdl, bool verbose);
extern double load_prob_mutation(argh::parser& cmdl, bool verbose);
extern double load_prob_crossover(argh::parser& cmdl, bool verbose);
extern err_t load_error_metric(argh::parser& cmdl, bool verbose);
extern double load_fitness_threshold(argh::parser& cmdl, bool verbose);
extern bool load_eletism(argh::parser& cmdl, bool verbose);
extern std::vector<std::vector<int>> load_seed_data(argh::parser& cmdl,
                                                    bool verbose);

extern double stan(double x);
extern double ssqrt(double x);
extern double sinv(double x);
extern double slog(double x);
extern double slog10(double x);
extern double scosh(double x);
extern double ssinh(double x);
extern double sdiv(double l, double r);
extern double spow(double l, double r);

int main(int argc, char **argv)
{
    auto cmdl = argh::parser(argc, argv,
                             argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);

    bool verbose = true;
    if (help_message(cmdl)) return 0;
    std::string filename = load_input(cmdl, verbose);        
    std::string outfile = output_file(cmdl, verbose);
    std::size_t n_vars = load_n_vars(cmdl, verbose);
    std::size_t population_sz = load_pop_sz(cmdl, verbose);
    std::size_t max_generation = load_max_pop(cmdl, verbose);
    gen_t generation_method = load_gen_method(cmdl, verbose);
    std::size_t max_depth = load_max_depth(cmdl, verbose);
    sel_t selection_method = load_sel_method(cmdl, verbose);
    std::size_t k = load_k_tournament(cmdl, verbose);  
    double prob_mutation = load_prob_mutation(cmdl, verbose);
    double prob_crossover = load_prob_crossover(cmdl, verbose);
    err_t error_metric = load_error_metric(cmdl, verbose);
    double fitness_threshold = load_fitness_threshold(cmdl, verbose);
    bool eletism = load_eletism(cmdl, verbose);
    std::vector<std::vector<int>> seeds = load_seed_data(cmdl, verbose);
    
    //    auto data = load_data(filename);
    // for (auto row: data) {
    //     for (auto x: row.first)
    //         std::cout << x << ", ";
    //     std::cout << "= " << row.second << std::endl;
    // }

    // for (auto seed: seeds) {
    //     for (auto s: seed)
    //         std::cout << s << ", ";
    //     std::cout << std::endl;
    // }
    
    params_t params{seeds};
    params.selection_method(selection_method);
    params.generation_method(generation_method);
    params.error_metric(error_metric);
    params.tournament(k);
    params.population_sz(population_sz);
    params.max_depth(max_depth);
    params.prob_mutation(prob_mutation);
    params.prob_crossover(prob_crossover);
    params.eletism(eletism);
    params.threshold(fitness_threshold);
    params.max_generation(max_generation);

    random_t rd{{0}};
    individual_handler_t hdl{rd};

    hdl.add_function("sin", std::sin);
    hdl.add_function("cos", std::cos);
    hdl.add_function("tan", stan);
    hdl.add_function("sinh", scosh);
    hdl.add_function("cosh", ssinh);
    hdl.add_function("tanh", std::tanh);
    hdl.add_function("log", slog);
    hdl.add_function("log10", slog10);
    hdl.add_function("inv", sinv);
    hdl.add_function("sqrt", ssqrt);

    hdl.add_operator("+", 1, [](double l, double r) { return l + r; });
    hdl.add_operator("-", 2, [](double l, double r) { return l - r; });
    hdl.add_operator("*", 3, [](double l, double r) { return l * r; });
    hdl.add_operator("/", 4, sdiv);
    hdl.add_operator("^", 5, spow);
    
    hdl.add_variable(n_vars);

    symbolic_regression_t sreg{params, rd, hdl};
    auto dataset = random_t::sample(load_data(filename), 500);
    sreg.train(dataset);

    std::string report = sreg.report();    

    if (!outfile.empty()) {
        std::ofstream ofs{outfile};
        ofs << report;
        std::cout << "report saved in: " << outfile << "!" << std::endl;
    } else {
        std::cout << report << std::endl;
    }
    
    return 0;
}
