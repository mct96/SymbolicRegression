#include <iostream>
#include <fstream>
#include <string>
#include <cctype>

#include "symbolic_regression.hpp"
#include "load_data.hpp"
#include "argh.h"

std::string to_lower(std::string str)
{
    std::transform(str.begin(), str.end(), str.begin(),
                       [&](auto c) { return std::tolower(c); });
    return str;
}

void usage_message() {
    std::ifstream ofs{"../docs/usage.txt"};
    std::string text{};
    
    while (ofs && !ofs.eof()) {
        char buffer[120] = "";
        ofs.read(buffer, 120);
        text += buffer;
    }
    std::cout << text << std::endl;
}

bool is_csv(std::string filename)
{
    auto ext = filename.substr(filename.size()-3);
    return ext == std::string{"csv"};
}

int main(int argc, char **argv)
{
    if (argc <= 1 && 0) {
        usage_message();
        return 1;
    }

    argh::parser cmdl{argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION};
    if (cmdl[{"-h", "--help"}]) {
        usage_message();
        return 0;
    }

    std::string filename = cmdl[1];
    if (!is_csv(filename)) {
        std::cerr << "Invalid file. Try a *.csv" << std::endl;
        return 1;
    }
    
    auto data = sr::load_data(filename); // loading dataset.
    training_set_t training = random_t{{1}}.sample(data, (int)data.size()/2);

    std::string outfile{};
    if (cmdl({"-o", "--output"}) >> outfile) {
        std::cout << "output-results: " << outfile << std::endl;
    } else {
        std::cout << "output-results: " << "NO REPORT" << std::endl;
    }

    std::size_t n_vars = 1;
    if ((cmdl({"-n", "--nvar"}, 1) >> n_vars)) {
        std::cout << "number of variables: " << n_vars << std::endl;
    } else {
        std::cout << "number of variables: " << n_vars;
    }

    int population_sz = 0;
    if (cmdl({"-p", "--population-size"}, 1000) >> population_sz) {
        std::cout << "population size: " << population_sz << std::endl;
    }

    int max_generation = 0;
    if (cmdl({"-G", "--generation-max"}, 50) >> max_generation) {
        std::cout << "maximum generation: " << max_generation << std::endl;
    }

    generation_method_t generation_method{};
    std::string g_met{};
    if (cmdl({"-g", "--generation-method"}, "full") >> g_met) {
        g_met = to_lower(g_met);
        if (g_met == "full")
            generation_method = generation_method_t::full;
        if (g_met == "grow")
            generation_method = generation_method_t::grow;
        if (g_met == "ramped_hh")
            generation_method = generation_method_t::ramped_hh;
        std::cout << "generation method: " << g_met << std::endl;
    }
    
    int max_depth = 5;
    if (cmdl({"-d", "--max-depth"}, 5) >> max_depth) {
        std::cout << "max-depth: " << max_depth << std::endl;
    }

    selection_method_t selection_method = selection_method_t::tournament;
    std::string s_met{};
    if (cmdl({"-s", "--selection-method"}, "tournament") >> s_met) {
        s_met = to_lower(s_met);
        if (s_met == "roulette")
            selection_method = selection_method_t::roulette_wheel;
        if (s_met == "tournament")
            selection_method = selection_method_t::tournament;
        std::cout << "selection method: " << s_met << std::endl;
    }

    int k = 2;
    if (cmdl({"-k", "--k"}, 2) >> k) {
        std::cout << "k-factor: " << k << std::endl;
    }

    double prob_mutation = 0.2;
    if (cmdl({"-M", "--prob-mutation"}, 0.2) >> prob_mutation) {
        std::cout << "mutation probability: " << prob_mutation << std::endl;
    }

    double prob_crossover = 0.8;
    if (cmdl({"-c", "--prob-crossover"}, 0.8) >> prob_crossover) {
        std::cout << "crossover probability: " << prob_crossover << std::endl;
    }
    
    error_metric_t error_metric; std::string e_metric;
    if (cmdl({"-e", "--error_metric"}, "mae") >> e_metric) {
        e_metric = to_lower(e_metric);
        if (e_metric == "mae") error_metric = error_metric_t::mae;
        if (e_metric == "mse") error_metric = error_metric_t::mse;            
        if (e_metric == "rmse") error_metric = error_metric_t::rmse;
        std::cout << "error metric: " << e_metric << std::endl;
    }

    double fitness_threshold = 0.2;
    if (cmdl({"-f", "--fitness-threshold"}, 0.3) >> fitness_threshold) {
        std::cout << "fitness threshold: " << fitness_threshold << std::endl;
    }

    bool eletism = false;
    eletism = (cmdl[{"-E", "--eletism"}]);
    std::cout << "eletism: " << (eletism ? "true" : "false") << std::endl;

    std::string seeds_file = "";
    std::vector<std::vector<int>> seeds{};
    if (cmdl({"-t", "--test"}) >> seeds_file) {
        seeds = sr::load_seeds(seeds_file);
        if (seeds.size() > 0)  {
            std::cout << "seeds loaded: " << seeds.size() << std::endl;
        } else {
            std::cout << "no seeds found in " << seeds_file << "!" << std::endl;
            return 1;
        }            
    } else {
        seeds.push_back({18, -1, -6, 12, -19, 15, 9, -13, -15, -5, 7, 17});
    }
    

    parameters_t params{n_vars, seeds[0]};
    
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

    // NOTE
    // If you are using Emacs, I created a script named "random_generator.el"
    // to generate a list of random numbers.

    int i = 1;
    for (auto seed: seeds) {
        std::string report{};
        params.seeds(seed);
        symbolic_regression_t sr{params, training};
        sr.train();
        report += sr.report();
    

        if (!outfile.empty()) {
            std::ofstream ofs{outfile + "_"+ std::to_string(i)};
            ofs << report;
            std::cout << "report saved in: " << outfile << "!" << std::endl;
        } else {
            std::cout << report << std::endl;
        }
        
        i++;
    }
    return 0;
}
