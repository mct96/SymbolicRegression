#include <iostream>
#include <fstream>
#include <string>

#include "symbolic_regression.hpp"
#include "load_data.hpp"
#include "argh.h"

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

    argh::parser cmdl(argv);
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
    training_set_t training = random_t{{1}}.sample(data, (int)data.size()/10);

    std::string outfile{};
    if (cmdl({"-o", "--output"}) >> outfile) {
        std::cout << "output-results: " << outfile << std::endl;
    } else {
        std::cout << "output-results: " << "NO REPORT" << std::endl;
    }

    std::size_t n_vars = 1;
    if (cmdl({"-n", "-nvar"}, 1) >> n_vars) {
        std::cout << "number of variables: " << n_vars << std::endl;
    }

    int population_sz = 0;
    if (cmdl({"-p", "--population-size"}, 1000) >> population_sz) {
        std::cout << "population size: " << population_sz << std::endl;
    }

    int max_generation = 0;
    if (cmdl({"-g", "--generation-max"}, 50) >> max_generation) {
        std::cout << "maximum generation: " << max_generation << std::endl;
    }
 
    int max_depth = 5;
    if (cmdl({"-d", "--max-depth"}, 5) >> max_depth) {
        std::cout << "max-depth: " << max_depth << std::endl;
    }

    int selection_method = 0;
    if (cmdl({"-s", "--selection-method"}, 1) >> selection_method) {
        std::cout << "selection method: ";
        auto s = (selection_method_t)selection_method;
        if (s == selection_method_t::roulette_wheel) std::cout  << "roulette";
        if (s == selection_method_t::tournament) std::cout << "tournament";
        std::cout << std::endl;
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
    
    int error_metric = 0;
    if (cmdl({"-e", "--error_metric"}, 1) >> error_metric) {
        std::cout << "error metric: ";
        auto e = (error_metric_t)error_metric;
        if (e == error_metric_t::mae) std::cout << "MAE";
        if (e == error_metric_t::mse) std::cout << "MSE";            
        if (e == error_metric_t::rmse) std::cout << "RMSE";
        std::cout << std::endl;
    }

    double fitness_threshold = 0.2;
    if (cmdl({"-f", "--fitness-threshold"}, 0.2) >> fitness_threshold) {
        std::cout << "fitness threshold: " << fitness_threshold << std::endl;
    }

    bool eletism = false;
    if (cmdl({"-E", "--eletism"}, false) >> eletism) {
        std::cout << "eletism: " << (eletism ? "true" : "false") << std::endl;
    }


    parameters_t params{n_vars, {1, 2, 6, 3}};
    params.selection_method((selection_method_t)selection_method);
    params.error_metric((error_metric_t)error_metric);
    params.tournament(k);
    params.population_sz(population_sz);
    params.max_depth(max_depth);
    params.prob_mutation(prob_mutation);
    params.prob_mutation(prob_crossover);
    params.eletism(eletism);
    params.threshold(fitness_threshold);
    
    symbolic_regression_t sr{params, training};



    std::cout << "Training started." << std::endl;
    sr.report();
    
    return 0;
}
