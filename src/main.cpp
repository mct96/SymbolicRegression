#include <iostream>
#include <fstream>
#include <string>

#include "symbolic_regression.hpp"
#include "load_data.hpp"

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

int main(int argc, char **argv)
{
    if (argc <= 1 && 0) {
        usage_message();
        return 1;
    }

    training_set_t training{};

    for (double x = -1; x < 1; x += 0.15)  {
        for (double y = -1; y < 1; y += 0.15) {
            training.emplace_back(std::vector<double>{x, y},
                                  cos(x)*x + sin(y)/2 + x * y);
        }
    }
       
    auto data = sr::load_data("../../data/concrete.txt");

    for (int i = 0; i < data.size(); ++i)
        if (i % 5 == 0)
            training.push_back(data[i]);
        
    parameters_t params{8};
    params.selection_method(selection_method_t::tournament);
    params.error_metric(error_metric_t::mae);
    params.tournament(12);
    params.max_depth(7);
    params.prob_one_point_mutation(0.2);
    params.eletism(false);
    
    symbolic_regression_t sm{params, training};
    sm.initialize_population();
    for (auto g = 0; g < 50; ++g)
        sm.next_generation();
    sm.report();
    
    return 0;
}
