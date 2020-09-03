#include "load_data.hpp"

#include <chrono>
#include <iostream>
#include "rapidcsv.h"

namespace sr
{

std::vector<Xy> load_data(std::string filename)
{
    using namespace std;
    using namespace std::chrono;

    using namespace rapidcsv;

    std::vector<Xy> output{};
    try {
        Document doc{filename, LabelParams{-1, -1}, SeparatorParams{'\t'}};
        auto rows = doc.GetRowCount(), cols = doc.GetColumnCount();

        steady_clock::time_point start = steady_clock::now();

        for (auto i = 0; i < rows; ++i) {
            data_t row = doc.template GetRow<double>(i);
            double y = row.back();
            row.pop_back();
            output.emplace_back(row, y);
        }

        steady_clock::time_point end = steady_clock::now();
        auto dur =  duration_cast<seconds>(end - start);

    } catch (...) {
        std::cerr << "Couldn't open file: \"" << filename << "\"!" << std::endl;
        throw;
    };
    return output;
}

std::vector<std::vector<int>> load_seeds(std::string filename)
{
    using namespace std;
    using namespace rapidcsv;
    vector<vector<int>> output{};
    try {
        Document doc{filename, LabelParams{-1, -1}, SeparatorParams{','}};
        auto rows = doc.GetRowCount(), cols = doc.GetColumnCount();

        for (auto i = 0; i < rows; ++i)
            output.push_back(doc.template GetRow<int>(i));
    } catch (...) {
        std::cerr << "Couldn't open file: \"" << filename << "\"!" << std::endl;
        throw;
    };
    return output;
}
    
}
