#include <chrono>
#include <iostream>

#include "load_data.hpp"
#include "rapidcsv.h"

std::vector<entry_t> load_data(std::string filename)
{
    using namespace std;
    using namespace std::chrono;

    using namespace rapidcsv;

    std::vector<entry_t> output{};
    try {
        Document doc{filename, LabelParams{-1, -1}, SeparatorParams{'\t'}};
        auto rows = doc.GetRowCount(), cols = doc.GetColumnCount();

        for (auto i = 0; i < rows; ++i) {
            vars_t row = doc.template GetRow<double>(i);
            label_t y = row.back();
            row.pop_back();
            output.emplace_back(row, y);
        }

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
