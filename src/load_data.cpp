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
    
    Document doc{filename, LabelParams{-1, -1}, SeparatorParams{'\t'}};
    auto rows = doc.GetRowCount(), cols = doc.GetColumnCount();
    cout << rows << "×" << cols << endl;
    cout << "starting read." << endl;

    steady_clock::time_point start = steady_clock::now();

    std::vector<Xy> output{};
    for (auto i = 0; i < rows; ++i) {
        data_t row = doc.template GetRow<double>(i);
        double y = row.back();
        row.pop_back();
        output.emplace_back(row, y);
    }

    steady_clock::time_point end = steady_clock::now();
    auto dur =  duration_cast<seconds>(end - start);
    cout << "finished [" << dur.count() << "s]" << endl;

    return output;
}

/*
int main(int argc, char **argv) {
    if (argc != 2) {
        std::cout << "Just one file name" << std::endl;
        return 1;
    }
    
    std::string filename = argv[1];
    std::vector<Xy> data = load_data(filename);

    for (auto row: data) {
        for (auto x: row.first)
            std::cout << x << " ";
        std::cout << " -> " << row.second << std::endl;
    }

    return 0;
}
*/
}
