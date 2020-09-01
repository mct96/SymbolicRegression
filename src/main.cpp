#include "symbolic_regression.hpp"

#include <iostream>
#include <fstream>
#include <string>

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
    if (argc <= 1) {
        usage_message();
        return 1;
    }

    
    return 0;
}
