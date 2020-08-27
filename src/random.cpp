#include "random.hpp"
#include <iostream>

using namespace std;

int main() {
    for (int i = 0; i < 10; ++i)
        cout << rd_function(10) << " ";

    cout << "\n";

    for (int i = 0; i < 10; ++i)
        cout << rd_operator(10) << " ";
    
    cout << "\n";

    for (int i = 0; i < 10; ++i)
        cout << rd_variable(10) << " ";

    cout << "\n";

    for (int i = 0; i < 10; ++i)
        cout << rd_value() << " ";

    return 0;
}
