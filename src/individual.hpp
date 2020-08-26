#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <functional>

// index from 0 to 15.
std::vector<std::string> _operators{"+", "-", "*", "/", "^"};

// index from 16 to 127.
std::vector<std::string> _functions{"sin", "cos", "tan", "sinh", "cosh", "tanh",
    "log", "ln"};

// index from 128 to 255.
std::vector<std::string> _variables{"x_0", "x_1", "x_2"};

/* DEFINITIONS */
constexpr char NULL_VALUE = 0b1111'1111;
constexpr double BIG_VALUE = 1000.0;
using individual_t = std::vector<char>;

/* HELPER FUNCTIONS */
bool is_operator(char cell) { return cell < 16; }
bool is_function(char cell) { return cell < 127 && cell > 15; }
bool is_variable(char cell) { return cell > 127; }

// Initialized all elements of individual to null value.
void init_individual(individual_t& individual)
{
    for (auto& element: individual)
        element = NULL_VALUE;
}

void substitute_variables_values(
                                 individual_t& individual,
                                 std::vector<double> variables_values)
{
    for (auto& element: individual) {
        if (is_variable(element)) {
            auto index = element - 127;
            element = variables_values[index];
        }
    }
}

double eval_function(std::string function_name, double x)
{
    if (function_name == "sin") {
        return sin(x);
    } else if (function_name == "cos") {
        return cos(x);
    } else if (function_name == "tan") {
        if (cos(x) == 0) return BIG_VALUE;
        return tan(x);
    } else if (function_name == "sinh") {
        return sinh(x);
    } else if (function_name == "cosh") {
        return cosh(x);
    } else if (function_name == "tanh") {
        return tanh(x);
    } else if (function_name == "log") {
        return log10(max(0.1, x));
    } else if (function_name == "ln") {
        return log(max(0.1, x));
    }
}

double eval_operator(std::string operator_symbol, double lop, double rop)
{
    if (operator_symbol == "+") {
        return lop + rop;
    } else if (operator_symbol == "-") {
        return lop - rop;
    } else if (operator_symbol == "*") {
        return lop * rop;
    } else if (operator_symbol == "/") {
        if (rop == 0) return BIG_VALUE;
        return lop / rop;
    } else if (operator_symbol == "^") {
        return pow(lop, rop);
    } 
}

double recursive_eval(individual_t& individual, std::size_t position, std::vector<double> variables)
{
    auto& element = individual[position];
    if (is_function(element)) {
        auto function_name = _functions[element - 16];
        double x = recursive_eval(individual, pow(2, position));
        return eval_function(op);
    } else if (is_operator(element)) {
        auto operator_symbol = _operators[element - 0];
        double lop = recursive_eval(individual, pow(2, position));
        double rop = recursive_eval(individual, pow(2, position) + 1);
        return eval_operator(lop, rop);
    } else if (is_variable(element)) {
        return element;
    } else { return 0; }
}

double eval(individual_t individual, std::vector<double> variables_values)
{
    substitute_variables_values(individual, variables_values);
}
