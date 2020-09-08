#include <string>
#include <vector>
#include <map>
#include <tuple>

#include "random.hpp"

using vars_t = std::vector<double>;

enum class class_t: unsigned char { oper, func, var, cons, null }; // operator, function, var, const.

struct gene_t
{
    class_t _code;
    unsigned char _value;
    bool operator<(const gene_t& other) const;
};

namespace std
{
    template <>
    struct hash<gene_t>
    {
        std::size_t operator()(const gene_t& gene)
        {
            std::size_t code = static_cast<std::size_t>(gene._code);
            std::size_t value = static_cast<std::size_t>(gene._value);
            return (code << 8) + value;
        };
    };
}

class individual_t
{
    friend class individual_handler_t;
public:
    individual_t(std::size_t depth);
    individual_t(const individual_t& other);
    ~individual_t() = default;
    individual_t& operator=(const individual_t& other);

    void clear();
    void clear_subtree(std::size_t node);
    void clear_lsubtree(std::size_t node);
    void clear_rsubtree(std::size_t node);
        
    std::size_t depth() const;        
    std::vector<gene_t> _gen;
private:

    bool has_child(std::size_t node) const;
    bool has_parent(std::size_t node) const;
    std::size_t lchild(std::size_t parent) const;
    std::size_t rchild(std::size_t parent) const;
    std::size_t parent(std::size_t child) const;
};

class individual_handler_t
{
    using pFunc = double(*)(double);
    using pOper = double(*)(double, double);

public:
    individual_handler_t(const random_t& rd);
    ~individual_handler_t() = default;
    
    void add_operator(std::string repr, std::size_t precedence, pOper oper);
    void add_function(std::string repr, pFunc func);    
    void add_variable(std::size_t n_vars);

    gene_t rd_func() const;
    gene_t rd_oper() const;
    gene_t rd_var() const;
    gene_t rd_cons() const;
    
    double eval(const individual_t& ind, vars_t vars) const;
    std::string str(const individual_t& ind) const;
    
private:
    double eval(const individual_t& ind, vars_t vars, std::size_t n) const;
    double eval_oper(const gene_t& g, double x, double y) const;
    double eval_func(const gene_t& g, double x) const;
    double eval_cons(const gene_t& g) const;

    std::string str(const individual_t& ind, std::size_t n) const;
    bool need_parentheses(const individual_t& ind, std::size_t n) const;
    std::string str_oper(const gene_t& g) const;
    std::string str_func(const gene_t& g) const;
    std::string str_cons(const gene_t& g) const;
    std::string str_var(const gene_t& g) const;

    std::map<gene_t, std::tuple<std::string, pFunc>> _func;
    std::map<gene_t, std::tuple<std::string, std::size_t, pOper>> _oper;
    std::size_t _n_vars;

    const random_t& _rd;
};
