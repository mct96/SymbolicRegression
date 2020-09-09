#include <string>
#include <vector>
#include <map>
#include <tuple>

#include "random.hpp"
#include "common_types.hpp"

enum class class_t: unsigned char { oper, func, var, cons, null }; // operator, function, var, const.

struct gene_t
{
    class_t _code;
    unsigned char _value;
    bool operator<(const gene_t& other) const;
    bool operator==(const gene_t& other) const;
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

    std::size_t size() const;
    
    const gene_t& operator[](std::size_t i) const;
    gene_t& operator[](std::size_t i);
    
    void clear();
    void clear_subtree(std::size_t node);
    void clear_lsubtree(std::size_t node);
    void clear_rsubtree(std::size_t node);

    void copy_subtree(std::size_t to,
                      const individual_t& other,
                      std::size_t from);

    static std::size_t depth(std::size_t node);
    
    std::size_t max_depth() const;
    std::size_t depth() const;        
    std::vector<gene_t> _gen;
private:
    std::size_t _depth;
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

    void random_generator(const random_t& rd);
    const random_t& random_generator() const;
    
    void add_operator(std::string repr, std::size_t precedence, pOper oper);
    void add_function(std::string repr, pFunc func);    
    void add_variable(std::size_t n_vars);

    gene_t rd_func() const;
    gene_t rd_oper() const;
    gene_t rd_var() const;
    gene_t rd_cons() const;
    gene_t rd_term() const;
    std::size_t rd_point(const individual_t& ind) const;
    
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
