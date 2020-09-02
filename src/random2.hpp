#include <random>
#include <string>
#include <vector>

std::mt19937 _gen; // <-- It's a singleton $%$#$%@%!!!

class random_t
{
public:
    random_t(std::vector<int> seeds);

    int variable(int from, int to) const;
    int binary() const;
    int function(int from, int to) const;
    bool choose() const;
        
    double real(double m1, double m2, double stddev = 1.0) const;
    int integer(int from, int to) const;
    operator std::string() const;

    template <typename T>
    std::vector<T> sample(const std::vector<T>& population,
                          std::size_t n_samples);
private:
    std::vector<int> _seeds_v;
    std::seed_seq _seeds;
};
