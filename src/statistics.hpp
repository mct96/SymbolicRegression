#include <vector>
#include <string>

class statistics_t
{
public:
    statistics_t(std::vector<double> values);

    void update_values();
    
    double mean() const;
    double median() const;
    double max() const;
    double min() const;
    double var() const;
    double stddev() const;
    double total() const;

    operator std::string() const;
private:
    std::vector<double> _values;
    double _mean;
    double _median;
    double _max;
    double _min;
    double _var;
    double _stddev;
    double _total;
};
          
