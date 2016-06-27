#ifndef PLUGIN_EQUIDISTANT_DELTAS
#define PLUGIN_EQUIDISTANT_DELTAS

#include "interface/plugin.hpp"
#include "interface/distribution.hpp"



/** \brief A sum of equidistant delta functions in one dimension
 *
 * It is configured with a setting group like
 * \code
 * {
 *  type = "equidistant_deltas";
 *  parameter = "p0";
 *  range = [0.0, 10.0];
 *  n = 101;
 * };
 * \endcode
 *
 * \c parameter specifies the parameter the normal distribution depends on
 *
 * \c range and \c n specify the range of the delta functions. The first and last value of the range will
 *  be both included as values, therefore, \c n has to be &gt;=2.
 *
 * 
 *
 */
class equidistant_deltas: public theta::Distribution{
public:
    /// \brief Constructor used by the plugin system to build an instance from settings in a configuration file
    equidistant_deltas(const theta::Configuration & cfg);
    
    virtual void sample(theta::ParValues & result, theta::Random & rnd) const;
    virtual void mode(theta::ParValues & result) const;
    virtual double eval_nl(const theta::ParValues & values) const;
    virtual const std::pair<double, double> & support(const theta::ParId&) const;
    
private:
    unsigned int n;
    double low, width_;
    std::pair<double, double> support_;
};


#endif
