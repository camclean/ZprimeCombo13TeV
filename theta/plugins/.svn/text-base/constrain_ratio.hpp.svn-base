#ifndef PLUGIN_CONSTRAIN_RATIO
#define PLUGIN_CONSTRAIN_RATIO

#include "interface/phys.hpp"

/** \brief Function to use to contrain a ratio of two model parameters with a gaussian
 *
 * Configure via a setting group like
 * \code
 *   {
 *     type = "constrain_ratio";
 *     nominator = "beta_zjets";
 *     denominator = "beta_wjets";
 *     mean = 1.0;
 *     width = 0.3;
 *   };
 * \endcode
 *
 * You should make sure that the denominator can never be equal to zero.
 */
class constrain_ratio: public theta::Function{
private:
    theta::ParId pid_nominator;
    theta::ParId pid_denominator;
    double mean, width;
public:
    /// constructor used by the plugin system
    constrain_ratio(const theta::Configuration & cfg);
    /// overloaded evaluation of the function
    virtual double operator()(const theta::ParValues & v) const;
};

#endif
