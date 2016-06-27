#ifndef PLUGIN_NL_ONE_OVER_SQRT
#define PLUGIN_NL_ONE_OVER_SQRT

#include "interface/phys.hpp"

/** \brief Function returning -log(1/sqrt(p)) = 0.5 * p where p is a configurable parameter
 *
 * This is the example class discussed on the page \subpage extend "Extending theta".
 */
class nl_one_over_sqrt: public theta::Function{
private:
    theta::ParId pid;

public:
    /// constructor for the plugin system
    nl_one_over_sqrt(const theta::Configuration & cfg);
    /// overloaded evaluation operator of theta::Function
    virtual double operator()(const theta::ParValues & v) const;
};


#endif
