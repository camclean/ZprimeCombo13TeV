#ifndef LBFGS_MINIMIZER_HPP
#define LBFGS_MINIMIZER_HPP

#include "interface/minimizer.hpp"
#include "liblbfgs/lbfgs.h"

//TODO:
// - new minimizer interface? Disentagle "setup" and "minimize" steps s.t.
//   anoter minimization with the same parameters and ranges will be faster ...

class lbfgs_minimizer: public theta::Minimizer{
private:
    lbfgs_parameter_t params;
public:
    lbfgs_minimizer(const theta::Configuration & cfg);
    virtual theta::MinimizationResult minimize(const theta::Function & f, const theta::ParValues & start,
                                               const theta::ParValues & step, const theta::Ranges & ranges);
    
};

#endif
