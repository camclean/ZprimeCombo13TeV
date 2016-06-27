#ifndef THETA_ROOT_MINUIT1_HPP
#define THETA_ROOT_MINUIT1_HPP

#include "interface/minimizer.hpp"

/** \brief Minimizer using the MINUIT minimizer from root
 * 
 * This is similar to \c root_minimizer, which uses the MINUIT2 interface; \c root_minuit1 uses the TMinuit
 * interface.
 *
 * Configuration with a setting like:
 * \code
 * {
 *  type = "root_minuit1";
 * 
 *  tolerance = 0.1; //optional
 *  strategy = 1; // optional
 *  n_retries = 10; // optional; default is 2
 *  infinity = 1e5; // optional; this is the default
 *  hesse = true; // optional, default is false
 * }
 * \endcode
 *
 * \c tolerance is the tolerance.
 * 
 * \c strategy is the minuit strategy ('SET STR'). From the Minuit documentation:
 *   0: MINIMIZE THE NUMBER OF CALLS TO FUNCTION
 *   1: TRY TO BALANCE SPEED AGAINST RELIABILITY
 *   2: MAKE SURE MINIMUM TRUE, ERRORS CORRECT
 *
 * \c n_retries is the maximum number of retries in case the first minimization attempt fails. In some cases,
 *    the success rate for the minimization can be increased by re-running the minimization starting at the parameter values
 *    of the current (failed) attempt.
 * 
 * \c infinity unfortunately, Minuit does not seem to support parameters with half-infinite bounds (i.e.: one interval end is finite, the other +-infinity).
 *  The value given for \c infinity here is multiplied with step size and used as the unbounded interval end.
 * 
 * \c hesse if true, run minuit HESse calculation at the minimum for a more accurate covariance matrix
 *
 */
class root_minuit1: public theta::Minimizer{
public:
    /// Constructor used by the plugin system
    root_minuit1(const theta::Configuration & cfg);
    virtual theta::MinimizationResult minimize(const theta::Function & f, const theta::ParValues & start,
            const theta::ParValues & step, const theta::Ranges & ranges);
private:
    
    double tolerance, infinity;
    int strategy;
    unsigned int n_retries;
    bool hesse;
};

#endif

