#ifndef PLUGINS_VARY_ONE_HPP
#define PLUGINS_VARY_ONE_HPP

#include "interface/distribution.hpp"
#include "interface/exception.hpp"
#include <vector>
#include <map>


/** \brief distribution varying one parameter at a time
 *
 *
 * Useful for scanning one parameter at a time, fixing others. Only for sampling, do not
 * use as part of a likelihood!
 *
 *
 * \code
 * dist = {
 *   type = "vary_one";
 *   p1 = {  // assuming p1 was declared as parameter
 *      default = 0.0;
 *      values = (-1.0, 1.0);
 *   };
 *   p2 = { // assuming p2 was declared as parameter
 *      default = 7.0;
 *      values = (3.0);
 *   };
 * };
 * \endcode
 *
 * \c values can be omitted; in this case, the behaviour is the same as specifying the empty list, i.e., the value for
 *  this parameter will always be the default.
 *
 * Sampling from this distribution is not random but will, in this order:
 * <ol>
 *  <li>return a set of values set all values to the defaults</li>
 *  <li>defaults for all parameters except the first, which is set to the first alternate value. Subsequently,
 *      all alternate values for this parameter are returned.</li>
 *  <li>defaults for all parameters except the <em>second</em>, which is set to the first alternate value. Subsequently,
 *      all alternate values for this parameter are returned.</li>
 *  <li>etc.</li>
 * </ol>
 *
 * For the example above, the values returned for (p1, p2) will be: (0, 7), (-1, 7), (1, 7), (0, 3), and then from the beginning again:
 * (0, 7), etc. The returned values are repeated after the sum of the length of values plus one (all defaults).
 */
class vary_one: public theta::Distribution{
public:
    virtual void sample(theta::ParValues & result, theta::Random & rnd) const;
    virtual void mode(theta::ParValues & result) const{
      throw std::invalid_argument("vary_one::mode is not implemented");
    }
    virtual double eval_nl(const theta::ParValues & values) const{
       throw std::invalid_argument("vary_one::eval_nl is not implemented");
    }
    virtual const std::pair<double, double> & support(const theta::ParId & p) const{
       throw std::invalid_argument("vary_one::support is not implemented");
    }
    virtual double width(const theta::ParId & p) const{
       throw std::invalid_argument("vary_one::width is not implemented");
    }
    vary_one(const theta::Configuration & cfg);
    
private:

    theta::ParValues default_values;
    std::vector<std::pair<theta::ParId, std::vector<double> > > other_values;
    mutable size_t next_index;
    //the total number of possibilities, i.e. other values lengths sum + 1
    size_t n_total;
};

#endif
