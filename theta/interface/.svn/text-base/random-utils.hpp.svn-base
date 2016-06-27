#ifndef RANDOM_UTILS_HPP
#define RANDOM_UTILS_HPP

#include "interface/decls.hpp"
#include <string>
#include <memory>

namespace theta{

/// \brief Base class for plugins using a random number generator.
class RandomConsumer{
public:
   virtual ~RandomConsumer();
protected:
   /** \brief Constructor to be used by derived classes
    *
    * Will save the random seed in the RndInfoTable of the cfg.pm, if this is set.
    */
   RandomConsumer(const Configuration & cfg, const std::string & name);
   
   /// random seed used
   int seed;
   
   /// random number generator instance to be used by derived classes
   std::auto_ptr<Random> rnd_gen;
};


/** \brief Replace data in the supplied DoubleVector by a Poisson value
 *
 * Replaces each data value by a Poisson with mean equal to the original value.
 */
void randomize_poisson(DoubleVector & d, Random & rnd);

/** \brief Smear the value in each bin by a Gaussian within its uncertainty
 * 
 * The returned Histogram1D consists of Gaussian random values drawn according to the values and uncertainties in histo.
 * The values are drawn according to truncated Gaussians, i.e., values are drawn until the value is >=0, unless the value itself
 * is < 0.
 */
Histogram1D randomize_gauss(const Histogram1DWithUncertainties & histo, Random & rnd);


}


#endif

