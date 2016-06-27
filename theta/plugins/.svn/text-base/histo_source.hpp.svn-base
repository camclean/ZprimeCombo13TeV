#ifndef PLUGIN_CORE_HPP
#define PLUGIN_CORE_HPP

#include "interface/plugin.hpp"
#include "interface/phys.hpp"
#include "interface/data.hpp"

/** \brief A data source using a list of constant Histograms
 *
 * It will always yield the same data: which is given as constant histogram functions (i.e.,
 * HistogramFunction instances which do not depend on any parameter).
 *
 * Configured via a setting group like
 * \code
 * data_source = {
 *   type = "histo_source";
 *   name = "source";
 *   obs1 = { // assuming "obs1" is an observable
 *        type = "root_histogram"; filename = "DATA.root"; histoname = "h_obs1";
 *    }; // or some other constant histogram specification
 *   obs2 = { ... }; // some constant histogram specification to use for observable obs2
 *
 *   rvobs-values = { rv1 = 2.0; rv2 = 3.0; }; //optional
 * };
 * \endcode
 *
 * \c type must be "histo_source" to select this DataSource type
 *
 * \c name is the name of the data source. (It is not actually used for this type but every DataSource
 * requires the name setting).
 *
 * \c rvobs-values is optional and specifies the values for the real-valued observables.
 *
 * For all observables, a setting specifying the (constant) HistogramFunction must be specified,
 * with the name of that observable.
 *
 * This DataSource does not save any values to the EventTable.
 */
class histo_source: public theta::DataSource{
public:

    /// Construct from a Configuration; required by the plugin system
    histo_source(const theta::Configuration & cfg);

    virtual void fill(theta::Data & dat);
    
private:
    theta::Data data;
};

#endif

