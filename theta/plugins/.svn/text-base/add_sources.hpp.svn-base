#include "interface/phys.hpp"
#include "interface/plugin.hpp"

#include <boost/ptr_container/ptr_vector.hpp>

/** \brief A data source adding the result of some other data sources
 *
 * Configured via a setting group like
 * \code
 * data_source = {
 *   type = "add_sources";
 *   sources = ("@source1", "@source2");
 * };
 *
 * source1 = { / * some DataSource definition* / };
 * source2 = { / * some DataSource definition* / };
 * \endcode
 *
 * \c type must be "add_sources" to select this DataSource type
 *
 * \c sources is a list of data source specifications.
 *
 * The data produced by this DataSource is -- for each observable -- the sum of all sources. The sources
 * may provide data for different observables; missing observables in a source are treated as if it
 * would return a histogram with zero entries only.
 *
 * This DataSource was created in an anlysis in which we were asked to check what happens
 * if we deliberately add very few signal-like events to the background model. The background of this question
 * was whether these (few) events where enough to explain some discrepancy we saw between the expected
 * and observed statistical result. This check can be done by using this plugin and specifying as 'sources',
 * (i) a model_source which produces background only data,
 * (ii) a histo_source which contains the (fixed) signal-like data to add.<br/>
 * I do not no whether it has other uses as well.
 *
 * add_sources will not create any columns in the products table.
 */
class add_sources: public theta::DataSource{
public:

    /// Construct from a Configuration; required by the plugin system
    add_sources(const theta::Configuration & cfg);

    /** \brief Fills the provided Data instance with data from the model
     */
    virtual void fill(theta::Data & dat);
    
private:
    boost::ptr_vector<theta::DataSource> sources;
};
