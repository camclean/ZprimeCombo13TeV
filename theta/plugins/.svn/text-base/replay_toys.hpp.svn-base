#include "interface/phys.hpp"
#include "interface/plugin.hpp"
#include "interface/database.hpp"

/** \brief A data source generating toys from a theta output file written by a pseudodata_writer
 *
 * Configured via a setting group like
 * \code
 * data_source = {
 *   type = "replay_toys";
 *   name = "source";
 *   input_database = { // some input daatbase specification, e.g.:
 *      type = "sqlite_database_in";
 *      filename = "toys.db";
 *   };
 *   pdw_name = "pdw"; // optional.
 *   observables = ("o1", "o2"); // optional, default is all observables
 * };
 * \endcode
 *
 * \c type must be "replay_toys" to select this DataSource type
 *
 * \c name is required to make column names unique
 *
 * \c input_database is a specification for a DatabaseInput plugin
 *
 * \c pdw_name is the name of the pseudodata_writer instance with which the input database has been created. This is used
 *  to construct the column names in the product table. If no \c pdw_name is given, it will be guessed from the column names of the
 *  input database, which usually works but might fail in some cases.
 *
 * \c observables is a list of observable names to fill. The default is to use all observables 
 * 
 * \c rvobservables is a list of real-values observable names. The default is to use all real-valued observables.
 *
 * replay_toys will not create any columns in the products table
 */
class replay_toys: public theta::DataSource{
public:
    replay_toys(const theta::Configuration & cfg);
    virtual void fill(theta::Data & dat);    
private:
    std::auto_ptr<theta::DatabaseInput> input_database;
    std::auto_ptr<theta::DatabaseInput::ResultIterator> res;
    std::vector<theta::ObsId> observables;
    std::vector<theta::ParId> rvobservables;
};
