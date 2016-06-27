#ifndef RUN_HPP
#define RUN_HPP

#include "interface/decls.hpp"
#include "interface/main.hpp"

#include <string>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>

namespace theta{


/** \brief A set of producers executed on (pseudo-)data pulled from a DataSource
 *
 * This is the central class in theta which "wires together" all required components and
 * implements the main algorithm of theta which consists of
 * <ol>
 *   <li>Get (pseudo-) data from a DataSource</li>
 *   <li>Pass this data to a list of producers which in turn write their results to the 'products' table</li>
 *   <li>Repeat steps 1 and 2 as many time as configured (\c n-events)</li>
 * </ol>
 *
 * Furthermore, this class serves as container for many objects shared across
 * one run, namly database and table objects, the random number generator, and the
 * (common) Model used by all producers.
 *
 * The configuration is done via a setting like:
 * \code
 * main = {
 *   data_source = {...}; //some data source definition
 *   output_database = {...}; // some database definition
 *
 *   producers = ("@hypotest", "@hypotest2");
 *   n-events = 10000;
 *   model = "@gaussoverflat";
 *
 *   //optional:
 *   log-level = "error";  //default is "warning"
 *   log-report = false;  //default is true
 * };
 *
 * hypotest = {...}; //some producer definition
 * hypotest2 = {...}; //some other producer definition
 * gaussoverflat = {...}; // some model definition
 * \endcode
 *
 * \c data_source is a definition of a \link theta::DataSource DataSource \endlink implementation to use
 *    as source of the data / pseudodata to run the producers on
 *
 * \c output_database is a definition of a \link theta::Database Database \endlink implementation.
 *
 * \c producers is a list of links within the configuration file to setting groups which contain the
 *    definition of the producers to run.
 *
 * \c n-events is the number of pseudo experiments to run. This is ignored for some subclasses.
 *
 * \c model is the model used for the producers: the producers will be given this for pseudo data creation.
 *
 * \c log-level controls the amount of logging information written to the log table: only log messages with a
 *      severity level equal to or exceeding the level given here are actually logged. Valid values
 *      are "error", "warning", "info" and "debug". Note that it is not possible to disable logging of messages with severity
 *      "error".
 *
 * \c log-report is a boolean specifying whether or not to print a log report to standard output at
 *       the end of the run. This report summarizes how many messages there have been from any non-suppressed
 *       level. This allows for a quick check by the user whether everything went Ok or whether there
 *       have been obvious errors.
 *
 *  Handling of result tables is done in the individual producers. Only run-wide tables
 *  are managed here, that is
 *  <ul>
 *   <li>A \link LogTable LogTable \endlink called 'log', where all log entries are stored.</li>
 *   <li>A \link RndInfoTable RndInfoTable \endlink called 'rndinfo', where information about used random seeds by all plugin
 *     which consume randomness are saved. This can be useful if one wants to reproduce the exact same run (e.g., the same
 *     pseudo data between different theta runs or to provoke a certain rarely ocurring bug). </li>
 *   <li>A \link ProductsTable ProductsTable \endlink called 'events' where the per-event results from all
 *     the producers are stored and other per-event information.</li>
 *  </ul>
 *
 *  For more information about these tables, refer to the documentation of the corresponding Table classes.
 */
class Run: public Main{
public:

    /** \brief Perform the actual run.
     * 
     * In a pseudo-experiment loop, ask the configured data_source for data and run
     * all the producers on it, using the configured model.
     */
    virtual void run();
    
    /** \brief Construct a Run from the given configuration
     */
    Run(const Configuration & cfg);
    
    virtual ~Run();
    
private:

    std::auto_ptr<Model> model;
    //note: the tables hold a shared_ptr to this database to ensure proper destruction order
    boost::shared_ptr<Database> db;
    std::auto_ptr<DataSource> data_source;

    std::auto_ptr<LogTable> logtable;
    bool log_report;

    //the producers to be run on the pseudo data:
    boost::ptr_vector<Producer> producers;
    boost::shared_ptr<ProductsTable> products_table;

    //the runid, and the total number of events to produce:
    int runid;
    int n_event;
};


}

#endif

