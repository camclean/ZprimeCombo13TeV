#ifndef RUN_MT_HPP
#define RUN_MT_HPP

#include "interface/decls.hpp"
#include "interface/main.hpp"
#include "interface/atomic.hpp"

#include <string>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>
#include <boost/exception_ptr.hpp>

namespace theta{

/** \brief A multithreaded version of theta::Run
 *
 * The configuration is the same as for theta::Run, the only additional
 * option is "n_threads" which specifies the number of therads to run. The default
 * is to use as many threads as the current machine has processors / cores.
 */
class run_mt: public Main{
public:
    virtual void run();
    run_mt(const Configuration & cfg);
    virtual ~run_mt();
private:

    struct worker{
        boost::shared_ptr<theta::ToyMaker> tm;
        boost::shared_ptr<theta::BufferingProductsSink> bs;
        boost::optional<boost::exception_ptr> exception;
        atomic_int done;
        void operator()(int n_event, atomic_int * toy_count, atomic_int * toy_error_count);
    };

    std::vector<worker> workers;

    boost::shared_ptr<Database> db;

    std::auto_ptr<LogTable> logtable;
    bool log_report;

    boost::shared_ptr<ProductsTable> products_table;

    //total number of events to produce, number of threads
    int n_event;
    int n_threads;
};


}

#endif

