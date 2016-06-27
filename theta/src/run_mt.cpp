#include "interface/run_mt.hpp"
#include "interface/database.hpp"
#include "interface/histogram.hpp"
#include "interface/phys.hpp"
#include "interface/model.hpp"
#include "interface/plugin.hpp"
#include "interface/toymaker.hpp"
#include "interface/atomic.hpp"
#include "interface/redirect_stdio.hpp"

#include <iomanip>
#include <iostream>

#include <boost/thread.hpp>

using namespace theta;

run_mt::~run_mt(){}

void run_mt::worker::operator()(int n_event, atomic_int * toy_count, atomic_int * toy_error_count){
    try{
        tm->run(n_event, toy_count, toy_error_count);
    }
    catch(...){
        exception = boost::current_exception();
    }
    atomic_set(&done, 1);
}


void run_mt::run(){
    logtable->append(0, 0, LogTable::info, "run start");
    theta::atomic_int toy_count, toy_error_count;
    atomic_set(&toy_count, 0);
    atomic_set(&toy_error_count, 0);
    boost::thread_group group;
    // in case n_events < n_threads, only create n_events threads:
    int n_threads = std::min(this->n_threads, n_event);
    for(int it=0; it<n_threads; ++it){
        int nev = n_event / n_threads;
        if(it < n_event % n_threads) nev += 1;
        atomic_set(&workers[it].done, 0);
        group.create_thread(boost::bind(&worker::operator(), boost::ref(workers[it]), nev, &toy_count, &toy_error_count));
    }
    while(true){
        boost::this_thread::sleep(boost::posix_time::millisec(20));
        bool all_done = true;
        for(int it=0; it<n_threads; ++it){
            uint64_t done = atomic_get(&workers[it].done);
            if(done == 0){
                all_done = false;
                break;
            }
        }
        if(progress_listener) progress_listener->progress(atomic_get(&toy_count), n_event, atomic_get(&toy_error_count));
        if(all_done) break;
    }
    group.join_all();
    // check for exceptions:
    bool exception = false;
    std::stringstream ss;
    for(int it=0; it<n_threads; ++it){
        if(workers[it].exception){
            exception = true;
            ss << "Thread " << it << " terminated with an exception: ";
            try{
                boost::rethrow_exception(*workers[it].exception);
            }
            catch(std::exception & ex){
                ss << ex.what();
            }
            catch(...){
                ss << " (unknown exception)";
            }
        }
    }
    if(exception) throw Exception("at least one thread had an exception:\n " + ss.str());
    for(int it=0; it<n_threads; ++it){
        if(it==0){
            workers[it].bs->declare_columns(*products_table);
        }
        int runid = 1 + it;
        workers[it].bs->replay_to(*products_table, boost::bind(&ProductsTable::add_row, boost::ref(*products_table), runid, _1));
        const std::vector<ToyMaker::LogMessage> & log_messages = workers[it].tm->get_log_messages();
        for(size_t i=0; i<log_messages.size(); ++i){
            logtable->append(runid, log_messages[i].eventid, log_messages[i].severity, log_messages[i].message);
        }
    }
    
    logtable->append(0, 0, LogTable::info, "run end");
    if(log_report){
        const int* n_messages = logtable->get_n_messages();
        LogTable::e_severity s = logtable->get_loglevel();
        theta::out << std::endl << std::endl << "Log report:" << std::endl;
        theta::out << "  errors:   " << std::setw(6) << n_messages[0] << std::endl;
        if(s > 0)
            theta::out << "  warnings: " << std::setw(6) << n_messages[1] << std::endl;
        if(s > 1)
            theta::out << "  infos:    " << std::setw(6) << n_messages[2] << std::endl;
        if(s > 2)
            theta::out << "  debug:    " << std::setw(6) << n_messages[3] << std::endl;
    }
}

run_mt::run_mt(const Configuration & cfg){
    Setting s = cfg.setting;
    log_report = true;
    n_event = s["n-events"];
    if(s.exists("n_threads")){
        n_threads = s["n_threads"];
    }
    else{
        n_threads = boost::thread::hardware_concurrency();
    }
    if(n_threads <= 0){
        throw ConfigurationException("n_threads <= 0 not allowed.");
    }
    
    //setup database and tables:
    db = PluginManager<Database>::build(Configuration(cfg, s["output_database"]));

    std::auto_ptr<Table> logtable_underlying = db->create_table("log");
    logtable.reset(new LogTable(logtable_underlying));
    
    std::auto_ptr<Table> rndinfo_table_underlying = db->create_table("rndinfo");
    boost::shared_ptr<RndInfoTable> rndinfo_table(new RndInfoTable(rndinfo_table_underlying));
    cfg.pm->set("default", rndinfo_table);
    
    std::auto_ptr<Table> products_table_underlying = db->create_table("products");
    products_table.reset(new ProductsTable(products_table_underlying));
    //cfg.pm->set<ProductsSink>("default", products_table);
    
    boost::shared_ptr<int> ptr_runid(new int(1));
    cfg.pm->set("runid", ptr_runid);

    size_t n_p = s["producers"].size();
    if (n_p == 0)
        throw ConfigurationException("no producers specified!");

    //create workers
    for(int i=0; i<n_threads; ++i){
        (*cfg.pm->get<int>("runid")) = i+1;
        
        workers.push_back(worker());
        workers.back().bs.reset(new BufferingProductsSink());
        cfg.pm->set<ProductsSink>("default", workers.back().bs);
        std::auto_ptr<Model> model = PluginManager<Model>::build(Configuration(cfg, s["model"]));
        std::auto_ptr<DataSource> data_source = PluginManager<DataSource>::build(Configuration(cfg, s["data_source"]));
        boost::ptr_vector<Producer> producers;
        for (size_t i = 0; i < n_p; i++) {
             producers.push_back(PluginManager<Producer>::build(Configuration(cfg, s["producers"][i])));
        }
        workers.back().tm.reset(new ToyMaker(data_source, model, producers, workers.back().bs));
    }
    
    //logging stuff
    LogTable::e_severity level = LogTable::warning;
    if(s.exists("log-level")){
        std::string loglevel = s["log-level"];
        if(loglevel=="error") level = LogTable::error;
        else if(loglevel=="warning") level = LogTable::warning;
        else if(loglevel=="info")level = LogTable::info;
        else if(loglevel=="debug")level = LogTable::debug;
        else{
            std::stringstream ss;
            ss << "log level given in " << s["log-level"].get_path() << " unknown (given '" << loglevel << 
                  "'; only allowed values are 'error', 'warning', 'info' and 'debug')";
            throw ConfigurationException(ss.str());
        }
    }
    logtable->set_loglevel(level);
    if(s.exists("log-report")){
        log_report = s["log-report"];
    }
}

REGISTER_PLUGIN(run_mt)

