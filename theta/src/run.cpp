#include "interface/run.hpp"
#include "interface/database.hpp"
#include "interface/histogram.hpp"
#include "interface/phys.hpp"
#include "interface/model.hpp"
#include "interface/plugin.hpp"
#include "interface/redirect_stdio.hpp"

#include <iomanip>
#include <iostream>

using namespace theta;
using namespace std;

Run::~Run(){}


void Run::run(){
    //log the start of the run:
    //use eventid = 0 to indicate a "run-scoped" entry
    logtable->append(runid, 0, LogTable::info, "run start");
   
    Data data;
    int n_errors = 0;
    if(progress_listener) progress_listener->progress(0, n_event, n_errors);
    //main event loop:
    for (int eventid = 1; eventid <= n_event; eventid++) {
        if(stop_execution)break;
        try{
            data_source->fill(data);
        }
        catch(theta::Exception & ex){
           ex.message += " (in Run::run while throwing toy data)";
           throw;
        }
        logtable->append(runid, eventid, LogTable::info, "start");
        bool error = false;
        for (size_t j = 0; j < producers.size(); j++) {
            try {
                producers[j].produce(data, *model);
            } catch (Exception & ex) {
                error = true;
                std::stringstream ss;
                ss << "Producer '" << producers[j].get_name() << "' failed: " << ex.message << ".";
                logtable->append(runid, eventid, LogTable::error, ss.str());
                ++n_errors;
                break;
            }
            catch(std::logic_error & f){
                stringstream ss;
                ss << "Producer '" << producers[j].get_name() << "': " << f.what();
                throw logic_error(ss.str());
            }
        }
        //only add a row if no error ocurred to prevent NULL values and similar things ...
        if(!error){
            products_table->add_row(runid, eventid);
        }
        logtable->append(runid, eventid, LogTable::info, "end");
        if(progress_listener) progress_listener->progress(eventid, n_event, n_errors);
    }
    
    logtable->append(runid, 0, LogTable::info, "run end");
    if(log_report){
        const int* n_messages = logtable->get_n_messages();
        LogTable::e_severity s = logtable->get_loglevel();
        theta::out << endl << endl << "Log report:" << endl;
        theta::out << "  errors:   " << setw(6) << n_messages[0] << endl;
        if(s > 0)
            theta::out << "  warnings: " << setw(6) << n_messages[1] << endl;
        if(s > 1)
            theta::out << "  infos:    " << setw(6) << n_messages[2] << endl;
        if(s > 2)
            theta::out << "  debug:    " << setw(6) << n_messages[3] << endl;
    }
}

Run::Run(const Configuration & cfg){
    Setting s = cfg.setting;
    
    //0. set default values for members:
    log_report = true;
    runid = 1;
    n_event = s["n-events"];
    if(n_event <= 0){
        throw ConfigurationException("n-events <= 0 not allowed");
    }
    
    //1. setup database and tables:
    db = PluginManager<Database>::build(Configuration(cfg, s["output_database"]));

    std::auto_ptr<Table> logtable_underlying = db->create_table("log");
    logtable.reset(new LogTable(logtable_underlying));
    
    std::auto_ptr<Table> rndinfo_table_underlying = db->create_table("rndinfo");
    boost::shared_ptr<RndInfoTable> rndinfo_table(new RndInfoTable(rndinfo_table_underlying));
    cfg.pm->set("default", rndinfo_table);
    
    std::auto_ptr<Table> products_table_underlying = db->create_table("products");
    products_table.reset(new ProductsTable(products_table_underlying));
    cfg.pm->set<ProductsSink>("default", products_table);
    
    boost::shared_ptr<int> ptr_runid(new int(runid));
    cfg.pm->set("runid", ptr_runid);
        
    //2. model and data_source
    model = PluginManager<Model>::build(Configuration(cfg, s["model"]));
    data_source = PluginManager<DataSource>::build(Configuration(cfg, s["data_source"]));
    
    //3. logging stuff
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
    
    //4. producers:
    size_t n_p = s["producers"].size();
    if (n_p == 0)
        throw ConfigurationException("no producers specified!");
    for (size_t i = 0; i < n_p; i++) {
         producers.push_back(PluginManager<Producer>::build(Configuration(cfg, s["producers"][i])));
    }
}

REGISTER_PLUGIN_DEFAULT(Run)


