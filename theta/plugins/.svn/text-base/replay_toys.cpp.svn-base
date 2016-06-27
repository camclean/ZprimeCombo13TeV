#include "plugins/replay_toys.hpp"
#include "interface/data.hpp"
#include "interface/histogram.hpp"

using namespace theta;
using namespace std;

namespace{
    
std::string guess_pdw_name(DatabaseInput & in){
    std::vector<std::pair<std::string, data_type> > cols = in.get_all_columns("products");
    //look which column names have "__data_" as part of their name and guess the initial part.
    for(size_t i=0; i<cols.size(); ++i){
        size_t p = cols[i].first.find("__data_");
        if(p!=string::npos){
            return cols[i].first.substr(0, p);
        }
    }
    throw ConfigurationException("could not guess pdw_name; please spcify it explicitly.");
}
    
    
}

void replay_toys::fill(theta::Data & dat){
    if(!res->has_data()){
        throw logic_error("replay_toys: no more input data available");
    }
    for(size_t i=0; i<observables.size(); ++i){
        dat[observables[i]] = res->get_histogram(i);
    }
    const size_t n_obs = observables.size();
    ParValues rvobs;
    for(size_t i=0; i<rvobservables.size(); ++i){
        rvobs.set(rvobservables[i], res->get_double(n_obs + i));
    }
    dat.set_rvobs_values(rvobs);
    ++(*res);
}

replay_toys::replay_toys(const theta::Configuration & cfg): DataSource(cfg){
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    if(cfg.setting.exists("observables")){
        size_t n = cfg.setting["observables"].size();
        observables.reserve(n);
        for(size_t i=0; i<n; i++){
            observables.push_back(vm->get_obs_id(cfg.setting["observables"][i]));
        }
    }
    else{
        ObsIds oids = vm->get_all_observables();
        observables.reserve(oids.size());
        for(ObsIds::const_iterator oit=oids.begin(); oit!=oids.end(); ++oit){
            observables.push_back(*oit);
        }
    }
    if(cfg.setting.exists("rvobservables")){
        size_t n = cfg.setting["rvobservables"].size();
        rvobservables.reserve(n);
        for(size_t i=0; i<n; i++){
            std::string rvobs_name = cfg.setting["rvobservables"][i];
            rvobservables.push_back(vm->get_par_id(rvobs_name));
            if(vm->get_type(rvobservables.back()) != "rvobs"){
                throw ConfigurationException("invalid type for '" + rvobs_name + "': expected 'rvobs'");
            }
        }
    }
    else{
        ParIds pids = vm->get_all_parameters();
        for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
            if(vm->get_type(*it) == "rvobs"){
                rvobservables.push_back(*it);
            }
        }
    }
    
    input_database = PluginManager<DatabaseInput>::build(Configuration(cfg, cfg.setting["input_database"]));
    
    std::vector<std::string> colnames; // column names in the same order as "observables"
    string pdw_name;
    if(cfg.setting.exists("pdw_name")){
        pdw_name = static_cast<string>(cfg.setting["pdw_name"]);
    }
    else{
        pdw_name = guess_pdw_name(*input_database);
    }
    for(size_t i=0; i<observables.size(); ++i){
        colnames.push_back(pdw_name + "__data_" + vm->get_name(observables[i]));
    }
    for(size_t i=0; i<rvobservables.size(); ++i){
        colnames.push_back(pdw_name + "__" + vm->get_name(rvobservables[i]));
    }
    res = input_database->query("products", colnames);
    if(!res->has_data()){
        throw ConfigurationException("no input data available");
    }
}

REGISTER_PLUGIN(replay_toys)
