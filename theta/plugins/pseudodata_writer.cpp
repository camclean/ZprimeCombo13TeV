#include "plugins/pseudodata_writer.hpp"
#include "interface/plugin.hpp"
#include "interface/phys.hpp"
#include "interface/data.hpp"
#include "interface/histogram.hpp"

using namespace theta;
using namespace std;

void pseudodata_writer::produce(const Data & data, const Model & model) {
    for(size_t i=0; i<observables.size(); ++i){
        const Histogram1D & h = data[observables[i]];
        double n_event = h.get_sum();
        products_sink->set_product(n_events_columns[i], n_event);
        if(write_data){
            products_sink->set_product(data_columns[i], h);
        }
    }
    const ParValues & rvobs_values = data.get_rvobs_values();
    for(size_t i=0; i<rvobservables.size(); ++i){
        double val = rvobs_values.get(rvobservables[i], numeric_limits<double>::quiet_NaN());
        products_sink->set_product(rvobs_columns[i], val);
    }
}

pseudodata_writer::pseudodata_writer(const theta::Configuration & cfg): Producer(cfg){
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
    write_data = cfg.setting["write-data"];
    for(size_t i=0; i<observables.size(); ++i){
        n_events_columns.push_back(products_sink->declare_product(*this, "n_events_" + vm->get_name(observables[i]), theta::typeDouble));
        if(write_data){
            data_columns.push_back(products_sink->declare_product(*this, "data_" + vm->get_name(observables[i]), theta::typeHisto));
        }
    }
    ParIds pars = vm->get_all_parameters();
    for(ParIds::const_iterator it=pars.begin(); it!=pars.end(); ++it){
        if(vm->get_type(*it) != "rvobs") continue;
        rvobservables.push_back(*it);
        rvobs_columns.push_back(products_sink->declare_product(*this, vm->get_name(*it), theta::typeDouble));
    }
}

REGISTER_PLUGIN(pseudodata_writer)
