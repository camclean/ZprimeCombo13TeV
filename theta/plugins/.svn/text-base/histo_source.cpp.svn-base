#include "plugins/histo_source.hpp"
#include "interface/histogram-function.hpp"

using namespace theta;
using namespace std;

histo_source::histo_source(const Configuration & cfg): DataSource(cfg){
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    ObsIds oids = vm->get_all_observables();
    for(ObsIds::const_iterator oit=oids.begin(); oit!=oids.end(); ++oit){
        string obs_name = vm->get_name(*oit);
        if(not cfg.setting.exists(obs_name)) continue;
        std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(Configuration(cfg, cfg.setting[obs_name]));
        if(hf->get_parameters().size() > 0){
            throw ConfigurationException("histo_source: given histogram depends on parameters, which is not allowed");
        }
        size_t nbins;
        double xmin, xmax;
        hf->get_histogram_dimensions(nbins, xmin, xmax);
        data[*oit].reset(nbins, xmin, xmax);
        hf->add_with_coeff_to(data[*oit], 1.0, ParValues());
    }
    if(cfg.setting.exists("rvobs-values")){
        ParValues rvobs_values;
        for(size_t i=0; i<cfg.setting["rvobs-values"].size(); ++i){
            string name = cfg.setting["rvobs-values"][i].get_name();
            ParId id = vm->get_par_id(name);
            // type checking:
            if(vm->get_type(id)!="rvobs"){
                throw ConfigurationException("type error: parameter '" + name + "' used as real-valued observable, but was not declared as such.");
            }
            rvobs_values.set(id, cfg.setting["rvobs-values"][i]);
        }
        data.set_rvobs_values(rvobs_values);
    }
}

void histo_source::fill(Data & dat){
    dat = data;
}

REGISTER_PLUGIN(histo_source)
