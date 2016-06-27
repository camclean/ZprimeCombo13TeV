#include "plugins/add_sources.hpp"
#include "interface/data.hpp"

using namespace theta;

add_sources::add_sources(const theta::Configuration & cfg): DataSource(cfg){
    Setting s = cfg.setting["sources"];
    size_t n = s.size();
    if(n==0) ConfigurationException("have to specify at least one source in 'sources'");
    for(size_t i=0; i<n; ++i){
        sources.push_back(PluginManager<DataSource>::build(Configuration(cfg, s[i])));
    }
}

void add_sources::fill(theta::Data & dat){
    dat.reset();
    sources[0].fill(dat);
    ParValues rvobs_values;
    rvobs_values.set(dat.get_rvobs_values());
    ObsIds obs_result = dat.get_observables();
    ObsIds obs_all = obs_result;
    for(size_t i=1; i<sources.size(); ++i){
        theta::Data data_tmp;
        sources[i].fill(data_tmp);
        rvobs_values.set(data_tmp.get_rvobs_values());
        ObsIds obs_tmp = data_tmp.get_observables();
        obs_all.insert_all(obs_tmp);
        for(ObsIds::const_iterator it=obs_all.begin(); it!=obs_all.end(); ++it){
            // add data_tmp[*it]  to  dat[*it]. However, some might be empty, so check:
            if(!obs_tmp.contains(*it)) continue;
            if(obs_result.contains(*it)){
               dat[*it] += data_tmp[*it];
            }
            else{
               dat[*it] = data_tmp[*it];
            }
        }
    }
    dat.set_rvobs_values(rvobs_values);
}

REGISTER_PLUGIN(add_sources)
