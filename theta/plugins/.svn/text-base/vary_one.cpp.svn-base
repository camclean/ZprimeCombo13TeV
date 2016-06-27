#include "plugins/vary_one.hpp"
#include "interface/plugin.hpp"

using namespace std;
using namespace theta;


void vary_one::sample(ParValues & result, Random &) const{
    result.set(default_values);
    if(next_index==0){
       next_index = 1 % n_total;
       return;
    }
    size_t k = next_index - 1;
    for(size_t i=0; i<other_values.size(); ++i){
       if(k < other_values[i].second.size()){
          result.set(other_values[i].first, other_values[i].second[k]);
          next_index = (next_index + 1) % n_total;
          return;
       }
       else{
          k-= other_values[i].second.size();
       }
    }
    theta_assert(false);
}

vary_one::vary_one(const theta::Configuration & cfg): next_index(0), n_total(1){
   size_t n = cfg.setting.size();
   boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
   for(size_t i=0; i<n; ++i){
       string parname = cfg.setting[i].get_name();
       if(parname=="type") continue;
       ParId pid = vm->get_par_id(parname);
       par_ids.insert(pid);
       default_values.set(pid, cfg.setting[i]["default"]);
       std::vector<double> values;
       if(cfg.setting[i].exists("values")){
           size_t N = cfg.setting[i]["values"].size();
           for(size_t k=0; k<N; ++k){
               values.push_back(cfg.setting[i]["values"][k]);
           }
           n_total += N;
       }
       other_values.push_back(make_pair(pid, values));
   }
}


REGISTER_PLUGIN(vary_one)
