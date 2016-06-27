#include "plugins/sys_rate_function.hpp"
#include "interface/plugin.hpp"

using namespace std;
using namespace theta;

sys_rate_function::sys_rate_function(const theta::Configuration & cfg){
    unsigned int npar_f = cfg.setting["factors"].size();
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    for(unsigned int i=0; i<npar_f; ++i){
       theta::ParId pid = vm->get_par_id(cfg.setting["factors"][i]);
       f_pids.push_back(pid);
       par_ids.insert(pid);
    }
    
    unsigned int npar_s = cfg.setting["sys_rates"].size();
    for(unsigned int i=0; i<npar_s; ++i){
       theta::ParId pid = vm->get_par_id(cfg.setting["sys_rates"][i][0]);
       s_pids.push_back(pid);
       par_ids.insert(pid);
       r_minus.push_back(cfg.setting["sys_rates"][i][1]);
       r_plus.push_back(cfg.setting["sys_rates"][i][2]);
    }
}

double sys_rate_function::operator()(const theta::ParValues & values) const{
    double result = 1.0;
    for(size_t i=0; i<s_pids.size(); ++i){
       double s = values.get(s_pids[i]);
       result *= (1 + fabs(s) * (s > 0.0?r_plus[i]:r_minus[i]));
       if(result <= 0.0) return 0.0;
    }
    for(size_t i=0; i<f_pids.size(); ++i){
       result *= values.get(f_pids[i]);
       if(result <= 0.0) return 0.0;
    }
    return result;
}

REGISTER_PLUGIN(sys_rate_function)
