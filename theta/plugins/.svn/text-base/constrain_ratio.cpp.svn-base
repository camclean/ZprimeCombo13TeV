#include "plugins/constrain_ratio.hpp"
#include "interface/exception.hpp"
#include "interface/plugin.hpp"

using namespace std;
using namespace theta;

constrain_ratio::constrain_ratio(const theta::Configuration & cfg): pid_nominator(cfg.pm->get<VarIdManager>()->get_par_id(cfg.setting["nominator"])),
   pid_denominator(cfg.pm->get<VarIdManager>()->get_par_id(cfg.setting["denominator"])){
    mean = cfg.setting["mean"];
    width = cfg.setting["width"];
    if(width <= 0.0) throw theta::ConfigurationException("width must be > 0");
    par_ids.insert(pid_nominator);
    par_ids.insert(pid_denominator);
}

double constrain_ratio::operator()(const theta::ParValues & values) const{
    double nominator = values.get(pid_nominator);
    double denominator = values.get(pid_denominator);
    if(denominator==0.0) throw invalid_argument("constrain_ratio: zero denominator");
    return 0.5 * pow((nominator / denominator - mean) / width, 2);
}

REGISTER_PLUGIN(constrain_ratio)
