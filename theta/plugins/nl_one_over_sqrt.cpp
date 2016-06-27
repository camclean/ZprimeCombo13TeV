#include "plugins/nl_one_over_sqrt.hpp"
#include "interface/exception.hpp"
#include "interface/plugin.hpp"

using namespace std;
using namespace theta;

nl_one_over_sqrt::nl_one_over_sqrt(const theta::Configuration & cfg): pid(cfg.pm->get<VarIdManager>()->get_par_id(cfg.setting["parameter"])){
    par_ids.insert(pid);
}

double nl_one_over_sqrt::operator()(const theta::ParValues & values) const{
    double val = values.get(pid);
    if(val < 0.0) throw invalid_argument("nl_one_over_sqrt: negative argument");
    return 0.5 * val;
}

REGISTER_PLUGIN(nl_one_over_sqrt)
