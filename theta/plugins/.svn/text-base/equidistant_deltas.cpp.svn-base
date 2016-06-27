#include "plugins/equidistant_deltas.hpp"
#include "interface/exception.hpp"
#include "interface/random.hpp"

using namespace theta;
using namespace std;

equidistant_deltas::equidistant_deltas(const theta::Configuration & cfg){
    n = cfg.setting["n"];
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    if(n<2) throw ConfigurationException("n>=2 must hold!");
    par_ids.insert(vm->get_par_id(cfg.setting["parameter"]));
    low = support_.first = cfg.setting["range"][0];
    support_.second = cfg.setting["range"][1];
    width_ = (support_.second - support_.first) / (n-1);
}

void equidistant_deltas::sample(theta::ParValues & result, theta::Random & rnd) const{
    unsigned int i = rnd.get_uniform_int(n);
    result.set(*par_ids.begin(), low + i * width_);
}

void equidistant_deltas::mode(theta::ParValues & result) const{
    result.set(*par_ids.begin(), (support_.second - support_.first) / 2 + support_.first);
}

double equidistant_deltas::eval_nl(const theta::ParValues & values) const{
    return 0;
}

const std::pair<double, double> & equidistant_deltas::support(const theta::ParId&) const{
    return support_;
}

REGISTER_PLUGIN(equidistant_deltas)
