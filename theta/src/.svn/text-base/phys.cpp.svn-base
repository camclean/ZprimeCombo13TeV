#include "interface/phys.hpp"
#include "interface/plugin.tcc"

using namespace theta;

REGISTER_PLUGIN_BASETYPE(Function);
REGISTER_PLUGIN_BASETYPE(DataSource);

Function::~Function(){}

void Function::fill_vpids() const {
    theta_assert(vpids.size()==0);
    vpids.reserve(par_ids.size());
    size_t i=0;
    for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it, ++i){
        vpids.push_back(*it);
    }
}

double Function::eval_with_derivative(const ParValues & v, ParValues & der) const{
    throw Exception("not implemented for " + demangle(typeid(*this).name()));
}

