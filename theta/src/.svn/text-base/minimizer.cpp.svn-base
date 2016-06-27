#include "interface/minimizer.hpp"
#include "interface/model.hpp"
#include "interface/asimov-utils.hpp"
#include "interface/plugin.tcc"

REGISTER_PLUGIN_BASETYPE(theta::Minimizer);

using namespace theta;
using namespace std;

void MinimizationResult::operator=(const MinimizationResult& rhs){
    fval = rhs.fval;
    values.set(rhs.values);
    errors_plus.set(rhs.errors_plus);
    errors_minus.set(rhs.errors_minus);
    covariance = rhs.covariance;
}

FunctionInfo::FunctionInfo(const ParValues & start_, const ParValues & step_, const Ranges & ranges_, const ParValues & fixed_parameters): start(start_), step(step_), ranges(ranges_){
    ParIds pids = start.get_parameters();
    for(ParIds::const_iterator pit=pids.begin(), it_end = pids.end(); pit!=it_end; ++pit){
        if(!step.contains(*pit)) throw invalid_argument("FunctionInfo: step does not contain all parameters from start");
        // 1. check whether parameter is fixed by step and range:
        const pair<double, double> & r = ranges.get(*pit);
        if(r.first==r.second){
            if(step.get_unchecked(*pit)>0.0){
                throw invalid_argument("FunctionInfo: inconsistent range/step given: range empty but step > 0");
            }
            fixed_parids.insert(*pit);
        }
        else{
            if(step.get_unchecked(*pit)<=0.0){
                throw invalid_argument("FunctionInfo: step <= 0.0 for non-empty range given");
            }
        }
        // 2. check whether parameter is fixed by fixed_parameters. Note that
        //  the value given in fixed_parameter overrides the value according to start.
        if(fixed_parameters.contains(*pit)){
            double val = fixed_parameters.get(*pit);
            start.set(*pit, val);
            ranges.set(*pit, make_pair(val, val));
            step.set(*pit, 0.0);
            fixed_parids.insert(*pit);
        }
    }
}

FunctionInfo::~FunctionInfo(){}

boost::shared_ptr<FunctionInfo> Minimizer::create_nll_function_info(const Model & model, const boost::shared_ptr<Distribution> & override_parameter_distribution, const ParValues & fixed_parameters){
    const Distribution & dist = override_parameter_distribution.get()? *override_parameter_distribution: model.get_parameter_distribution();
    ParValues start;
    dist.mode(start);
    Ranges ranges(dist);
    ParIds pids = fixed_parameters.get_parameters();
    for(ParIds::const_iterator pit = pids.begin(); pit!=pids.end(); ++pit){
        double val = fixed_parameters.get(*pit);
        start.set(*pit, val);
        ranges.set(*pit, make_pair(val, val));
    }
    ParValues step = asimov_likelihood_widths(model, override_parameter_distribution);
    return boost::shared_ptr<FunctionInfo>(new DefFunctionInfo(start, step, ranges, fixed_parameters));
}

MinimizationResult Minimizer::minimize2(const Function & f, const FunctionInfo & info, const ParValues & fixed_parameters){
    dynamic_cast<const DefFunctionInfo&>(info); // throws bad_cast
    ParIds pids = fixed_parameters.get_parameters();
    if(pids.size()==0){
        return minimize(f, info.get_start(), info.get_step(), info.get_ranges());
    }
    else{
        ParValues start(info.get_start());
        ParValues step(info.get_step());
        Ranges ranges(info.get_ranges());
        const ParIds & info_fixed = info.get_fixed_parameters();
        for(ParIds::const_iterator pit = pids.begin(); pit!=pids.end(); ++pit){
            if(!info_fixed.contains(*pit)){
                throw invalid_argument("fixed parameter in minimize2 which is not fixed in info. This is not allowed.");
            }
            double val = fixed_parameters.get(*pit);
            start.set(*pit, val);
            step.set(*pit, 0.0);
            ranges.set(*pit, make_pair(val, val));
        }
        return minimize(f, start, step, ranges);
    }
}

Minimizer::~Minimizer(){}

Minimizer::DefFunctionInfo::~DefFunctionInfo(){}

