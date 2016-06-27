#include "plugins/minimizer_chain.hpp"
#include "interface/plugin.hpp"
#include "interface/phys.hpp"

using namespace theta;
using namespace std;

MinimizationResult minimizer_chain::minimize(const Function & f, const ParValues & start,
                const ParValues & step, const Ranges & ranges){
    MinimizationResult res;
    bool success = false;
    for(size_t i=0; i < minimizers.size(); ++i){
        try{
            res = minimizers[i].minimize(f, start, step, ranges);
            success = true;
        }
        catch(MinimizationException & ex){
            // if this was the last attempt: re-throw, otherwise silently ignore and try the next minimizer ...
            if(i+1==minimizers.size()) throw;
        }
        catch(std::logic_error & ex){
            stringstream ss;
            ss << ex.what();
            ss << " (in minimizer_chain, minimizer " << i << ")";
            throw logic_error(ss.str());
        }
        if(success) break;
    }
    if(last_minimizer.get()){
        ParValues step2 = step;
        // set step2 to the errors from the minimization, if available:
        const ParIds & pids = f.get_parameters();
        for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
            if(res.errors_plus.contains(*it)){
                double width = res.errors_plus.get(*it);
                if(width > 0){
                    step2.set(*it, width);
                }
            }
        }
        res = last_minimizer->minimize(f, res.values, step2, ranges);
    }
    return res;
}


minimizer_chain::minimizer_chain(const theta::Configuration & cfg){
    Setting s = cfg.setting;
    const size_t n = s["minimizers"].size();
    minimizers.reserve(n);
    size_t n_minimizers = 0;
    for(size_t i=0; i<n; ++i){
        minimizers.push_back(PluginManager<Minimizer>::build(Configuration(cfg, s["minimizers"][i])));
        ++n_minimizers;
    }
    if(s.exists("last_minimizer")){
        last_minimizer = PluginManager<Minimizer>::build(Configuration(cfg, s["last_minimizer"]));
        ++n_minimizers;
    }
    if(n_minimizers==0) throw ConfigurationException("no minimizers specified; required is at least one (counting last_minimizer)");
}

REGISTER_PLUGIN(minimizer_chain)

