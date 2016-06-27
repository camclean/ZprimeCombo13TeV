#include "interface/plugin.hpp"
#include "interface/cfg-utils.hpp"
#include "interface/phys.hpp"
#include "interface/minimizer.hpp"
#include "test/utils.hpp"

#include <iostream>

#include <boost/test/unit_test.hpp> 

using namespace std;
using namespace theta;
using namespace libconfig;

using namespace theta::utils;

BOOST_AUTO_TEST_SUITE(minimizer_tests)


class ImpossibleFunction: public Function{
public:
    ImpossibleFunction(const ParIds & pids){
        par_ids = pids;
    }

    virtual double operator()(const ParValues & v) const{
        return 0.0;
    }

    virtual double gradient(const ParValues & v, const ParId & pid) const{
        return 0.0;
    }
};

BOOST_AUTO_TEST_CASE(minuit){
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    if(!load_root_plugins()){
      std::cout << "In test minuit: root plugin could not be loaded, not executing root tests";
      return;
    }
    
    ParId p0 = vm->create_par_id("p0");
    ParId p1 = vm->create_par_id("p1");
    ParIds pars;
    pars.insert(p0);
    pars.insert(p1);
    
    ConfigCreator cc2("type = \"root_minuit\";", vm);
    BOOST_REQUIRE(true);//create checkpoint
    std::auto_ptr<Minimizer> min = PluginManager<Minimizer>::build(cc2.get());
    BOOST_REQUIRE(min.get());
    ImpossibleFunction f(pars);
    bool exception;
    ParValues start, step;
    std::map<ParId, pair<double, double> > ranges;
    start.set(p0, 0.0).set(p1, 0.0);
    step.set(p0, 1.0).set(p1, 1.0);
    ranges[p0] = make_pair(-100.0, 100.0);
    ranges[p1] = make_pair(-100.0, 100.0);
    try{
        BOOST_REQUIRE(true);//create checkpoint
        min->minimize(f, start, step, ranges);
    }
    catch(MinimizationException & ex){
        exception = true;
    }
    catch(std::exception & ex){ //should not happen ...
        std::cout << ex.what() << endl;
        BOOST_REQUIRE(false);
    }
    BOOST_REQUIRE(exception);
}

BOOST_AUTO_TEST_SUITE_END()

