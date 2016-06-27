#include "interface/phys.hpp"
#include "test/utils.hpp"

#include <iostream>
#include <boost/test/unit_test.hpp>

using namespace theta;
using namespace std;

BOOST_AUTO_TEST_SUITE(nl_gauss)

BOOST_AUTO_TEST_CASE(basic){
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    ParId p1 = vm->create_par_id("p1");
    ParId p2 = vm->create_par_id("p2");
    load_core_plugins();
    ConfigCreator cc("f = {type = \"nl_gauss\"; mu = [1.0, 2.0]; rows = (\"p1\", \"p2\"); covariance = \n"
            "   ([0.5, 0.0], [0.0, 0.3]); };"
            , vm);
    const theta::Configuration & cfg = cc.get();
    std::auto_ptr<Function> f_ptr;
    try{
        f_ptr = PluginManager<Function>::build(Configuration(cfg, cfg.setting["f"]));
    }
    catch(ConfigurationException & ex){
        cout << ex.message << endl;
        BOOST_CHECK(false);
        return;
    }
    const Function & f = *f_ptr;
    ParValues values;
    values.set(p1, 0.0);
    values.set(p2, 0.0);
    BOOST_CHECK(values.contains_all(f.get_parameters()));
    double expected = 0.5 * (pow(1.0, 2) / 0.5 + pow(2.0, 2) / 0.3);
    BOOST_CHECK(close_to_relative(f(values),expected));
    values.set(p1, 2.7);
    values.set(p2, -3.2);
    expected = 0.5 * (pow(2.7 - 1.0, 2)  / 0.5 + pow(2.0 + 3.2, 2) / 0.3);
    BOOST_CHECK(close_to_relative(f(values), expected));
}

BOOST_AUTO_TEST_SUITE_END()

