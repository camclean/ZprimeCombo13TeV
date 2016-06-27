#include "interface/phys.hpp"
#include "interface/utils.hpp"
#include "test/utils.hpp"

#include <iostream>
#include <boost/test/unit_test.hpp>

using namespace theta;
using namespace std;

BOOST_AUTO_TEST_SUITE(functions)

BOOST_AUTO_TEST_CASE(add){
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    ParId p1 = vm->create_par_id("p1");
    ParId p2 = vm->create_par_id("p2");
    ParId p3 = vm->create_par_id("p3");
    load_core_plugins();
    ConfigCreator cc("f = {type = \"add\"; addends = (\"p1\", \"p2\", 1.7, 3.2, \"@f2\"); };\n"
                     "f2 = {type = \"add\"; addends = (\"p3\"); };"
            , vm);
    const theta::Configuration & cfg = cc.get();
    std::auto_ptr<Function> f_ptr = PluginManager<Function>::build(Configuration(cfg, cfg.setting["f"]));
    const Function & f = *f_ptr;
    ParValues values;
    values.set(p1, 0.0);
    values.set(p2, 0.0);
    values.set(p3, 0.0);
    BOOST_REQUIRE(values.contains_all(f.get_parameters()));
    double expected = 1.7 + 3.2;
    BOOST_CHECK(close_to_relative(f(values),expected));
    values.set(p1, 24.7);
    values.set(p2, -25.2);
    values.set(p3, 42.0);
    expected = (1.7 + 3.2) + 24.7 - 25.2 + 42.0;
    BOOST_CHECK(close_to_relative(f(values), expected));
}

BOOST_AUTO_TEST_CASE(test_exp){
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    ParId p1 = vm->create_par_id("p1");
    ParId p2 = vm->create_par_id("p2");
    ParId p3 = vm->create_par_id("p3");
    load_core_plugins();
    ConfigCreator cc("f = {type = \"exp_function\"; parameters = (\"p1\", \"p2\", \"p3\"); lambdas_plus = (0.1, 0.2, 0.3); lambdas_minus = (0.1, 0.19, 0.3); };\n" , vm);
    const theta::Configuration & cfg = cc.get();
    std::auto_ptr<Function> fp = PluginManager<Function>::build(Configuration(cfg, cfg.setting["f"]));
    const Function & f = *fp;
    BOOST_CHECK(f.get_parameters().contains(p1));
    BOOST_CHECK(f.get_parameters().contains(p2));
    BOOST_CHECK(f.get_parameters().contains(p3));
    BOOST_CHECK_EQUAL(f.get_parameters().size(), 3);
    ParValues values;
    values.set(p1, 1.0).set(p2, 0.0).set(p3, 0.0);
    double fval = f(values);
    double fval_exp = theta::utils::exp(0.1 * 1.0);
    BOOST_CHECK_EQUAL(fval_exp, fval);
    values.set(p1, -1.7);
    fval = f(values);
    fval_exp = exp(0.1 * -1.7);
    BOOST_CHECK_EQUAL(fval_exp, fval);
    values.set(p2, -0.9).set(p3, -1.2);
    fval = f(values);
    fval_exp = theta::utils::exp(0.1 * -1.7 - 0.9 * 0.19 - 1.2 * 0.3);
    BOOST_CHECK_EQUAL(fval_exp, fval);
    
    
    values.set(p1, 1.0).set(p2, 1.0).set(p3, 1.0);
    ParValues der;
    fval = f.eval_with_derivative(values, der);
    fval_exp = theta::utils::exp(0.1 + 0.2 + 0.3);
    BOOST_CHECK_EQUAL(fval_exp, fval);
    BOOST_REQUIRE(der.contains_all(f.get_parameters()));
    BOOST_CHECK_EQUAL(der.get(p1), 0.1 * fval_exp);
    BOOST_CHECK_EQUAL(der.get(p2), 0.2 * fval_exp);
    BOOST_CHECK_EQUAL(der.get(p3), 0.3 * fval_exp);
    
    values.set(p1, -1.0).set(p2, -1.0).set(p3, -1.0);
    fval = f.eval_with_derivative(values, der);   
    fval_exp = theta::utils::exp(-0.1 - 0.19 - 0.3);
    BOOST_CHECK_EQUAL(fval_exp, fval);
    BOOST_REQUIRE(der.contains_all(f.get_parameters()));
    BOOST_CHECK_EQUAL(der.get(p1), 0.1 * fval_exp);
    BOOST_CHECK_EQUAL(der.get(p2), 0.19 * fval_exp);
    BOOST_CHECK_EQUAL(der.get(p3), 0.3 * fval_exp);
}

BOOST_AUTO_TEST_SUITE_END()

