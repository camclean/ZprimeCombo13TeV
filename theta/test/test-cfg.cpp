#include "test/utils.hpp"
#include "interface/variables.hpp"
#include "interface/exception.hpp"

#include <boost/test/unit_test.hpp>

using namespace theta;

// Test some of the features of the configuratio system, most notably the "@" indirection

BOOST_AUTO_TEST_SUITE(testcfg)

BOOST_AUTO_TEST_CASE(cfg){
    boost::shared_ptr<VarIdManager> vm(new VarIdManager);
    PropertyMap pm;
    pm.set("default", vm);
    ConfigCreator cc(
            "setting1 = {s1k1 = 1.0; s1k2 = 2; s1k3 = \"@value3\"; s1k4 = \"@nested1.nested2.value4\"; };\n"
            "value3 = 3;\n"
            "nested1 = { nested2 = { value4 = 4.0; }; };\n"
            "\n"
            "setting2 = \"@setting3\";\n"
            "setting3 = \"@setting2\";\n"
            "\n"
            "setting4 = [1.0, 2.0, 3.0];\n"
            , vm);
    const theta::Configuration & cfg = cc.get();
    
    // some basic stuff: readin name, type and value:
    Setting s1 = cfg.setting;
    BOOST_CHECK(s1.exists("setting1"));
    Setting s2 = s1["setting1"]["s1k1"];
    BOOST_CHECK(s2.get_name() == "s1k1");
    BOOST_CHECK(s2.get_type() == Setting::TypeFloat);
    double dvalue = s2;
    BOOST_CHECK(dvalue == 1.0);
    
    // same for s1k2:
    s2 = s1["setting1"]["s1k2"];
    BOOST_CHECK(s2.get_type() == Setting::TypeInt);
    int ivalue = static_cast<int>(s2);
    BOOST_CHECK(ivalue == 2);
    // int and double are not the same:
    bool exception = false;
    try{
        dvalue = s2;
    }
    catch(const ConfigurationException &){
        exception = true;
    }
    BOOST_CHECK(exception);
    
    // basic indirection:
    s2 = s1["setting1"]["s1k3"];
    ivalue = s2;
    BOOST_CHECK(ivalue == 3);
    // we should keep the original name in indirection:
    BOOST_CHECK(s2.get_name() == "s1k3");
    
    // multiple levels of indirection:
    s2 = s1["setting1"]["s1k4"];
    dvalue = s2;
    BOOST_CHECK(dvalue == 4.0);
    BOOST_CHECK(s2.get_name() == "s1k4");
    
    // endless redirections:
    exception = false;
    try{
        s2 = s1["setting2"];
    }
    catch(ConfigurationException & ){
        exception = true;
    }
    BOOST_CHECK(exception);
    
    // looping through settings and reading the names:
    s2 = s1["setting1"];
    BOOST_REQUIRE(s2.size() == 4);
    BOOST_REQUIRE(s2[0].get_name() == "s1k1");
    BOOST_REQUIRE(s2[1].get_name() == "s1k2");
    BOOST_REQUIRE(s2[2].get_name() == "s1k3");
    BOOST_REQUIRE(s2[3].get_name() == "s1k4");
    
    // arrays and names:
    s2 = s1["setting4"];
    BOOST_CHECK(s2.get_type() == Setting::TypeList);
    BOOST_REQUIRE(s2.size() == 3);
    BOOST_CHECK(s2[0].get_name() == "");
    BOOST_CHECK(s2[2].get_name() == "");
    dvalue = s2[0];
    BOOST_CHECK(dvalue = 1.0);
}

BOOST_AUTO_TEST_SUITE_END()
