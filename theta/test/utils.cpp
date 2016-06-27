#include "test/utils.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"
#include "interface/variables.hpp"

#include <string>
#include <iostream>
#include <boost/test/unit_test.hpp>

using namespace theta;
using namespace std;

bool histos_equal(const Histogram1D & h1, const Histogram1D & h2){
    if(h1.get_nbins()!=h2.get_nbins()) return false;
    if(h1.get_xmin()!=h2.get_xmin()) return false;
    if(h1.get_xmax()!=h2.get_xmax()) return false;
    const size_t n = h1.get_nbins();
    for(size_t i=0; i<n; i++){
        if(h1.get(i)!=h2.get(i)) return false;
    }
    return true;
}


ConfigCreator::ConfigCreator(const std::string & cfg_string, const boost::shared_ptr<theta::VarIdManager> & vm):
      rec(new SettingUsageRecorder()), cfg(setup_cfg(cfg_string)){
    cfg.pm->set("default", vm);
}

Configuration ConfigCreator::setup_cfg(const string & cfg_string){
    try{
        Setting root = LibconfigSetting::parse(cfg_string, rec);
        return Configuration(root);
    }
    catch(Exception & ex){
        std::cerr << "ConfigCreator: " << ex.what() << endl;
        throw;
    }
}

void load_core_plugins(){
    /*static bool loaded(false);
    if(loaded) return;
    BOOST_TEST_CHECKPOINT("loading core plugin");
    try{
        PluginLoader::load("lib/core-plugins.so");
    }
    catch(exception & ex){
      std::cout << "std::exception in load_core_plugins: " << ex.what() << std::endl;
      throw;
    }
    BOOST_TEST_CHECKPOINT("loaded core plugin");
    loaded = true;*/
}

bool load_root_plugins(){
    static bool loaded(false);
    if(loaded) return true;
    BOOST_TEST_CHECKPOINT("loading root plugin");
    try{
        PluginLoader::load("lib/root.so");
    }
    catch(exception & ex){
        return false;
    }
    BOOST_TEST_CHECKPOINT("loaded root plugin");
    loaded = true;
    return true;
}

bool load_llvm_plugins(){
    static bool loaded(false);
    static bool load_err(false);
    if(loaded) return true;
    if(load_err) return false;
    BOOST_TEST_CHECKPOINT("loading llvm plugins");
    try{
        PluginLoader::load("lib/llvm-plugins.so");
    }
    catch(exception & ex){
        std::cout << "error loadin llvm plugins: " << ex.what() << endl;
    	load_err = true;
        return false;
    }
    BOOST_TEST_CHECKPOINT("loaded llvm plugin");
    loaded = true;
    return true;
}


