#include "interface/plugin.hpp"
#include "interface/utils.hpp"
#include "interface/variables.hpp"

#include <dlfcn.h>

using namespace std;
using namespace theta;

Configuration::Configuration(const Setting & setting_): pm(new theta::PropertyMap()), setting(setting_){}

Configuration::Configuration(const Configuration & cfg, const Setting & setting_): pm(cfg.pm), setting(setting_){}

void PluginLoader::execute(const Configuration & cfg) {
    Setting files = cfg.setting["plugin_files"];
    size_t n = files.size();
    for (size_t i = 0; i < n; i++) {
        std::string filename = utils::replace_theta_dir(files[i]);
        load(filename);
    }
}

void PluginLoader::load(const std::string & soname) {
    void* handle = 0;
    try {
        handle = dlopen(soname.c_str(), RTLD_NOW | RTLD_GLOBAL);
    } catch (std::exception & ex) {
        std::stringstream ss;
        ss << ex.what() << " (in PluginLoader::load while loading plugin file '" << soname << "')";
        throw Exception(ss.str());
    }
    if (handle == 0) {
        std::stringstream s;
        const char * error = dlerror();
        if (error == 0) error = "0";
        s << "PluginLoader::load: error loading plugin file '" << soname << "': " << error << std::endl;
        throw Exception(s.str());
    }
}

theta::ParId theta::get_parameter(const Configuration & cfg, const std::string & key){
    return cfg.pm->get<VarIdManager>()->get_par_id(cfg.setting[key]);
}
