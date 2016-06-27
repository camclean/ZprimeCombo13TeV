#ifndef PLUGIN_TCC
#define PLUGIN_TCC

// This is the file implementing the class templates from plugin.hpp.
// It needs to be included only by those .cpp files which want to make a template
// instantiation via REGISTER_PLUGIN_BASETYPE, while most other only need to include
// plugin.hpp (including those calling REGISTER_PLUGIN).
// Splitting this section from the header file saves some compile time and space.

#include "interface/plugin.hpp"
#include "interface/exception.hpp"

//we need to make explicit template instantiations of the PluginManager registry,
// otherwise, the central registry associated to the PluginManager does not work reliably.
// This can happen if theta core does not instantiate the PluginManager class for a given type but
// some plugin does. In this case, the code behaves as if there are separate PluginManager
// classes for the same type and plugins are not found ...
//
// It seems to be a linker issue, and is platform-dependent. It does not affect modern
// Ubuntu distributions, but e.g. Mac OS X.
#define REGISTER_PLUGIN_BASETYPE(type) template class theta::PluginManager<type>;

namespace{
    //to increase the plugin_build_depth in an exception-safe manner, use the build_depth_sentinel,
    // which decrements depth count at destruction:
    struct plugin_build_depth_sentinel{
        int & i;
        plugin_build_depth_sentinel(int & i_):i(i_){
            ++i;
            if(i > 30){
                throw std::overflow_error("PluginManager::build: detected too deep plugin building (circular reference?)");
            }
        }
        ~plugin_build_depth_sentinel(){ --i; }
    };
}

namespace theta{

template<typename product_type>
PluginManager<product_type> & PluginManager<product_type>::instance(){
    static PluginManager pm;
    return pm;
}

template<typename product_type>
PluginManager<product_type>::PluginManager(): plugin_build_depth(0) {}

template<typename product_type>
std::auto_ptr<product_type> PluginManager<product_type>::build(const Configuration & ctx){
    PluginManager & pm = instance();
    plugin_build_depth_sentinel b(pm.plugin_build_depth);
    std::string type;
    if(!ctx.setting.exists("type")) type = "default";
    else type = static_cast<std::string>(ctx.setting["type"]);
    if(type=="") throw ConfigurationException("Error while constructing plugin: empty 'type' setting given in path '" + ctx.setting.get_path() + "'");
    std::auto_ptr<product_type> result;
    if(pm.pb.get()){
        result = pm.pb->build(ctx, type);
    }
    else{
        result = build_type(ctx, type);
    }
    return result;
}

template<typename product_type>
std::auto_ptr<product_type> PluginManager<product_type>::build_type(const Configuration & ctx, const std::string & type){
    PluginManager & pm = instance();
    std::auto_ptr<product_type> result;
    plugin_build_depth_sentinel b(pm.plugin_build_depth);
    if(type=="") throw std::invalid_argument("empty type");
    if(result.get()==0){
       for (size_t i = 0; i < pm.factories.size(); ++i) {
           if (pm.factories[i]->get_typename() != type) continue;
           try {
               result = pm.factories[i]->build(ctx);
           }catch (theta::Exception & ex) {
               std::stringstream ss;
               ss << "Error while constructing plugin in path '" << ctx.setting.get_path()
                   << "' (type='" << type << "', base_type = '" << demangle(typeid(product_type).name()) << "'): " << ex.message;
               throw ConfigurationException(ss.str());
           }
           catch (std::logic_error & ex) {
               std::stringstream ss;
               ss << "Error while constructing plugin in path '" << ctx.setting.get_path()
                   << "' (type='" << type << "', base_type = '" << demangle(typeid(product_type).name()) << "') logic error: " << ex.what();
               throw ConfigurationException(ss.str());
           }
       }
    }
    if(result.get()==0){
       std::stringstream ss;
       ss << "Error while constructing plugin according to configuration path '" << ctx.setting.get_path()
          << "': no plugin registered to create type='" << type << "', base_type = '" << demangle(typeid(product_type).name()) <<
          "'. Check spelling of the type and make sure to load all necessary plugin files via the setting 'options.plugin_files'";
       throw ConfigurationException(ss.str());
    }
    return result;
}

template<typename product_type>    
void PluginManager<product_type>::set_plugin_builder(std::auto_ptr<PluginBuilder<product_type> > & b){
    PluginManager & pm = instance();
    pm.pb = b;
}
    

template<typename product_type>
void PluginManager<product_type>::reset_plugin_builder(){
    PluginManager & pm = instance();
    pm.pb.reset();
}

template<typename product_type>
void PluginManager<product_type>::register_factory(factory_type * new_factory) {
    PluginManager & pm = instance();
    for (size_t i = 0; i < pm.factories.size(); i++) {
        if (pm.factories[i]->get_typename() == new_factory->get_typename()) {
            std::stringstream ss;
            ss << "PluginManager::register_factory: there is already a plugin registered for type '" << pm.factories[i]->get_typename() << "'";
            throw Exception(ss.str());
        }
    }
    pm.factories.push_back(new_factory);
}

}

#endif

