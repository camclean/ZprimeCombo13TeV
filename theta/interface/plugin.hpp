#ifndef PLUGIN_HPP
#define PLUGIN_HPP

#include "interface/decls.hpp"
#include "interface/pm.hpp"
#include "interface/cfg-utils.hpp"

#include <boost/shared_ptr.hpp>

#include <typeinfo>
#include <sstream>

namespace theta {

    /** \brief A container class which is used to construct conrete types managed by the plugin system
     *
     * An instance of this class is passed to the constructor of classes managed by the plugin system.
     * It contains the required information for plugins to construct an instance of the plugin class,
     * the most important being \c setting, which is the setting group from the configuration file
     * for which this plugin class should be created.
     *
     * The standard properties in the PropertyMap are:
     * <ul>
     *   <li>VarIdManager "default" to get information about variables</li>
     *   <li>ProductsSink "default" where the producers  and other ProductSources write their results</li>
     *   <li>int "runid"</li>
     *   <li>RndInfoTable "default"</li>
     * </ul>
     */
    class Configuration{
    public:

        /// A property map giving access to various shared objects
        boost::shared_ptr<PropertyMap> pm;

        /// The setting in the configuration file from which to build the instance
        Setting setting;

        /** \brief Construct Configuration from the root setting
         */
        explicit Configuration(const Setting & root_setting);

        /** \brief Copy elements from another Configuration, but replace Configuration::setting
         *
         * Copy all from \c cfg but \c cfg.setting which is replaced by \c setting_.
         */
        Configuration(const Configuration & cfg, const Setting & setting_);
    };
    
    // convert the string cfg.setting[key] to a ParId, using cfg.pm->get<VarIdManager>().
    // This is useful for many producers, as ParId does not have a default constructor.
    theta::ParId get_parameter(const Configuration & cfg, const std::string & key);


    template<typename> class PluginManager;
    /** \brief Class used internally for the PluginManager
     *
     * You usually do not have to care about it, as this is a detail handled by the REGISTER_PLUGIN macro.
     *
     * This is the abstract factory class for a certain \c base_type. For each \c base_type, there
     * is an instance of PluginManager&lt;base_type&gt; which will save pointers to all currently registered
     * factories of type factory&lt;base_type&gt;
     *
     * By use of the REGISTER_PLUGIN(some_concrete_type), a derived class of factory&lt;some_concrete_type::base_type&gt;
     * will be created which handles the construction of \c some_concrete_type. This factory will be registered
     * at the PluginManager&lt;some_conrete_type::base_type&gt; as soon as the shared object file is loaded.
     */
     template<typename base_type>
     class factory{
     public:
         /// build an instance from a Configuration object
         virtual std::auto_ptr<base_type> build(const Configuration & cfg) = 0;

         /// the type of the object this factory is responsible for; it corresponds to the type="..." configuration file setting
         virtual std::string get_typename() = 0;

         virtual ~factory(){}
     protected:
         /// register this factory at the correct PluginManager
         void reg(){
             PluginManager<base_type>::register_factory(this);
         }
     };

     /* define a template specialization of abstract_factory<base_type> which constructs the desired type
     */
     #define REGISTER_PLUGIN_NAME(type,name) namespace { struct factory__##type##name : public theta::factory<type::base_type>{ \
     virtual std::auto_ptr<type::base_type> build(const theta::Configuration & cfg){return std::auto_ptr<type::base_type>(new type(cfg)); }\
     virtual std::string get_typename(){ return #name ;}\
     factory__##type##name (){reg();}\
     } factory_instance__##type##name;}

     #define REGISTER_PLUGIN(type) REGISTER_PLUGIN_NAME(type, type)
     #define REGISTER_PLUGIN_DEFAULT(type) REGISTER_PLUGIN_NAME(type, default) 
     
     template<typename T>
     class PluginBuilder{
     public:
         virtual std::auto_ptr<T> build(const Configuration & cfg, const std::string & type) = 0;
         virtual ~PluginBuilder(){}
     };

    /** \brief Central registry class for plugins.
     *
     * This class serves two purposes:
     * <ol>
     * <li>Build an instance from a plugin class using a theta::Configuration as only argument</li>
     * <li>Serve as central registry for all plugins of a certain type.</li>
     * </ol>
     * Registration in 2. should be done via the REGISTER_PLUGIN macro, so the user
     * usually does not have to care about the details of plugin registration.
     *
     * Building an instance of a certain plugin basis type T, use PluginManager<T>::build(const Configuration & ). This
     * returns an auto_ptr<T> which contains the requested instance or throws an exception in case of failure.
     */
    template<typename product_type>
    class PluginManager {
    public:

        /** \brief Use the registered factories to build an instance from a configuration settings block.
         *
         * This will go through all registered plugin and use the factory with the matching name and use this
         * factory to construct the instance.
         *
         * The lookup rules for finding out the name are:
         * <ol>
         *   <li>If using \c build_type, the name given in \c type ise used as plugin name directly.</li>
         *   <li>If using \c build, the string in cfg.setting["type"] is used if set; otherwise, use the string "default" is used.
         *      Using this string, the currently set PluginBuilder is called. Usually, this will eventually call \c build_type.
         * </li>
         * </ol>
         *
         * The default PluginBuilder just uses the typename as determined in 2. to lookup the plugin.
         */
        static std::auto_ptr<product_type> build(const Configuration & cfg);
        static std::auto_ptr<product_type> build_type(const Configuration & cfg, const std::string & type);

        static void set_plugin_builder(std::auto_ptr<PluginBuilder<product_type> > & b);
        
        // revert to default builder
        static void reset_plugin_builder();

    private:
        typedef typename product_type::base_type base_type;
        typedef factory<base_type> factory_type;
        friend class factory<base_type>;

        static PluginManager & instance();

        //prevent instance construction from "outside" by making constructor private:
        PluginManager();
        PluginManager(const PluginManager & rhs);// not implemented

        /** \brief Register a new factory.
         *
         * An InvalidArgumentException is thrown if there already is
         * a registered plugin for the same type string.
         *
         * Used by the REGISTER_PLUGIN macro. Not for direct call.
         */
        static void register_factory(factory_type * new_factory);

        std::vector<factory_type*> factories;
        std::auto_ptr<PluginBuilder<product_type> > pb;
        
        // to prevent endless recursion within PluginManager::build, use a depth counter:
        int plugin_build_depth;
    };

    /** \brief Class responsible to load the shared object files containing plugins
     */
    class PluginLoader {
    public:

        /** \brief Run the loader according to the configuration file setting
         *
         * The shared-object files in the \c plugins_files setting are loaded. This
         * setting must be a list of strings of the filenames of the .so files.
         */
        static void execute(const Configuration & cfg);

        /** \brief load a single plugin file
         *
         * The given shared object file will be loaded which will trigger the plugin registration of all plugins
         * defined via the REGISTER_PLUGIN macro automagically.
         *
         * \param soname is the filename of the plugin (a .so file), including the path from the current directory
         */
        static void load(const std::string & soname);
    };
}

#endif

