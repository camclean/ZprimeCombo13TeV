#ifndef CFG_UTILS_HPP
#define CFG_UTILS_HPP

#include "interface/decls.hpp"
#include "libconfig/libconfig.h++"

#include <string>
#include <set>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_map.hpp>


namespace theta{
    
    class SettingImplementation;
   
    /** \brief A class to keep record of the used configuration parameters
     * 
     * In order to warn users about unused configuration file flags,
     * plugins should call call the SettingUsageRecorder::markAsUsed methods.
     *
     * The user is warned by the %theta main program about unused warnings based on the
     * recorded usage using this class.
     */
    class SettingUsageRecorder{
    public:
        /** \brief Mark the setting at the given path as used
         */
        void markAsUsed(const std::string & path);
        
        /** \brief Get a list of paths of unused settings
         *
         * Searches recursively for unused settings.
         *
         * \param[out] unused The unused paths will be stored here
         * \param[in] aggregate_setting The setting to start the recursion; usually the root setting of the configuration file
         */
        void get_unused(std::vector<std::string> & unused, const Setting & aggregate_setting) const;

    private:
        std::set<std::string> used_paths;
    };
    
    /** \brief A Setting in a configuration file
     * 
     * A configuration file is a hierarchy of Settings. Each Setting is either a
     * primitive type (boolean, int float, string) or a compund type of key/value pairs (group)
     * or a list.
     * 
     * The Setting provides a read-only view on a setting/setting hierarchy.
     *
     * A setting is considered <em>used</em> one of the following methods are called on it:
     * <ul>
     *  <li>a casting operator</li>
     *  <li>access to a member in a compound with operator[]</li>
     *  <li>calling size()</li>
     * </ul>
     * The last definition is required to treat empty lists and groups correctly as used.
     * 
     * Note that the actual implementation is done in derived classes of SettingImplementation which actually reflects the hierarchally
     * structured, typed data. This class takes care of link resolving (with the "@..." notation) and the tracking
     * of which settings have been used. It also converts library-specific exceptions of cases in which a type is not found or type mismatches
     * to theta::ConfigurationException.
     */
    class Setting{
    public:
        
        enum Type{
            TypeBoolean, TypeInt, TypeFloat, TypeString, TypeGroup, TypeList
        };
        
        /// return a string representation of the Type enum (for error messages and debugging)
        static std::string type_to_string(const Type & type);
        
        ///@{
        /** \brief Convert the current setting to the given type
         *
         * If the setting has not the correct type, a ConfugurationException is thrown.
         */
        operator bool() const;
        operator std::string() const;
        operator double() const;
        operator long int() const;
        operator int() const;
        operator unsigned int() const;
        ///@}
        
       /** \brief Return a double, but allow the special strings "inf", "-inf" for infinity
        *
        * At some places in the configuration, it is allowed to use "inf" or "-inf" instead of a double,
        * and some imeplmentations can't express that natively as double type.
        * This function is used to parse these settings.
        */
        double get_double_or_inf() const;
        
        
       /** \brief Get number of entries in a setting group, list or array
        *
        * Returns 0 if the setting is not a group, list or array, and the number of
        * sub-settings otherwise.
        */
        size_t size() const;
        
       /** \brief Get entry by index in a setting group, list or array
        *
        * same as libconfig::Setting::operator[](int)
        */
        Setting operator[](int i) const;
        
       /** \brief Get entry by name in a setting group
        *
        * same as libconfig::Setting::operator[](string), with the only difference
        * that links of the form "@&lt;path&gt;" are followed
        */
        ///@{
        Setting operator[](const std::string & name) const;
        
        //This is defined in addition to operator[](string). Otherwise,
        // an expression in Setting s; s["abc"]; could also be read as
        // (int)s ["abc"]  (which is the same as "abc"[(int)s]) and the compiler would not
        // have a clear priority ...
        Setting operator[](const char * name) const{
            return operator[](std::string(name));
        }
        ///@}        
        
        /** \brief Returns whether a setting of this key exists in the current setting group
        */
        bool exists(const std::string & key) const;
        
        /** \brief Returns the name of the current setting within its parent setting
        *
        * In case the setting has no name, the string "&lt;noname&gt;" is returned.
        */
        std::string get_name() const;
        
        /** \brief Returns the configuration path of the current setting
         *
         * The path is a string built of all keys of the settings in the current hierarchy, separated by dots.
         */
        std::string get_path() const;
        
        /** \brief Returns the type of the setting
         */
        Type get_type() const;
        
        /** construct from a root setting.
         * 
         * Note that ownership of impl_root transferred.
         */
        Setting(std::auto_ptr<SettingImplementation> & impl_root, const boost::shared_ptr<SettingUsageRecorder> & rec);
        
    private:
        // memory is managed by the toplevel SettingImplementation (e.g. hierarchally, but that's up to the implementation)
        // and this Setting class takes care of deleting impl_root. This means that the destructor of this class
        // should call delete on impl_root if and only if:
        // 1. impl_this == impl_root   and
        // 2. the current instance is the last one referencing to this SettingImplementation
        // While 1. can be easily checked, 2. is handled by putting the root setting into a shared_ptr container, but only if this is the toplevel setting.
        const SettingImplementation * impl_this;
        const SettingImplementation * impl_root;
        boost::shared_ptr<SettingImplementation> sp_impl_root;
        
        // save the orginal name in case of links:
        std::string original_name;
        
        boost::shared_ptr<SettingUsageRecorder> rec;
        
        static Setting resolve_link(const SettingImplementation * setting, const SettingImplementation * root, const boost::shared_ptr<SettingUsageRecorder> & rec);
        Setting(const SettingImplementation * impl_this, const SettingImplementation * impl_root, const boost::shared_ptr<SettingUsageRecorder> & rec, const std::string & original_name_);
    };    
    
    /// The methods are the same as in Setting.
    /** \brief Base class for actual implementations of Settings
     * 
     * Not copyable; always passed around pointers. Note that Config manages the memory of a
     * SettingImplementation; derived classes should make sure that direct construction is not possible.
     */
    class SettingImplementation{
    public:
        virtual operator bool() const = 0;
        virtual operator std::string() const = 0;
        // note: implementations imlement long int only. Setting will convert it to int and unsigned int as required and perform range checking.
        virtual operator long int() const = 0;
        virtual operator double() const = 0;
        /// size is always 0 for types other than TypeList / TypeGroup
        virtual size_t size() const = 0;
        
        /// throws exception if type is not TypeList
        virtual const SettingImplementation & operator[](int i) const = 0;
        
        /// throws exception if type is not TypeGroup
        virtual const SettingImplementation & operator[](const std::string & name) const = 0;
        
        /** return whether a setting with the name exists in the current TypeGroup.
         * 
         * throws exception if type is not TypeGroup
         */
        virtual bool exists(const std::string & path) const = 0;
        
        /** return the name of the current setting (the key)
         * 
         * empty string for the top-level setting or for lists/arrays
         */        
        virtual std::string get_name() const = 0;
        
        /** return the path, i.e., dot-separated names in the current hierarchy of keys
         */
        virtual std::string get_path() const = 0;
        
        /// return the type of this setting
        virtual Setting::Type get_type() const = 0;
        
        virtual ~SettingImplementation();
    protected:
        SettingImplementation(){}
    private:
        SettingImplementation(const SettingImplementation & rhs); // not implemented
    };    
    
    class LibconfigSetting: public SettingImplementation {
        const libconfig::Setting * setting;
        std::auto_ptr<libconfig::Config> cfg; // filled only if this is the toplevel setting
        mutable boost::ptr_map<libconfig::Setting*, LibconfigSetting> children;
        
        LibconfigSetting(const LibconfigSetting &); // not implemented

        explicit LibconfigSetting(const libconfig::Setting * setting_);
    public:
        LibconfigSetting();
        
        static Setting parse(const std::string & cfg_string, const boost::shared_ptr<SettingUsageRecorder> & rec);
        virtual operator bool() const;
        virtual operator std::string() const;
        virtual operator long int() const;
        virtual operator double() const;
        virtual size_t size() const;
        virtual const SettingImplementation & operator[](int i) const;
        virtual const SettingImplementation & operator[](const std::string & name) const;
        virtual bool exists(const std::string & path) const;
        virtual std::string get_name() const;
        virtual std::string get_path() const;
        virtual Setting::Type get_type() const;
        virtual ~LibconfigSetting();
    };
}

#endif
