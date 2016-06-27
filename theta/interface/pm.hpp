#ifndef PM_HPP
#define PM_HPP

#include <boost/shared_ptr.hpp>
#include <map>
#include <string>
#include <typeinfo>

#include <sstream>

#include "interface/exception.hpp"

namespace theta{

/** \brief Property map mapping strings to shared_ptr<T> where T is an arbitrary type
 *
 * Main use of this class is as member of theta::Configuration. See there for standard
 * properties plugins can access.
 * 
 * Elements in the container are indexed by both a supplied instance name and the type
 * passed to the get and set function templates, i.e., both the name and the type have to be the same.
 *
 * Setting a property to a NULL pointer effectively deletes the property.
 */
class PropertyMap{
public:

   /** \brief Retrieve a value previously set
    *
    * Throws a std::range_error if no value is stored with this combination of instance name and type T.
    */
   template<typename T>
   boost::shared_ptr<T> get(const std::string & instance_name = "default") const;

   /** \brief Store a value
    *
    * Stores (and possibly overwrites) a value using the instance_name and given type T.
    */
   template<typename T>
   void set(const std::string & instance_name, const boost::shared_ptr<T> & value);

   /** Get the use count os the given instance.
    *
    * Returns 0 if the instance is not part of the property map. For existing instances,
    * the value is at least 1, the reference from the PropertyMap.
    */
   template<typename T>
   size_t use_count(const std::string & instance_name = "default") const;
   
   /** Returns whether the property indicated by type T and instance_name exists, i.e. is non-null
    *
    */
   template<typename T>
   bool exists(const std::string & instance_name = "default") const;
   
private:
   struct nametype{
       std::string name;
       const std::type_info & type;
       nametype(const std::string & name_, const std::type_info & type_):name(name_), type(type_){}
       bool operator<(const nametype & rhs) const{
           return (name < rhs.name) || (name == rhs.name && type.before(rhs.type));
       }
   };
   std::map<nametype, boost::shared_ptr<void> > instances;
};


template<typename T>
boost::shared_ptr<T> PropertyMap::get(const std::string & instance_name) const {
    nametype nt(instance_name, typeid(T));
    std::map<nametype, boost::shared_ptr<void> >::const_iterator it=instances.find(nt);
    if(it==instances.end()){
        std::stringstream ss;
        ss << "PropertyMap: instance with name '" << instance_name << "' with type " << typeid(T).name() << " not found";
        throw std::range_error(ss.str());
    }
    return boost::static_pointer_cast<T>(it->second);
}

template<typename T>
size_t PropertyMap::use_count(const std::string & instance_name) const {
    nametype nt(instance_name, typeid(T));
    std::map<nametype, boost::shared_ptr<void> >::const_iterator it=instances.find(nt);
    if(it==instances.end()) return 0;
    return it->second.use_count();
}

template<typename T>
bool PropertyMap::exists(const std::string & instance_name) const {
    nametype nt(instance_name, typeid(T));
    return instances.find(nt) != instances.end();
}

template<typename T>
void PropertyMap::set(const std::string & instance_name, const boost::shared_ptr<T> & value){
    nametype nt(instance_name, typeid(T));
    if(value.get()) instances[nt] = value;
    else instances.erase(nt);
}

}

#endif

