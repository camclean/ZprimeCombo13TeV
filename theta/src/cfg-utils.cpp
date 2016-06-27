#include "interface/cfg-utils.hpp"
#include "interface/exception.hpp"
#include <sstream>
#include <limits>

using namespace std;
using namespace theta;

namespace{
    std::string format_error(const Setting & s, const string & name, const string & message){
        return "While accessing setting at path '" + s.get_path() + "." + name + "': " + message;
    }
    
    std::string format_wrongtype(const Setting & s,  const Setting::Type & expected_type){
        return "expected type " + (Setting::type_to_string(expected_type)) + " for setting path '" + s.get_path() + "', but found type " + Setting::type_to_string(s.get_type());
    }
}

void SettingUsageRecorder::markAsUsed(const string & path){
    used_paths.insert(path);
}

void SettingUsageRecorder::get_unused(std::vector<std::string> & unused, const Setting & aggregate_setting) const{
    size_t n = aggregate_setting.size();
    for(size_t i=0; i<n; ++i){
        std::string path = aggregate_setting[i].get_path();
        bool a_unused = false;
        if(used_paths.find(path) == used_paths.end()){
            unused.push_back(path);
            a_unused = true;
        }
        //don't descend if already aggregate was reported as unused ...
        if((aggregate_setting[i].get_type() == Setting::TypeList || aggregate_setting[i].get_type() == Setting::TypeGroup) && !a_unused){
            get_unused(unused, aggregate_setting[i]);
        }
    }
}

/* Setting */
std::string Setting::type_to_string(const Type & type){
    switch(type){
        case TypeBoolean: return "boolean";
        case TypeInt: return "int";
        case TypeFloat: return "float";
        case TypeString: return "string";
        case TypeGroup: return "group";
        case TypeList: return "list";
    }
    throw logic_error("type_to_string internal error");
}

Setting::Setting(const SettingImplementation * impl_this_, const SettingImplementation * impl_root_, const boost::shared_ptr<SettingUsageRecorder> & rec_, const string & original_name_):
   impl_this(impl_this_), impl_root(impl_root_), original_name(original_name_), rec(rec_){}
   
Setting::Setting(std::auto_ptr<SettingImplementation> & impl_root_, const boost::shared_ptr<SettingUsageRecorder> & rec_): sp_impl_root(impl_root_), rec(rec_){
    impl_this = impl_root = sp_impl_root.get();
}

Setting::operator bool() const{
    if(rec)rec->markAsUsed(get_path());
    try{
        return *impl_this;
    }
    catch(...){
        throw ConfigurationException(format_wrongtype(*this, TypeBoolean));
    }
}

Setting::operator std::string() const{
    if(rec)rec->markAsUsed(get_path());
    try{
        return *impl_this;
    }
    catch(...){
        throw ConfigurationException(format_wrongtype(*this, TypeString));
    }
}

Setting::operator long int() const{
    if(rec)rec->markAsUsed(get_path());
    try{
        return *impl_this;
    }
    catch(...){
        throw ConfigurationException(format_wrongtype(*this, TypeInt));
    }
}

Setting::operator int() const{
    if(rec)rec->markAsUsed(get_path());
    try{
        long int res = static_cast<long int>(*impl_this);
        if(res > numeric_limits<long int>::max()){
            throw ConfigurationException("int overflow at setting path '" + get_path() + "'");
        }
        return res;
    }
    catch(...){
        throw ConfigurationException(format_wrongtype(*this, TypeInt));
    }
}

Setting::operator unsigned int() const{
    if(rec)rec->markAsUsed(get_path());
    try{
        long int res = static_cast<long int>(*impl_this);
        if(res > numeric_limits<unsigned int>::max() || res < 0){
            throw ConfigurationException("unsigned int overflow/underflow at setting path '" + get_path() + "'");
        }
        return res;
    }
    catch(...){
        throw ConfigurationException(format_wrongtype(*this, TypeInt));
    }
}

Setting::operator double() const{
    if(rec)rec->markAsUsed(get_path());
    try{
        return *impl_this;
    }
    catch(...){
        throw ConfigurationException(format_wrongtype(*this, TypeFloat));
    }
}

double Setting::get_double_or_inf() const {
    if(rec)rec->markAsUsed(get_path());
    if(get_type()==TypeFloat) return *this;
    string infstring = *this;
    if(infstring == "inf" || infstring == "+inf") return numeric_limits<double>::infinity();
    if(infstring == "-inf") return -numeric_limits<double>::infinity();
    throw ConfigurationException(format_wrongtype(*this, TypeFloat));
}

size_t Setting::size() const{
    if(rec)rec->markAsUsed(get_path());
    return impl_this->size();
}

Setting Setting::resolve_link(const SettingImplementation * setting, const SettingImplementation * root, const boost::shared_ptr<SettingUsageRecorder> & rec){
    string original_name = setting->get_name();
    try{
        std::string next_path;
        //hard-code maximum redirection level of 10:
        for(int i=0; i <= 10; ++i){
            const SettingImplementation * s = setting;
            if(i>0){
                s = root;
                do{
                    size_t dotpos = next_path.find('.');
                    if(dotpos==string::npos) dotpos = next_path.size();
                    s = &((*s)[next_path.substr(0, dotpos)]);
                    next_path.erase(0, dotpos+1);
                }while(next_path.size());
            }
            if(s->get_type() != Setting::TypeString) return Setting(s, root, rec, original_name);
            std::string link = *s;
            if(link.size()==0 || link[0]!='@'){
                return Setting(s, root, rec, original_name);
            }
            link.erase(0, 1);
            next_path = link;
            //mark any intermediate link as used:
            if(rec)rec->markAsUsed(s->get_path());
        }
    }
    catch(Exception & ex){
        std::stringstream ss;
        ss << "Exception while trying to resolve link at " << setting->get_path() << ": " << ex.message;
        ex.message = ss.str();
        throw;
    }
    std::stringstream ss;
    ss << "While trying to resolve link at " << setting->get_path() << ": link level is too deep";
    throw ConfigurationException(ss.str());
}


Setting Setting::operator[](int i) const{
    try{
        if(rec)rec->markAsUsed(get_path());
        return Setting::resolve_link(&((*impl_this)[i]), impl_root, rec);
    }
    catch(const Exception & ex){
        stringstream ss;
        ss << "[" << i << "]";
        throw ConfigurationException(format_error(*this, ss.str(), ex.message));
    }
    catch(std::exception & ex){
        stringstream ss;
        ss << "[" << i << "]";
        throw ConfigurationException(format_error(*this, ss.str(), ex.what()));
    }
    catch(...){
        stringstream ss;
        ss << "[" << i << "]";
        throw ConfigurationException(format_error(*this, ss.str(), "unknown exception"));
    }
}

Setting Setting::operator[](const std::string & name) const{
    try{
        if(rec)rec->markAsUsed(get_path());
        return Setting::resolve_link(&((*impl_this)[name]), impl_root, rec);
    }
    catch(const Exception & ex){
        throw ConfigurationException(format_error(*this, name, ex.message));
    }
    catch(std::exception & ex){
        throw ConfigurationException(format_error(*this, name, ex.what()));
    }
    catch(...){
        throw ConfigurationException(format_error(*this, name, "unknown exception"));
    }
}

bool Setting::exists(const std::string & path) const{
    return impl_this->exists(path);
}

std::string Setting::get_name() const{
    if(original_name != "") return original_name;
    return impl_this->get_name();
}

std::string Setting::get_path() const{
    return impl_this->get_path();
}

Setting::Type Setting::get_type() const{
    return impl_this->get_type();
}

SettingImplementation::~SettingImplementation(){}

/* LibconfigSetting */
LibconfigSetting::LibconfigSetting(): setting(0){}

LibconfigSetting::LibconfigSetting(const libconfig::Setting * setting_): setting(setting_){}

Setting LibconfigSetting::parse(const std::string & cfg_string, const boost::shared_ptr<SettingUsageRecorder> & rec){
    std::auto_ptr<LibconfigSetting> s(new LibconfigSetting());
    s->cfg.reset(new libconfig::Config());
    try{
        s->cfg->readString(cfg_string);
    }
    catch(libconfig::ParseException & p){
        stringstream ss;
        ss << "Error parsing configuration in line " << p.getLine() << ": " << p.getError();
        throw ConfigurationException(ss.str());
    }
    s->setting = &(s->cfg->getRoot());
    std::auto_ptr<SettingImplementation> si(s.release());
    return Setting(si, rec);
}

LibconfigSetting::operator bool() const{
    return *setting;
}

LibconfigSetting::operator std::string() const{
    return *setting;
}

LibconfigSetting::operator long int() const{
    libconfig::Setting::Type typ = setting->getType();
    switch(typ){
        case libconfig::Setting::TypeInt64: return *setting;
        case libconfig::Setting::TypeInt: return static_cast<int>(*setting);
        default:
            throw ConfigurationException("wrong type"); // will be converted to better message in Setting
    }
}

LibconfigSetting::operator double() const{
    return *setting;
}

size_t LibconfigSetting::size() const{
    return setting->getLength();
}

const SettingImplementation & LibconfigSetting::operator[](int i) const{
    libconfig::Setting * s = &((*setting)[i]);
    if(children.find(s)==children.end()){
        children.insert(s, new LibconfigSetting(s));
    }
    return children[s];
}

const SettingImplementation & LibconfigSetting::operator[](const std::string & name) const{
    libconfig::Setting * s = &((*setting)[name]);
    if(children.find(s)==children.end()){
        children.insert(s, new LibconfigSetting(s));
    }
    return children[s];
}

bool LibconfigSetting::exists(const std::string & path) const{
    return setting->exists(path);
}

std::string LibconfigSetting::get_name() const{
    const char * res = setting->getName();
    return res == 0 ? "": res;
}

std::string LibconfigSetting::get_path() const{
    return setting->getPath();
}

Setting::Type LibconfigSetting::get_type() const{
    libconfig::Setting::Type typ = setting->getType();
    switch(typ){
        case libconfig::Setting::TypeInt:
        case libconfig::Setting::TypeInt64: return Setting::TypeInt;
        case libconfig::Setting::TypeFloat: return Setting::TypeFloat;
        case libconfig::Setting::TypeString: return Setting::TypeString;
        case libconfig::Setting::TypeBoolean: return Setting::TypeBoolean;
        case libconfig::Setting::TypeArray:
        case libconfig::Setting::TypeList: return Setting::TypeList;
        case libconfig::Setting::TypeGroup:  return Setting::TypeGroup;
        case libconfig::Setting::TypeNone:
            throw logic_error("LibconfigSetting::get_type: TypeNone encountered");
            
    }
    throw logic_error("LibconfigSetting::get_type");
}

LibconfigSetting::~LibconfigSetting(){}
