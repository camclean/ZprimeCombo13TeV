#include "interface/exception.hpp"

#include <sstream>
#include <typeinfo>
#include <stdexcept>
#include <cstdlib>

#if __GNUC__
#include <cxxabi.h>
#endif

using namespace theta;

void fail_assert(const char * filename, int lineno, const char * expression){
    std::stringstream ss;
    ss << "Assertion '" << expression << "' failed in " << filename << ":" << lineno;
    throw std::logic_error(ss.str());
}

std::string theta::demangle(const std::string & s){
    std::string result(s);
    int status = 1;
    char * realname;
#if __GNUC__
    realname = abi::__cxa_demangle(s.c_str(), 0, 0, &status);
#endif
    if(status==0){
        result = realname;
        std::free(realname);
    }
    return result;
}

const char* Exception::what() const throw(){
    std::stringstream ss;
    ss << demangle(typeid(*this).name()) << ": " << message;
    whatstring = ss.str();
    return whatstring.c_str();
}

Exception::Exception(const std::string & m): runtime_error(m), message(m){}

ConfigurationException::ConfigurationException(const std::string & msg): Exception(msg){}

DatabaseException::DatabaseException(const std::string & s): Exception(s){}

MinimizationException::MinimizationException(const std::string & s): Exception(s){}

ExitException::ExitException(const std::string & message_): message(message_){}

