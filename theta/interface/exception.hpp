#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <string>
#include <stdexcept>

void fail_assert(const char * filename, int lineno, const char * expression);

// define our own assert macro because the standard assert does an abort which does not destroy object properly
// theta_assert raises a FatalException, which -- as part of exception handling -- destroys all objects.
#define theta_assert(expression) if(!(expression)) { ::fail_assert(__FILE__, __LINE__, #expression); }

namespace theta {

std::string demangle(const std::string & typename_);

/** \brief Base class for the runtime exceptions used in %theta
 */
class Exception: public std::runtime_error {
public:
    /// The human-readable message for the user
    std::string message;
    
    /// Constructor taking a message intended for the user which will be written to Exception::message
    Exception(const std::string & msg);
    
    /// Declare destructor virtual
    virtual ~Exception() throw(){}
    
    /// override std::exception::what to print out the message
    virtual const char* what() const throw();
    
protected:
    /// buffer for the const char * returned by what()
    mutable std::string whatstring;
};

/** \brief Thrown during configuration file processing and initial object construction
 */
class ConfigurationException: public Exception{
public:
    /// Constructor taking a message intended for the user which will be written to Exception::message
    ConfigurationException(const std::string & msg);
};

/// \brief Thrown in case of database errors.
class DatabaseException: public Exception{
public:
    /// Constructor taking a message intended for the user which will be written to Exception::message
    DatabaseException(const std::string & s);
};

/** \brief Thrown in case of minimization errors.
 *
 * \sa theta::Minimizer
 */
class MinimizationException: public Exception{
public:
    /// Constructor taking a message intended for the user which will be written to Exception::message
    MinimizationException(const std::string & s);
};


/** \brief Thrown in case an immediate exit of the main program is requested. 
 *
 * As this exception is caught at a different place, it is not part of the usual Exception
 * hierarchy.
 */
class ExitException {
public:
   std::string message;

   /// Constructor taking a message intended for the user which will be written to Exception::message
   ExitException(const std::string & message_);
};

}

#endif

