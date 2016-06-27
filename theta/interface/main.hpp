#ifndef MAIN_HPP
#define MAIN_HPP

#include <boost/shared_ptr.hpp>

namespace theta{

/** \brief A global flag indicating whether to stop execution as soon as possible
 *
 * All Run classes should check this flag regulary. If set to true, the Main object
 * should return as soon as possible without violating the consistency of the result, if possible.
 *
 * This variable is set to true once theta receives SIGINT.
 */
extern volatile bool stop_execution;

/** \brief Install the SIGINT signal handler which sets stop_execution
 *
 * Calling this function will use the sigaction function to install a signal handler for
 * SIGINT which sets stop_execution to true to indicate to the Main::run
 * to terminate as soon as possible. If stop_execution is already set to true,
 * exit(1) will be called immediately by the handler.
 */
void install_sigint_handler();

/** \brief A callback class called by \link Main \endlink to indicate progress
 *
 * User interfaces can derive their classes from  ProgressListener and implement
 * some sort of feedback for the user ("progress bar"). An instance
 * of this user-defined class should be registered with the
 * Main::set_progress_listener() method.
 */
class ProgressListener{
public:
   /** \brief Indicate progress.
    *
    * This function will be called by Main to indicate the current progress.
    * \param done how many units of work have been done
    * \param total total units of work
    * \param how many failures there have been in 'done'
    *
    * In some cases, the 'total units of work' are not known. In this case, it is still useful
    * to indicate some activity / progress and total should be set to a value <= 0 in that case.
    */
    virtual void progress(int done, int total, int errors) = 0;
    
    /// Make destructor virtual as we have polymorphic access
    virtual ~ProgressListener();
};


/** \brief Class Representing the "main" setting in the configuration file
 *
 *
 */
class Main{
public:
    /// Define this class the base_type for derived classes; required for the plugin system
    typedef Main base_type;
    
   /** \brief Register progress listener.
    *
    * The registered progress listener will be informed about the current progress of the run during
    * execution of Run::run().
    */
    void set_progress_listener(const boost::shared_ptr<ProgressListener> & l);
    
    /** \brief Perform the actual run
     *
     * The meaning depends on derived classes.
     */
     virtual void run() = 0;
     
     
     /// Declare destructor virtual as we will have polymorphic access to derived classes
     virtual ~Main();
    
protected:
    boost::shared_ptr<ProgressListener> progress_listener;

};

}

#endif

