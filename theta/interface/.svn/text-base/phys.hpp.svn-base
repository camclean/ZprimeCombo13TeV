#ifndef PHYS_HPP
#define PHYS_HPP

#include "interface/decls.hpp"
#include "interface/variables.hpp"
#include "interface/histogram.hpp"
#include "interface/producer.hpp"
#include "interface/histogram-with-uncertainties.hpp"

#include <vector>


namespace theta {
    
    /** \brief A real-valued function which depends on one or more parameters
     *
     * This is the base class for function plugins.
     */
    class Function{
    public:
        /// Define this as the base_type for derived classes; required for the plugin system
        typedef Function base_type;
        
        /** \brief Evaluate the function at the given parameter values.
         *
         * @return The function value at \c v.
         */
        virtual double operator()(const ParValues & v) const = 0;

        /** \brief Evaluate the function, using the parameter values given as array of doubles.
         *
         * This does the same as operator()(const ParValues&), however, it takes a pointer
         * to an array of doubles instead.
         *
         * This function is provided to make it easier to provide an interface to other libraries, which
         * do not need to know about the theta-style of variable handling (i.e.,
         * \c ParId, \c ParValue, \c VarIdManager classes).
         *
         * The translation from this "external" array format into the internal (i.e.,
         * \c ParValues) format is done by simultaneously iterating over the
         * ParIds as returned by get_parameters() and the given array \c x.
         * (As iteration over a \c ParIds object always yields the same order of parameters,
         * this mapping is well-defined.)
         */
        double operator()(const double * x) const{
            if(vpids.size()==0){
                fill_vpids();
            }
            typedef std::vector<ParId>::const_iterator it_type;
            const it_type end = vpids.end();
            size_t i = 0;
            for(it_type it=vpids.begin(); it!=end; ++it, ++i){
                theta_assert(!std::isnan(x[i]));
                pv.set(*it, x[i]);
            }
            return operator()(pv);
        }


        /** \brief Returns the parameters this function depends on
         */
        const ParIds & get_parameters() const{
            return par_ids;
        }
        
        /// Declare destructor virtual as polymorphic access to derived classes will happen.
        virtual ~Function();
        
        virtual double eval_with_derivative(const ParValues & v, ParValues & der) const;
        
    protected:
        /** \brief The parameters this function depends on
         *
         * Has to be set by derived classes in their constructor
         */
        ParIds par_ids;
        
    private:
        mutable ParValues pv; //saving this class-wide and not in operator()(const double*) saves quiet some time ...
        mutable std::vector<ParId> vpids; // will be filled in operator()(const double *)
        
        void fill_vpids() const;
    };
    
    /** \brief Function which returns the value of one specific parameter
     */
    class IdFunction: public Function{
    private:
        ParId pid;
    public:
        IdFunction(const ParId & pid_): pid(pid_){
            par_ids.insert(pid);
        }
        virtual double operator()(const ParValues & values)const{
            return values.get(pid);
        }
    };
    
    
    /** \brief A data-providing class, can be used as base class in the plugin system
     *
     * DataSource classes are used as part of a run, which, for each pseudo
     * experiment, calls the DataSource::fill function to get the pseudo data.
     */
    class DataSource: public ProductsSource{
    public:
        
        /// Define this as the base_type for derived classes; required for the plugin system
        typedef DataSource base_type;
        
        /** \brief Fill the provided Data instance with data
         *
         * The method sets the histograms of the observables it has data for and
         * does not change any other Data in \c dat.
         */
        virtual void fill(Data & dat) = 0;
        
        /// Declare destructor virtual as polymorphic access to derived classes will happen.
        virtual ~DataSource(){}
        
    protected:
        /// proxy to ProductsSource constructor for derived classes
        DataSource(const theta::Configuration & cfg): ProductsSource(cfg){}
        DataSource(const std::string & name, const boost::shared_ptr<ProductsSink> & sink): ProductsSource(name, sink){}
    };
    

}

#endif
