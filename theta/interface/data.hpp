#ifndef DATA_HPP
#define DATA_HPP

#include "interface/decls.hpp"
#include "interface/variables.hpp"
#include "interface/histogram.hpp"
#include "interface/histogram-with-uncertainties.hpp"
#include <vector>

namespace theta{

/** \brief Contains data more histogram observables and real-values observables
 * 
 * HT is the Histogram type and must be assignable, copy-contructible, and provide a "reset" method.
 */
template<typename HT>
class DataT {
    void fail_get(const ObsId & oid) const{
        throw std::invalid_argument("Data::operator[]() const: no data found for given ObsId");
    }
public:
    /** \brief Returns all observable ids for which there is data
     * 
     * Returns those observable ids for which there was previously
     * a Histogram1D saved via the operator[].
     */
    ObsIds get_observables() const{
        ObsIds result;
        typename std::vector<HT>::const_iterator it = data.begin();
        size_t i=0;
        for(;it!=data.end(); ++it, ++i){
            if(it->get_nbins()!=0) result.insert(ObsId(i));
        }
        return result;
    }
    
    /** \brief Access the histogram with an observable
     * 
     * The const version is usually only used to read a previously set
     * Histogram. If no Histogram is saved for the supplied observable id,
     * a invalid_argument exception will be thrown from the const version.
     */
    ///@{
    HT & operator[](const ObsId & id){
        if(id.id >= data.size()) data.resize(id.id + 1);
        return data[id.id];
    }
    
    const HT & operator[](const ObsId & id) const{
        if(id.id >= data.size() || data[id.id].get_nbins()==0) fail_get(id);
        return data[id.id];
    }
    ///@}
    
    const theta::ParValues & get_rvobs_values() const{
        return rvobs_values;
    }
    
    void set_rvobs_values(const theta::ParValues & values){
        rvobs_values = values;
    }
    
    /// \brief reset all current Histograms, i.e., set to zero entry
    void reset(){
        rvobs_values.clear();
        std::vector<Histogram1D>::iterator it = data.begin();
        for(; it!=data.end(); ++it){
            it->reset();
        }
    }
        
private:
    std::vector<HT> data;
    theta::ParValues rvobs_values;
};

}

#endif
