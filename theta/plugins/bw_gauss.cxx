#include "interface/histogram-function.hpp"
#include "interface/plugin.hpp"

#include <complex.h>
#include <fftw3.h>
#include <complex>

// NOTE: this source file is not compiled per default, as only .cpp files are compiled, not .cxx
// To compile it:
//  - make sure you have fftw library available
//  - rename to .cpp
//  - change Makefile in this directory to link against fftw3

using namespace std;
using namespace theta;

namespace {

// helper class dor bw_gauss, see below
class GaussConvolution {
public:
    
    // spectrum will be used as both input and output. spectrum should outlive this object
    // (at least all calls to convolve!).
    explicit GaussConvolution(size_t nbins, double * spectrum, double binwidth_): gauss(nbins), binwidth(binwidth_){
        gauss_trafo = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nbins/2 + 1));
        s_trafo = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nbins/2 + 1));
        
        p_gauss = fftw_plan_dft_r2c_1d(nbins, &gauss[0], gauss_trafo, FFTW_ESTIMATE);
        p_spectrum = fftw_plan_dft_r2c_1d(nbins, spectrum, s_trafo, FFTW_ESTIMATE);
        p_spectrum_back = fftw_plan_dft_c2r_1d(nbins, s_trafo, spectrum, FFTW_ESTIMATE);
    }
    
    
    void convolve(double gauss_width){
        const size_t nbins = gauss.size();
        for(size_t i=0; i<nbins/2; ++i){
            double x = binwidth * (i + 0.5);
            gauss[i] = exp(-x*x / (2 * gauss_width * gauss_width));
        }
        for(size_t i=nbins/2; i<nbins; ++i){
            double x = -binwidth * ((nbins - i) - 0.5);
            gauss[i] = exp(-x*x / (2 * gauss_width * gauss_width));
        }
        // transform both spectrum and gauss:
        fftw_execute(p_gauss);
        fftw_execute(p_spectrum);
        const size_t nbins_trafo = nbins/2+1;
        // multiply pointwise:
        for(size_t i=0; i<nbins_trafo; ++i){
            complex<double> gauss = reinterpret_cast<complex<double>*>(gauss_trafo)[i];
            //cout << i << " " << std::abs(gauss) << " " << arg(gauss) << endl;
            reinterpret_cast<complex<double>*>(s_trafo)[i] *= gauss;
        }
        // transform back:
        fftw_execute(p_spectrum_back);
    }
    
    ~GaussConvolution(){
        fftw_destroy_plan(p_spectrum);
        fftw_destroy_plan(p_spectrum_back);
        fftw_destroy_plan(p_gauss);
        fftw_free(s_trafo);
        fftw_free(gauss_trafo);
    }
    
private:
    vector<double> gauss;
    double binwidth;
    
    fftw_complex * gauss_trafo;
    fftw_complex * s_trafo;
    
    fftw_plan p_gauss, p_spectrum_back, p_spectrum;
};

}

/** \brief A HistogramFunction which is the convolution of a Breit-Wigner with a Gaussian resolution function
 *
 * This can be useful for modeling a Z->ll mass peak.
 *
 * \code
 *  {
 *     type = "landau_c_gauss";
 *     bw_mean = "mean";
 *     bw_width = "zmass";
 *     g_width = "eeres";
 *     nbins = 60;
 *     xrange = [10.0, 20.0];
 *     oversample = 5;
 *  }
 * \endcode
 *
 *
 * Oversample determines the internal accuracy: the convolution is done with nbins * oversample bins;
 * for each bin in the output oversample bins are summed.
 */
class bw_gauss : public theta::HistogramFunction {
public:
    
    /** \brief Constructor used by the plugin system to build an instance from settings in a configuration file
     */
    bw_gauss(const theta::Configuration & cfg);
    
    virtual void add_with_coeff_to(theta::Histogram1DWithUncertainties & hres, double coeff, const theta::ParValues & values) const;
    virtual void add_with_coeff_to(theta::Histogram1D & hres, double coeff, const theta::ParValues & values) const;
    virtual void get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const;

    virtual ~bw_gauss(){}

private:
    void fill_h(const theta::ParValues & values) const;
    
    vector<theta::ParId> parameters;

    mutable theta::Histogram1D h_oversampled;
    mutable theta::Histogram1D h;
    std::auto_ptr<GaussConvolution> gc;
};


void bw_gauss::fill_h(const ParValues & values) const {
    // fill Breit-Wigner in oversampled:
    const double binwidth_o = (h_oversampled.get_xmax() - h_oversampled.get_xmin()) / h_oversampled.get_nbins();
    const double bw_mean = values.get(parameters[0]);
    const double bw_width = values.get(parameters[1]);
    for(size_t i=0; i<h_oversampled.get_nbins(); ++i){
        const double x = h_oversampled.get_xmin() + (i + 0.5) * binwidth_o;
        h_oversampled.set(i, 1.0 / (pow(x*x - bw_mean*bw_mean, 2) + bw_mean*bw_mean*bw_width*bw_width));
    }
    // make convolustion with Gaussian:
    const double g_width = values.get(parameters[2]);
    gc->convolve(g_width);
    // fill h_oversampled into h, and normalize to 1.0:
    double sum = h_oversampled.get_sum();
    size_t oversample = h_oversampled.get_nbins() / h.get_nbins();
    for(size_t i=0; i<h.get_nbins(); ++i){
        double value = h_oversampled.get(i*oversample);
        for(size_t k=1; k<oversample; ++k){
            value += h_oversampled.get(i*oversample + k);
        }
        h.set(i, value / sum);
    }
}

void bw_gauss::add_with_coeff_to(theta::Histogram1DWithUncertainties & hres, double coeff, const theta::ParValues & values) const{
    fill_h(values);
    hres.add_with_coeff(coeff, h);
}

void bw_gauss::add_with_coeff_to(theta::Histogram1D & hres, double coeff, const theta::ParValues & values) const{
    fill_h(values);
    hres.add_with_coeff(coeff, h);
}


void bw_gauss::get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const{
    nbins = h.get_nbins();
    xmin = h.get_xmin();
    xmax = h.get_xmax();
}


bw_gauss::bw_gauss(const Configuration & cfg){
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    parameters.reserve(3);
    parameters.push_back(vm->get_par_id(cfg.setting["bw_mean"]));
    parameters.push_back(vm->get_par_id(cfg.setting["bw_width"]));
    parameters.push_back(vm->get_par_id(cfg.setting["g_width"]));
    for(size_t i=0; i<3; ++i){
        par_ids.insert(parameters[i]);
    }

    int nbins = cfg.setting["nbins"];
    double xmin = cfg.setting["xmin"];
    double xmax = cfg.setting["xmax"];
    
    int oversample = 1;
    if(cfg.setting.exists("oversample")){
        oversample = cfg.setting["oversample"];
        if(oversample <= 0){
            throw ConfigurationException("oversample <= 0 not allowed");
        }
    }
    h.reset(nbins, xmin, xmax);
    h_oversampled.reset(nbins * oversample, xmin, xmax);
    
    gc.reset(new GaussConvolution(nbins * oversample, h_oversampled.get_data(), (xmax - xmin) / h_oversampled.get_nbins()));
}

REGISTER_PLUGIN(bw_gauss)

