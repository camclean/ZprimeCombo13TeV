#include "plugins/mcmc-strategies.hpp"
#include "interface/phys.hpp"
#include "interface/asimov-utils.hpp"
#include "interface/random.hpp"
#include "interface/plugin.hpp"
#include "interface/model.hpp"
#include "interface/distribution.hpp"

using namespace theta;

//TODO: rewrite more generally by using policy pattern with templates ...
//TODO: also implement the combine-like ortho in which some parameters are always stepped, but only one of the others ...
///   this could be done by a more general framework which would make steps in different directions with different probability.
//    Current ortho would be 1/n for each of the n axes. New would be 1.0 for beta_signal and 1/n for each of the n nuisance parameters.

// oneatatime:
// - if true, the parameter value is changed only along one (randomly selected) axis in which the step size is non zero. In this case, the
//    jump kernel is a Gaussian along that direction with default std.-dev of 3 * width
// - if false, all parameters are varied simultanesouly. The jump kernel is a multivariate gaussian with a diagonal covariance
//   given by 2.38 * width / sqrt(n)
void metropolisHastings_diag(const theta::Function & nllikelihood, const MCMCOptions & opts, MCMCResult &res, theta::Random & rand, const std::vector<double> & widths,
                             bool oneatatime){
    const size_t npar = opts.startvalues.size();
    theta_assert(npar == res.getnpar());
    theta_assert(npar == widths.size());
    const ParIds & parameters = nllikelihood.get_parameters();
    theta_assert(parameters.size() == npar);

    //0. initialization:
    std::vector<double> x(npar);
    std::vector<double> x_new(npar);
    

    //set the starting point:
    std::copy(opts.startvalues.begin(), opts.startvalues.end(), &x[0]);
    double nll = nllikelihood(&x[0]);
    if(!std::isfinite(nll) and not opts.ignore_inf_nll_start) throw Exception("first nll value was inf");
    // if this exception is thrown, it means that the likelihood function at the start values was inf or NAN.
    // One common reason for this is that the model produces a zero prediction in some bin where there is >0 data which should
    // be avoided by the model (e.g., by filling in some small number in all bins with content zero).
    
    const double f = 2.38 / sqrt(npar);

    const size_t iter = opts.burn_in + opts.iterations;
    size_t weight = 1;
    for (size_t it = 1; it < iter; ++it) {
        if(oneatatime){
            // randomly choose a dimension:
            size_t index;
            do{
                index = rand.uniform() * npar;
            }while(widths[index] == 0.0);
            for (size_t i = 0; i < npar; i++) {
                x_new[i] = x[i];
            }
            x_new[index] += rand.gauss(widths[index] * opts.factor * 3.0);
        }
        else{
            for (size_t i = 0; i < npar; i++) {
                x_new[i] = x[i] + rand.gauss(widths[i] * opts.factor * f);
            }
        }
        double nll_new = nllikelihood(&x_new[0]);
        if ((nll_new <= nll) || (rand.uniform() < exp(nll - nll_new))) {
            if(it > opts.burn_in){
                res.fill(&x[0], nll, weight);
                weight = 1;
            }
            swap(x, x_new);
            nll = nll_new;
        } else if(it > opts.burn_in){
            ++weight;
        }
    }
    res.fill(&x[0], nll, weight);
}


/* asimov_der_cov */
asimov_der_cov::asimov_der_cov(const theta::Configuration & cfg): MCMCStrategy(cfg), epsilon(1e-4){
    Setting s = cfg.setting;
    if(s.exists("epsilon")){
        epsilon = s["epsilon"];
    }
}

void asimov_der_cov::init(const theta::Model & model, const boost::shared_ptr<theta::Distribution> & override_parameter_distribution){
    sqrt_cov = asimov_likelihood_matrix(model, override_parameter_distribution, epsilon);
    sqrt_cov.cholesky_decomposition();
    const Distribution & dist = override_parameter_distribution.get()? *override_parameter_distribution : model.get_parameter_distribution();
    // fill start:
    ParValues pv_start;
    dist.mode(pv_start);
    options.startvalues.resize(dist.get_parameters().size());
    pv_start.fill(&options.startvalues[0], dist.get_parameters());
}

void asimov_der_cov::run_mcmc(const theta::Function & nllikelihood, MCMCResult &res) const{
    metropolis_hastings_multigauss(nllikelihood, options, res, *rnd_gen, sqrt_cov);
}


/* asimov_widths_1d */
asimov_widths_1d::asimov_widths_1d(const theta::Configuration & cfg): MCMCStrategy(cfg){
}

void asimov_widths_1d::init(const theta::Model & model, const boost::shared_ptr<theta::Distribution> & override_parameter_distribution){
    ParValues pv_widths = asimov_likelihood_widths(model, override_parameter_distribution);
    const Distribution & dist = override_parameter_distribution.get()? *override_parameter_distribution : model.get_parameter_distribution();
    widths.resize(dist.get_parameters().size());
    pv_widths.fill(&widths[0], dist.get_parameters());
    ParValues pv_start;
    dist.mode(pv_start);
    options.startvalues.resize(dist.get_parameters().size());
    pv_start.fill(&options.startvalues[0], dist.get_parameters());
}

void asimov_widths_1d::run_mcmc(const theta::Function & nllikelihood, MCMCResult &res) const{
    metropolisHastings_diag(nllikelihood, options, res, *rnd_gen, widths, true);
}

/* asimov_widths */
asimov_widths::asimov_widths(const theta::Configuration & cfg): MCMCStrategy(cfg){
}

void asimov_widths::init(const theta::Model & model, const boost::shared_ptr<theta::Distribution> & override_parameter_distribution){
    ParValues pv_widths = asimov_likelihood_widths(model, override_parameter_distribution);
    const Distribution & dist = override_parameter_distribution.get()? *override_parameter_distribution : model.get_parameter_distribution();
    widths.resize(dist.get_parameters().size());
    pv_widths.fill(&widths[0], dist.get_parameters());
    ParValues pv_start;
    dist.mode(pv_start);
    options.startvalues.resize(dist.get_parameters().size());
    pv_start.fill(&options.startvalues[0], dist.get_parameters());
}

void asimov_widths::run_mcmc(const theta::Function & nllikelihood, theta::MCMCResult &res) const{
    metropolisHastings_diag(nllikelihood, options, res, *rnd_gen, widths, false);
}

REGISTER_PLUGIN(asimov_widths_1d);
REGISTER_PLUGIN(asimov_der_cov);
REGISTER_PLUGIN(asimov_widths);
