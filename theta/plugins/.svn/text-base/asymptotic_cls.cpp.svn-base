#include "plugins/asymptotic_cls.hpp"
#include "interface/minimizer.hpp"
#include "interface/secant.hpp"
#include "interface/model.hpp"
#include "interface/phys.hpp"
#include "interface/asimov-utils.hpp"
#include "interface/plugin.hpp"
#include "interface/database.hpp"
#include "interface/redirect_stdio.hpp"

#include <cmath>

#include <boost/math/distributions/normal.hpp>

using namespace theta;
using namespace std;

namespace{

const boost::math::normal norm;

double phi(double x){
    return boost::math::cdf(norm, x);
}

// Calculate sigma from the likelihood ratio from asimov data (this is the method called "sigma_A" in the paper).
// Also calculate the likelihood ratio for data.
// These calculations are encapsulated into a class here to do the setup of parameter start values, ranges, etc. for the minimizer only
// once.
class SigmaCalculator{
public:
    SigmaCalculator(const Model & m_, Minimizer & min_, const ParId & pid_, const boost::shared_ptr<NLLikelihood> & data_nll_):
          model(m_), minimizer(min_), pid(pid_), data_nll(data_nll_), data_init(false){
        model.get_parameter_distribution().mode(start);
        const ParIds & pids = model.get_parameter_distribution().get_parameters();
        for(ParIds::const_iterator pit=pids.begin(); pit!=pids.end(); ++pit){
            ranges[*pit] = model.get_parameter_distribution().support(*pit);
        }
        widths = asimov_likelihood_widths(model, boost::shared_ptr<Distribution>());
    }

    double get_poi_width() const{
        return widths.get(pid);
    }
    
    void set_data_nll(const boost::shared_ptr<NLLikelihood> & data_nll_){
        data_nll = data_nll_;
        data_init = false;
    }

    // beta_toy is the value of mu assumed for toys (if there were toys). This is called mu' in the paper.
    // beta_q is the value of mu in the test statistic calculation of q_mu-tilde
    double sigma(double beta_toy, double beta_q) const {
        ParValues start_(start);
        start_.set(pid, beta_toy);
        asimov_data_nll nll(model, boost::shared_ptr<Distribution>(), start_);
        const double nll0 = nll(start_);
        // restict poi to beta_q and fit:
        start_.set(pid, beta_q);
        map<ParId, pair<double, double> > ranges_(ranges);
        ranges_[pid] = make_pair(beta_q, beta_q);
        ParValues widths_(widths);
        widths_.set(pid, 0.0);
        MinimizationResult minres = minimizer.minimize(nll, start_, widths_, ranges_);
        const double nll1 = minres.fval;

        // found minimum should be worse than the parameter value used to construct the asimov dataset with:
        if(nll1 <= nll0){
            theta::err << "Warning: asimov nll weird: found better value than at true input of asimov data. Asymptotic will be unreliable." << endl;
            return numeric_limits<double>::quiet_NaN();
        }
        return fabs(beta_toy - beta_q) / sqrt(2 * (nll1 - nll0));
    }

    double data_tsvalue(double beta_q) const{
        theta_assert(data_nll.get());
        if(!data_init){
            // determine global minimum:
            MinimizationResult minres = minimizer.minimize(*data_nll, start, widths, ranges);
            nll_data_min = minres.fval;
            beta_min = minres.values.get(pid);
            data_init = true;
        }
        if(beta_q < beta_min){
            return 0.0;
        }
        // restrict poi to beta_q and fit:
        ParValues start_(start);
        start_.set(pid, beta_q);
        map<ParId, pair<double, double> > ranges_(ranges);
        ranges_[pid] = make_pair(beta_q, beta_q);
        ParValues widths_(widths);
        widths_.set(pid, 0.0);
        MinimizationResult minres = minimizer.minimize(*data_nll, start_, widths_, ranges_);
        const double nll1 = minres.fval;
        return 2 * (nll1 - nll_data_min);
    }
private:
    const Model & model;
    Minimizer & minimizer;
    ParId pid;

    // minimizer start / step settings:
    ParValues start, widths;
    map<ParId, pair<double, double> > ranges;

    // for data ts value calculation:
    boost::shared_ptr<NLLikelihood> data_nll; // can be 0, if only expected limits are calculated
    // cache values for global minimum:
    mutable double beta_min, nll_data_min;
    mutable bool data_init;
};

namespace asymptotic_ts_dist{
    // quantile function; inverse of cumulative distribution function
    double qf(double beta_toy, double beta_q, double sigma, double q){
        theta_assert(q > 0.0 and q < 1.0);
        double F = phi(beta_toy / sigma);
        if(q < F){
            double res = theta::utils::phi_inverse(q) + (beta_q - beta_toy) / sigma;
            return pow(max(res, 0.0), 2);
        }
        else{
            return 2 * utils::phi_inverse(q) * beta_q / sigma + (beta_q*beta_q - 2 * beta_q * beta_toy) / (sigma*sigma);
        }
    }

    // survival function = 1 - cumulative distribution function
    // tsval is 2 * LR
    //note: can be done numerically more stable for small values of sf.
    // But on the other hand, we usually don't have clb / clsb values THAT small.
    double sf(double beta_toy, double beta_q, double sigma, double tsval){
        double cdf;
        if(tsval <= pow(beta_q / sigma,2)){
            cdf = phi(sqrt(tsval) - (beta_q - beta_toy) / sigma);
        }
        else{
            cdf = phi((tsval - (beta_q*beta_q - 2*beta_q * beta_toy) / pow(sigma, 2)) / (2 * beta_q / sigma));
        }
        return 1.0 - cdf;
    }

};

// asymptotic (expected cls value) - (target cls value) as a function of beta_q.
// Suitable to pass to a root finding routine that expected a function double -> double.
class cls_expected{
public:
    cls_expected(double beta_expected_, double target_cls_, double q_, const SigmaCalculator & sc_):
        beta_expected(beta_expected_), target_cls(target_cls_), q(q_), sc(sc_){}

    double operator()(double beta_q)const{
        const double w = sc.get_poi_width();
        if(beta_q < beta_expected + 1e-2 * w) return 1.0;
        double sigma_expected = sc.sigma(beta_expected, beta_q);
        if(std::isnan(sigma_expected)){
            return 1.0;
        }
        double ts_expected = asymptotic_ts_dist::qf(beta_expected, beta_q, sigma_expected, q);
        double sigma0 = beta_expected == 0.0 ? sigma_expected : sc.sigma(0.0, beta_q);
        double clb = asymptotic_ts_dist::sf(0.0, beta_q, sigma0, ts_expected); // TODO: should be 1 - q. Can check and / or save this call.
        double clsb = asymptotic_ts_dist::sf(beta_q, beta_q, sigma0, ts_expected);
        return clsb / clb - target_cls;
    }

private:
    double beta_expected, target_cls, q;
    const SigmaCalculator & sc;
};

// same as cls_expected, but for observed CLs value.
class cls_observed{
public:
    cls_observed(const SigmaCalculator & sc_, double target_cls_): sc(sc_), target_cls(target_cls_){
    }

    double operator()(double beta_q)const{
        double ts_data = sc.data_tsvalue(beta_q);
        // the test statistic as defined here is expected to be >= 0.0 always by definition.
        // However, because of numerical inaccuracies, it can be < 0.0, especially if beta_q is close to
        // the minimum found for data. In this case, cls is 1.0. The cutoff of 1e-2
        // is reached at a value of beta_q only a little above the best fit value, which will be well
        // below the limit. Note that this assumes confidence levels > 50%.
        if(ts_data < 1e-2) return 1.0 - target_cls;
        double sigma0 = sc.sigma(0.0, beta_q);
        double clb = asymptotic_ts_dist::sf(0.0, beta_q, sigma0, ts_data);
        if(clb <= 0.0) return 1.0; // actually infinity, but it does not matter as long as it is > target_cls ...
        double clsb = asymptotic_ts_dist::sf(beta_q, beta_q, sigma0, ts_data);
        return clsb / clb - target_cls;
    }

private:
    const SigmaCalculator & sc;
    double target_cls;
};


// convenience function for getting the ParId from the name in cfg_name in the current Setting cfg.setting
ParId getpar(const Configuration & cfg, const std::string & cfg_name){
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    return vm->get_par_id(cfg.setting[cfg_name]);
}


// find the CLs limit by finding the root of the function f on the interval [0, inf)
// width should be the "typical scale" of the parameter.
// xtol_rel is the requested precision of the limit, relative to the width.
template<typename FT>
double find_limit(double width, const FT & f, double xtol_rel){
    double beta_low = 0.0;
    double beta_high = width;
    double cls_low = f(beta_low);
    theta_assert(cls_low > 0.0);
    double cls_high = f(beta_high);
    theta_assert(isfinite(cls_high));
    int n_iterations = 0;
    double scale = 1.0;
    while(cls_high >= 0.0){
        if(cls_high > 0){
            beta_low = beta_high;
            cls_low = cls_high;
        }
        beta_high += width * scale;
        scale *= 1.05;
        cls_high = f(beta_high);
        theta_assert(isfinite(cls_high));
        ++n_iterations;
        if(n_iterations > 100){
            stringstream ss;
            ss << "Could not find the CLs limit: CLs value depends too little on parameter or initial width was too small; width=" << width << "; latest tested interval: " << beta_low << " -- " << beta_high << " (f = cls - target_cls: " << cls_low << " -- " << cls_high << ")";
            throw invalid_argument(ss.str());
        }
    }
    return brent(f, beta_low, beta_high, xtol_rel * width, cls_low, cls_high, 1e-5);
}

}


void asymptotic_cls::run(){
    int total = quantiles_expected.size() + n;
    int done = 0;
    if(progress_listener){
        progress_listener->progress(done, total, 0);
    }
    // get mode, ranges and width:
    ParValues mode;
    model->get_parameter_distribution().mode(mode);
    const ParIds & pids = model->get_parameter_distribution().get_parameters();
    map<ParId, pair<double, double> > ranges;
    for(ParIds::const_iterator pit=pids.begin(); pit!=pids.end(); ++pit){
        ranges[*pit] = model->get_parameter_distribution().support(*pit);
    }
    ParValues widths = asimov_likelihood_widths(*model, boost::shared_ptr<Distribution>());
    boost::shared_ptr<NLLikelihood> data_nll;
    SigmaCalculator calc(*model, *minimizer, parameter, data_nll);

    // the expected limits:
    Column c_index = limits_table->add_column("index", typeInt);
    Column c_q = limits_table->add_column("q", typeDouble);
    Column c_limit = limits_table->add_column("limit", typeDouble);
    for(size_t i=0; i<quantiles_expected.size(); ++i){
        cls_expected ex(parameter_value_expected, 1 - cl, 1 - quantiles_expected[i], calc);
        double limit = find_limit(calc.get_poi_width(), ex, limit_reltol);
        Row r;
        r.set_column(c_index, static_cast<int>(i));
        r.set_column(c_q, quantiles_expected[i]);
        r.set_column(c_limit, limit);
        limits_table->add_row(r);
        if(progress_listener){
            progress_listener->progress(++done, total, 0);
        }
    }
    // the "observed" limit(s):
    int errors = 0;
    Data observed_data;
    for(int i=0; i<n; ++i){
        data->fill(observed_data);
        data_nll = model->get_nllikelihood(observed_data);
        calc.set_data_nll(data_nll);
        cls_observed obs(calc, 1 - cl);
        try{
            double limit = find_limit(calc.get_poi_width(), obs, limit_reltol);
            Row r;
            r.set_column(c_q, 0.0);
            r.set_column(c_index, static_cast<int>(quantiles_expected.size() + i));
            r.set_column(c_limit, limit);
            limits_table->add_row(r);
        }
        catch(MinimizationException & ex){
            ++errors;
        }
        if(progress_listener){
            progress_listener->progress(++done, total, errors);
        }
    }
}

asymptotic_cls::asymptotic_cls(const Configuration & cfg): 
  parameter(getpar(cfg, "parameter")), cl(0.95), n(1), parameter_value_expected(0.0), limit_reltol(1e-3){
    // add some config for submodules:
    cfg.pm->set("runid", boost::shared_ptr<int>(new int(1)));
    output_database = PluginManager<Database>::build(Configuration(cfg, cfg.setting["output_database"]));
    std::auto_ptr<Table> rnd_table = output_database->create_table("rndinfo");
    boost::shared_ptr<RndInfoTable> rnd_info(new RndInfoTable(rnd_table));
    cfg.pm->set("default", rnd_info);
    cfg.pm->set("default", boost::shared_ptr<ProductsSink>(new NullProductsSink()));
    limits_table = output_database->create_table("limits");
    
    Setting s = cfg.setting;
    if(s.exists("cl")) cl = s["cl"];
    // we assume cl > 0.5 at some places, but sometimes, we need a small margin, so test for 0.6 here:
    if(cl < 0.6){
        throw ConfigurationException("cl < 0.6 not supported.");
    }
    model = PluginManager<Model>::build(Configuration(cfg, cfg.setting["model"]));
    minimizer = PluginManager<Minimizer>::build(Configuration(cfg, cfg.setting["minimizer"]));
    if(s.exists("data")){
        data = PluginManager<DataSource>::build(Configuration(cfg, cfg.setting["data"]));
        if(s.exists("n")) n = s["n"];
    }
    else{
        n = 0;
    }
    if(s.exists("limit_reltol")) limit_reltol = s["limit_reltol"];
    if(limit_reltol <= 1e2 * numeric_limits<double>::epsilon()){
        throw ConfigurationException("limit_reltol is too small");
    }
    if(s.exists("quantiles_expected")){
        Setting se = s["quantiles_expected"];
        const unsigned int n = se.size();
        quantiles_expected.reserve(n);
        for(unsigned int i=0; i<n; ++i){
            quantiles_expected.push_back(se[i]);
        }
    }
    else{
        quantiles_expected.reserve(5);
        quantiles_expected.push_back(0.025);
        quantiles_expected.push_back(0.16);
        quantiles_expected.push_back(0.5);
        quantiles_expected.push_back(0.84);
        quantiles_expected.push_back(0.975);
    }
    if(s.exists("parameter_value_expected")) parameter_value_expected = s["parameter_value_expected"];
}

REGISTER_PLUGIN(asymptotic_cls)

