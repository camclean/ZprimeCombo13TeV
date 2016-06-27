#include "plugins/cls_limits.hpp"

#include "interface/cfg-utils.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"
#include "interface/variables-utils.hpp"
#include "interface/minimizer.hpp"
#include "interface/random.hpp"
#include "interface/database.hpp"
#include "interface/main.hpp"
#include "interface/phys.hpp"
#include "interface/model.hpp"
#include "interface/distribution.hpp"
#include "interface/random-utils.hpp"
#include "interface/asimov-utils.hpp"

#include <fstream>
#include <iomanip>
#include <boost/optional.hpp>

using namespace std;
using namespace theta;

// runid conventions for the products tables:
// * runid = 0 is used to calculate the ts values for the current dataset; the eventid is set to idata
// * runid >= 1 is used to make the toys at different truth values for the CLs construction


// the negative log likelihood to fit the function
// y = target_cls * exp(lambda * (x - limit)), with fixed target_cls and parameters lambda, limit
// to a dataset with correlated uncertainties for y.
class exp_nll: public Function{
    ParId lambda, limit;
    double target_cls;
    vector<double> x_values, y_values;
    Matrix inverse_covariance;
public:
    exp_nll(const ParId & pid_lambda, const ParId & pid_limit, double target_cls_, const vector<double> & x_values_, const vector<double> & y_values_, const Matrix & inverse_covariance_):
      lambda(pid_lambda), limit(pid_limit), target_cls(target_cls_), x_values(x_values_), y_values(y_values_), inverse_covariance(inverse_covariance_){
        par_ids.insert(limit);
        par_ids.insert(lambda);
        theta_assert(inverse_covariance.get_n_rows() == inverse_covariance.get_n_cols());
        theta_assert(x_values.size()==inverse_covariance.get_n_rows() && x_values.size()==y_values.size());
    }
    
    double operator()(const ParValues & values) const{
        double limit_v = values.get(limit);
        double lambda_v = values.get(lambda);
        double result = 0.0;
        const size_t N = x_values.size();
        for(size_t i=0; i < N; ++i){
            for(size_t j=0; j < N; ++j){
                result += (y_values[i] - target_cls * exp(lambda_v * (x_values[i] - limit_v))) * (y_values[j] - target_cls * exp(lambda_v * (x_values[j] -limit_v))) * inverse_covariance(i,j);
            }
        }
        return 0.5 * result;
    }
};

// forwards everything to multiple other sinks
class MultiplexingProductsSink: public ProductsSink{
private:
    std::vector<boost::shared_ptr<ProductsSink> > sinks;
    std::map<Column, std::vector<Column> > sink_columns; // for each Column we return, this holds the columns of the sinks
    int next_column_id;
public:
    MultiplexingProductsSink(const std::vector<boost::shared_ptr<ProductsSink> > & sinks_): sinks(sinks_), next_column_id(0){}

    virtual Column declare_column_impl(const std::string & full_product_name, const data_type & type){
        Column result(next_column_id++);
        for(size_t i=0; i<sinks.size(); ++i){
            sink_columns[result].push_back(sinks[i]->declare_column(full_product_name, type));
        }
        return result;
    }
    
    virtual void set_product(const Column & c, double d){
        for(size_t i=0; i<sinks.size(); ++i){
            sinks[i]->set_product(sink_columns[c][i], d);
        }
    }
    virtual void set_product(const Column & c, int k){
        for(size_t i=0; i<sinks.size(); ++i){
            sinks[i]->set_product(sink_columns[c][i], k);
        }
    }
    virtual void set_product(const Column & c, const std::string & s){
        for(size_t i=0; i<sinks.size(); ++i){
            sinks[i]->set_product(sink_columns[c][i], s);
        }
    }
    virtual void set_product(const Column & c, const Histogram1D & h){
        for(size_t i=0; i<sinks.size(); ++i){
            sinks[i]->set_product(sink_columns[c][i], h);
        }
    }
};


class BlackholeProductsSink: public ProductsSink{
public:
    virtual Column declare_column_impl(const std::string & product_name, const data_type & type){
        return Column(0);
    }
    virtual void set_product(const Column & c, double d){ }
    virtual void set_product(const Column & c, int k){ }
    virtual void set_product(const Column & c, const std::string & s){ }
    virtual void set_product(const Column & c, const Histogram1D & h){ }
};

class SaveDoubleColumn: public ProductsSink{
private:
    string column_name;
    double value;
    bool value_set;
public:
    virtual Column declare_column_impl(const std::string & product_name, const data_type & type){
        return Column(column_name == product_name ? 1 : 0);
    }
    
    virtual void set_product(const Column & c, double d){
        if(c.get_id()==1){
            value_set = true;
            value = d;
        }
    }
    virtual void set_product(const Column & c, int i){
    }
    virtual void set_product(const Column & c, const std::string & s){
    }
    virtual void set_product(const Column & c, const Histogram1D & h){
    }
    
    double consum_value() {
        if(!value_set) throw invalid_argument("column name '" + column_name + "' not set");
        value_set = false;
        return value;
    }
    
    explicit SaveDoubleColumn(const string & name): column_name(name), value(NAN), value_set(false){}
};

class data_filler: public theta::RandomConsumer, public theta::DataSource{
 public:
     enum t_toy_mode {
        ttm_prior, ttm_datafit
     };
     
     data_filler(const theta::Configuration & cfg, const boost::shared_ptr<Model> & model_, const ParId & truth_parameter, const t_toy_mode & mode_, const Data & real_data_, const boost::shared_ptr<Minimizer> & minimizer_, std::auto_ptr<ostream> & debug_out_);
     virtual void fill(Data & dat);
     void set_truth_value(double truth_value_){
         truth_value = truth_value_;
     }
     
private:
    boost::shared_ptr<Model> model;
    Column truth_column;
    ParId truth_parameter;
    double truth_value;
    
    // for debugging:
    boost::shared_ptr<VarIdManager> vm;
    std::auto_ptr<ostream> & debug_out;
    
    t_toy_mode mode;
    // all from here is only relevant for mode = ttm_datafit
    map<double, ParValues> truth_to_nuisancevalues;
    Data real_data;
    boost::shared_ptr<Minimizer> minimizer;
    
    // for the maximum likelihood fit:
    bool start_step_ranges_init;
    theta::ParValues start, step;
    theta::Ranges ranges;
};

// contains the (numerically derived) cls values and uncertainties as function of the truth value for a fixed ts value.
class cls_vs_truth_data {
private:
    vector<double> ptruth_values, pcls_values, pcls_uncertainties_n0, pcls_uncertainties_n;
public:
    cls_vs_truth_data(const vector<double> & truth_values_, const vector<double> & cls_values_, const vector<double> & cls_uncertainties_n0_,
       const vector<double> & cls_uncertainties_n_): ptruth_values(truth_values_),
       pcls_values(cls_values_), pcls_uncertainties_n0(cls_uncertainties_n0_), pcls_uncertainties_n(cls_uncertainties_n_){
    }
    const vector<double> & truth_values() const{
        return ptruth_values;
    }
    const vector<double> & cls_values() const{
        return pcls_values;
    }
    
    // absolute uncertainty on the CLs values due to finite number of s+b toys ("n")
    const vector<double> & cls_uncertainties_n() const{
        return pcls_uncertainties_n;
    }
    
    // absolute uncertainty on the CLs values due to finite number of b only toys ("n0")
    const vector<double> & cls_uncertainties_n0() const{
        return pcls_uncertainties_n0;
    }
    void write_txt(ostream & out) const{
        out << "# truth cls  cls_error_n0  cls_error_n" << endl;
        for(size_t i=0; i<ptruth_values.size(); ++i){
            out << ptruth_values[i] << "   " << pcls_values[i] << "   " << pcls_uncertainties_n0[i]  << "   " << pcls_uncertainties_n[i] << endl;
        }
    }
};

ostream & operator<<(ostream & out, const map<double, double> & m){
    for(map<double, double>::const_iterator it=m.begin(); it!=m.end(); ++it){
        out << it->first << " " << it->second << "\n";
    }
    return out;
}


void binom_with_error(size_t nominator, size_t denominator, double & p, double & p_error){
    p = 1.0 * nominator / denominator;
    p_error = max(sqrt(p*(1-p) / denominator), 1.0 / denominator);
}

// container for toy outcomes, i.e., pairs of (truth, test statistic for s+b) and (truth, test statistic for b only) tuples,
// where many truth values are the same. Note that this construction allows for the test statistic to be truth-dependent.
class truth_ts_values{
private:
    map<double, multiset<double> > truth_to_ts_sb; // signal plus background test statistic values
    map<double, multiset<double> > truth_to_ts_b;  // background only test statistic values


    static double get_quantile(const map<double, multiset<double> > & truth_to_ts, double truth, double q){
        theta_assert(q>=0.0 && q < 1.0);
        map<double, multiset<double> >::const_iterator it = truth_to_ts.find(truth);
        if(it==truth_to_ts.end()) throw invalid_argument("truth_ts_value::get_quantile");
        const multiset<double> & ts_values = it->second;
        multiset<double>::const_iterator itt = ts_values.begin();
        advance(itt, q*ts_values.size());
        return *itt;
    }
public:
    // bulk insertion (I guess this has much better performance ...)
    template <class InputIterator>
    void add_points_b(double truth, InputIterator first, InputIterator last){
        truth_to_ts_b[truth].insert(first, last);
    }
    template <class InputIterator>
    void add_points_sb(double truth, InputIterator first, InputIterator last){
        truth_to_ts_sb[truth].insert(first, last);
    }

    void add_point_b(double truth, double ts){
        truth_to_ts_b[truth].insert(ts);
    }

    void add_point_sb(double truth, double ts){
        truth_to_ts_sb[truth].insert(ts);
    }

    // truth_to_to for sb and b should cnotain same truth values. If not, complain with an invalid_argument:
    void check_consistency(){
        map<double, multiset<double> >::const_iterator it_sb = truth_to_ts_sb.begin();
        map<double, multiset<double> >::const_iterator it_b = truth_to_ts_b.begin();
        for(;it_sb!=truth_to_ts_sb.end() && it_b!=truth_to_ts_b.end(); ++it_sb, ++it_b){
            if(it_sb->first != it_b->first){
                stringstream ss;
                ss << "sb has " << it_sb->first << "; while b has " << it_b->first;
                throw invalid_argument(ss.str());
            }
        }
        if(it_sb!=truth_to_ts_sb.end()){
            stringstream ss;
            ss << "sb has " << it_sb->first << "; while b has no more values";
            throw invalid_argument(ss.str());
        }
        if(it_b!=truth_to_ts_b.end()){
            stringstream ss;
            ss << "b has " << it_b->first << "; while sb has no more values";
            throw invalid_argument(ss.str());
        }
    }

    // get the test statistic quantile q (0 <= q < 1) for the given truth value for the b-only case.
    double get_ts_b_quantile(double truth, double q){
        return get_quantile(truth_to_ts_b, truth, q);
    }

    // avoid calling this function often; it is expensive
    bool is_outlier(double truth, double ts_value, double ts_epsilon) const{
        double ts_sb_median = get_quantile(truth_to_ts_sb, truth, 0.5);
        double ts_b_median = get_quantile(truth_to_ts_b, truth, 0.5);
        double ts_sb_width = get_quantile(truth_to_ts_sb, truth, 0.84) - get_quantile(truth_to_ts_sb, truth, 0.16);
        double ts_b_width = get_quantile(truth_to_ts_b, truth, 0.84) - get_quantile(truth_to_ts_b, truth, 0.16);
        // this is a 5sigma cutoff ...
        double ts_diff_cutoff = 2.5 * max(ts_b_width + ts_epsilon, ts_sb_width + ts_epsilon);
        return fabs(ts_value - ts_sb_median) > ts_diff_cutoff && fabs(ts_value - ts_b_median) > ts_diff_cutoff;
    }
    
    bool contains_truth_value(double truth) const{
        return truth_to_ts_sb.find(truth) != truth_to_ts_sb.end() && truth_to_ts_b.find(truth) != truth_to_ts_b.end();
    }

    // get all truth values we CLs infos for, i.e., all truth values we have toys for and 0.0
    set<double> truth_values() const{
        set<double> result;
        result.insert(0.0);
        for(map<double, multiset<double> >::const_iterator it=truth_to_ts_sb.begin(); it!=truth_to_ts_sb.end(); ++it){
            result.insert(it->first);
        }
        return result;
    }
    
    size_t get_n(double truth) const{
        map<double, multiset<double> >::const_iterator it = truth_to_ts_sb.find(truth);
        if(it==truth_to_ts_sb.end()) return 0;
        return it->second.size();
    }
    
    size_t get_n0(double truth) const{
        map<double, multiset<double> >::const_iterator it = truth_to_ts_b.find(truth);
        if(it==truth_to_ts_b.end()) return 0;
        return it->second.size();
    }
    
    struct cls_info{
        double clsb, clsb_uncertainty;
        double clb, clb_uncertainty;
        double cls, cls_uncertainty_n0, cls_uncertainty_n;
    };
    
    cls_info get_cls(double ts_value, double truth_value) const{
        map<double, multiset<double> >::const_iterator it_b = truth_to_ts_b.find(truth_value);
        map<double, multiset<double> >::const_iterator it_sb = truth_to_ts_sb.find(truth_value);
        if(it_b == truth_to_ts_b.end() || it_sb == truth_to_ts_sb.end()){
            throw invalid_argument("truth_ts_values::get_cls: truth_to_ts_sb or truth_to_b do not contain given truth value");
        }
        cls_info result;
        multiset<double>::const_iterator it_ts = it_sb->second.lower_bound(ts_value);
        double p, p_error;
        size_t nominator = distance(it_ts, it_sb->second.end());
        binom_with_error(nominator, it_sb->second.size(), p, p_error);
        result.clsb = 1 - p;
        result.clsb_uncertainty = p_error;
        it_ts = it_b->second.lower_bound(ts_value);
        double p0, p0_error;
        nominator = distance(it_ts, it_b->second.end());
        binom_with_error(nominator, it_b->second.size(), p0, p0_error);
        result.clb = 1 - p0;
        result.clb_uncertainty = p0_error;
        //treat the case p0 == 1.0 by regularising the denominator with p0_error
        result.cls = (1 - p) / max(1 - p0, p0_error);
        result.cls_uncertainty_n = p_error / max(1 - p0, p0_error);
        result.cls_uncertainty_n0 = result.cls * p0_error / max(1 - p0, p0_error);
        // If estimated cls value is zero, uncertainty_n0 is not acurrate. Therefore, "iterate" the error propagation, using
        // a new, regularised cls value. For usual situations where cls > 0 and cls_uncertainty small, cls_reg and cls are almost identical.
        double cls_up = (1 - p + p_error) / max((1 - p0 - p0_error), p0_error);
        double cls_down = max((1 - p - p_error) / (1 - p0 + p0_error), 0.0);
        double cls_reg = 0.5 * (cls_up + cls_down);
        result.cls_uncertainty_n0 = cls_reg * p0_error / max(1 - p0, p0_error);
        return result;
    }
    
    // truth_to_ts is a map which specifies -- for each truth value -- the test statistic value.
    // (This is important for cases in which the test statistic value depends on the truth parameter).
    // The caller has to make sure that the ts value in this map is set for all truth values, see truth_values, including
    // truth = 0.0 (!).
    //
    // For the function "cls versus truth", there are two error sources for the error on cls:
    // * the finite number of toys for truth=0  (n0), correlated between different truth values
    // * the finite number of toys for the current truth value (n) which is uncorrelated
    // the two covariance matrices (where the latter is diagonal ...) are calculated and stored separately. This is helpful for the decision on
    // where to perform more toys, at truth=0, or at the intersection ...
    cls_vs_truth_data get_cls_vs_truth(const map<double, double> & truth_to_ts) const{
         theta_assert(truth_to_ts_b.size() == truth_to_ts_sb.size());
         vector<double> truth_values;
         vector<double> cls_values;
         vector<double> cls_uncertainties_n;
         vector<double> cls_uncertainties_n0;
         truth_values.push_back(0);
         cls_values.push_back(1);
         cls_uncertainties_n.push_back(0.001);
         cls_uncertainties_n0.push_back(0.001);
         map<double, multiset<double> >::const_iterator it_b = truth_to_ts_b.begin();
         map<double, multiset<double> >::const_iterator it_sb = truth_to_ts_sb.begin();
         for(size_t i=0; it_b!=truth_to_ts_b.end(); ++it_b, ++it_sb, ++i){
            if(it_b->first != it_sb->first){
                throw invalid_argument("truth_ts_values::get_cls_vs_truth: truth_to_ts_sb and truth_to_b contain different truth values");
            }
            double t = it_b->first;
            if(t==0.0) continue;
            truth_values.push_back(t);
            map<double, double>::const_iterator it_t_ts = truth_to_ts.find(t);
            theta_assert(it_t_ts != truth_to_ts.end());
            double ts = it_t_ts->second;
            cls_info res = get_cls(ts, t);
            cls_values.push_back(res.cls);
            cls_uncertainties_n0.push_back(res.cls_uncertainty_n0);
            cls_uncertainties_n.push_back(res.cls_uncertainty_n);
         }
         return cls_vs_truth_data(truth_values, cls_values, cls_uncertainties_n0, cls_uncertainties_n);
    }
};


template<typename T>
const std::auto_ptr<std::ostream> & operator<<(const std::auto_ptr<std::ostream> & p_out, const T & t){
    if(p_out.get()){
        (*p_out) << t;
    }
    return p_out;
}

void flush(const std::auto_ptr<std::ostream> & p_out){
    if(p_out.get()){
        (*p_out) << std::flush;
    }
}



struct fitexp_result{
    double limit, limit_error;
    fitexp_result():limit(NAN), limit_error(NAN){}
};

struct fitexp_parameters{
    Minimizer & minimizer;
    ParId pid_limit, pid_lambda;
    fitexp_parameters(Minimizer & min_, const ParId & limit_, const ParId & lambda_): minimizer(min_), pid_limit(limit_), pid_lambda(lambda_){}
};

// fit from xmin to (including) xmax. xmin and xmax must be truth_values in data.
fitexp_result fitexp(const cls_vs_truth_data & data, double target_cls, fitexp_parameters & pars, double xmin, double xmax, std::auto_ptr<std::ostream> & debug_out){
    theta_assert(xmax > xmin);
    const vector<double> & truth_values = data.truth_values();
    const vector<double> & cls_values = data.cls_values();
    const vector<double> & cls_uncertainties_n0 = data.cls_uncertainties_n0();
    const vector<double> & cls_uncertainties_n = data.cls_uncertainties_n();
    
    size_t imin = find(truth_values.begin(), truth_values.end(), xmin) - truth_values.begin();
    size_t imax = find(truth_values.begin(), truth_values.end(), xmax) - truth_values.begin();
    if(imin==truth_values.size() || imax==truth_values.size()){
        throw invalid_argument("fitexp: xmin, xmax not found in data");
    }
    if(imin >= imax){
        throw invalid_argument("fitexp: xmin, xmax define empty range");
    }
    const size_t N_range = imax - imin + 1;
    vector<double> x_values(truth_values.begin() + imin, truth_values.begin() + imax + 1);
    vector<double> y_values(cls_values.begin() + imin, cls_values.begin() + imax + 1);
    debug_out << "fitexp: fitting to imin, imax = (" << imin << ", " << imax << ");  xmin, xmax = "  << x_values[0] << ", " << x_values[x_values.size()-1] << ")\n";
    Matrix cov_inv(N_range, N_range);
    // build inverse covariance matrix:
    for(size_t i=0; i<N_range; ++i){
        cov_inv(i, i) = 1.0 / (pow(cls_uncertainties_n0[imin+i], 2) + pow(cls_uncertainties_n[imin+i], 2));
    }
    // use first and last point in range to determine initial lambda value; take care not to use log(0), so always add the
    // uncertainty which is always > 0.
    double lambda0 = min((log(y_values[N_range-1] + sqrt(1.0 / cov_inv(N_range-1, N_range-1))) - log(y_values[0] + sqrt(1 / cov_inv(0,0)))) / (x_values[N_range-1] - x_values[0]), -1e-10);
    double limit0 = 0.5 * (x_values[0] + x_values[N_range-1]);
    debug_out << "fitexp: limit0, lambda0 = " << limit0 << ", " << lambda0 << "\n";
    exp_nll nll(pars.pid_lambda, pars.pid_limit, target_cls, x_values, y_values, cov_inv);
    ParValues start, step;
    start.set(pars.pid_lambda, lambda0);
    start.set(pars.pid_limit, limit0);
    step.set(pars.pid_lambda, fabs(lambda0) * 0.05);
    step.set(pars.pid_limit, fabs(x_values.back() - x_values[0]) * 0.05);
    const double inf = numeric_limits<double>::infinity();
    map<ParId, pair<double, double> > ranges;
    ranges[pars.pid_lambda] = make_pair(-inf, 0.0);
    ranges[pars.pid_limit] = make_pair(0.0, inf);
    fitexp_result result;
    result.limit = result.limit_error = NAN;
    try{
        MinimizationResult res = pars.minimizer.minimize(nll, start, step, ranges);
        result.limit = res.values.get(pars.pid_limit);
        result.limit_error = 0.5 * (res.errors_plus.get(pars.pid_limit) + res.errors_minus.get(pars.pid_limit));
    }
    catch(MinimizationException & ex){
        debug_out << "fitexp_range: MinimizationException: " << ex.message << "\n";
        //result will be NAN
    }
    return result;
}

data_filler::data_filler(const Configuration & cfg, const boost::shared_ptr<Model> & model_, const ParId & par, const t_toy_mode & mode_, const Data & real_data_,
                         const boost::shared_ptr<Minimizer> & minimizer_, std::auto_ptr<ostream> & debug_out_): RandomConsumer(cfg, "source"),
    DataSource("source", cfg.pm->get<ProductsSink>()), model(model_), truth_parameter(par), debug_out(debug_out_), mode(mode_), real_data(real_data_), minimizer(minimizer_), start_step_ranges_init(false){
    truth_column = products_sink->declare_product(*this, "truth", theta::typeDouble);
    if(mode==ttm_datafit){
        theta_assert(model->get_observables() == real_data.get_observables());
    }
    vm = cfg.pm->get<VarIdManager>();
}

void data_filler::fill(Data & dat){
    products_sink->set_product(truth_column, truth_value);
    dat.reset();
    ParValues values;
    Random & rnd = *rnd_gen;
    if(mode==ttm_prior){
        model->get_parameter_distribution().sample(values, rnd);
    }
    else{
        map<double, ParValues>::const_iterator it = truth_to_nuisancevalues.find(truth_value);
        if(it != truth_to_nuisancevalues.end()){
            values.set(it->second);
        }
        else{
            std::auto_ptr<NLLikelihood> nll = model->get_nllikelihood(real_data);
            if(not start_step_ranges_init){
                const Distribution & d = nll->get_parameter_distribution();
                d.mode(start);
                ranges.set_from(d);
                step.set(asimov_likelihood_widths(*model, boost::shared_ptr<theta::Distribution>()));
                step.set(truth_parameter, 0.0);
                start_step_ranges_init = true;
            }
            ParValues mystart(start);
            map<double, ParValues>::const_iterator it_higher = truth_to_nuisancevalues.lower_bound(truth_value);
            map<double, ParValues>::const_iterator it_lower = it_higher;
            if(it_lower != truth_to_nuisancevalues.begin()) --it_lower;
            if(it_higher == truth_to_nuisancevalues.end()){
                if(it_lower != truth_to_nuisancevalues.end()){
                    mystart.set(it_lower->second);
                }
            }else{
                if(it_lower == it_higher){
                    mystart.set(it_lower->second);
                }
                else{
                    //TODO: could interpolate here ...
                    mystart.set(it_lower->second);
                }
            }
            mystart.set(truth_parameter, truth_value);
            ranges.set(truth_parameter, make_pair(truth_value, truth_value));
            try{
                MinimizationResult minres = minimizer->minimize(*nll, mystart, step, ranges);
                truth_to_nuisancevalues[truth_value].set(minres.values);
                values.set(minres.values);
                if(debug_out.get()){
                    debug_out << "datafit minimization for truth_value=" << truth_value << ": ";
                    const ParIds & pids = model->get_parameters();
                    for(ParIds::const_iterator pit=pids.begin(); pit!=pids.end(); ++pit){
                        debug_out << vm->get_name(*pit) << " = " <<  values.get(*pit) << "; ";
                    }
                    debug_out << "\n";
                }
            }
            catch(Exception & ex){
                debug_out << "datafit minimization for truth_value=" << truth_value << " failed.\n";
                ex.message = "in datafiller::fill: " + ex.message;
                throw;
            }
        }
    }
    values.set(truth_parameter, truth_value);
    const Distribution * rvobsdist = model->get_rvobservable_distribution();
    if(rvobsdist){
        rvobsdist->sample(values, rnd);
        dat.set_rvobs_values(ParValues(values, model->get_rvobservables()));
    }
    DataWithUncertainties data_wu;
    model->get_prediction(data_wu, values);
    const ObsIds & observables = model->get_observables();
    for(ObsIds::const_iterator it=observables.begin(); it!=observables.end(); ++it){
        dat[*it] = randomize_gauss(data_wu[*it], rnd);
        randomize_poisson(dat[*it], rnd);
    }
}


void cls_limits::run_single_truth(double truth, bool bkg_only, int n_event){
    if(!isfinite(truth)) throw invalid_argument("run_single_truth: truth not finite");
    ++runid;
    //debug_out << "starting run_single_truth(truth = " << truth << ", bkg_only = " << bkg_only << ", n_event = " << n_event << ")\n";
    // if the producer is parameter dependent, pass the truth value to the producer:
    if(pp_producer){
        ParValues values;
        values.set(truth_parameter, truth);
        pp_producer->set_parameter_values(values);
    }
    else if(bkg_only && input_bonly_ts_pool.size() > 0){
        int n = min<int>(n_event, input_bonly_ts_pool.size());
        tts->add_points_b(truth, input_bonly_ts_pool.begin(), input_bonly_ts_pool.begin() + n);
        input_bonly_ts_pool.erase(input_bonly_ts_pool.begin(), input_bonly_ts_pool.begin() + n);
        n_event -= n;
    }
    if(n_event == 0) return;
    Data data;
    vector<double> ts_values;
    ts_values.reserve(n_event);
    for (int eventid = 1; eventid <= n_event; eventid++) {
        if(stop_execution) break;
        ++n_toys;
        source->set_truth_value(bkg_only?0.0:truth);
        source->fill(data);
        bool error = false;
        try {
            producer->produce(data, *model);
        } catch (Exception & ex) {
            error = true;
            std::stringstream ss;
            ss << "Producer '" << producer->get_name() << "' failed: " << ex.message << ".";
            logtable->append(runid, eventid, LogTable::error, ss.str());
            ++n_toy_errors;
        }
        catch(std::logic_error & f){
              stringstream ss;
              ss << "Producer '" << producer->get_name() << "': " << f.what();
              throw logic_error(ss.str());
        }
        if(!error){
            if(sdc.get()) ts_values.push_back(sdc->consum_value());
            products_table->add_row(runid, eventid);
        }
        if(progress_listener){
            progress_listener->progress(n_toys, n_toys_total, n_toy_errors);
        }
    }
    if(tts.get()){
        if(bkg_only){
           tts->add_points_b(truth, ts_values.begin(), ts_values.end());
        }
        else{
           tts->add_points_sb(truth, ts_values.begin(), ts_values.end());
        }
        // check success rate and fail in case this is < 0.8. Do not check in case of stop_execution flag set:
        if(ts_values.size() * 1.0 / n_event < 0.8 && !stop_execution){
            throw logic_error("cls_limits: ts_producer fails in more than 20% of the cases");
        }
    }
}

class OutlierException{
public:
    string message;
    explicit OutlierException(const string & message_): message(message_){};
};


// Make toys at the given truth value until CLs limit accuracy reached, depending on the mode.
//
// run_single_truth_adaptive is THE place where new truth values are added.
// Apart from initialization, it is also the only place which calls run_single_truth.
//
// propagates the OutlierExcpetion from update_truth_to_ts, if the truth value is a new one.
// throws an OutlierException setting the exception limit to +infinity if clb + clb_uncertainty is smaller than clb_cutoff
//
// In mode==0 (seeding), the stopping criterion is that CLs uncertainty should be smaller than tol_cls OR CLs value is > 3sigma
// away from target cls. Calling the method twice with the same truth and ts values does nothing(!).
// In mode==1 (improvement), the number of toys at the given value is doubled, but not more than 2000 are made.
// If no toys are at this value, the behaviour is the same as mode==0.
void cls_limits::run_single_truth_adaptive(map<double, double> & truth_to_ts, double ts_epsilon, double truth, int idata, e_mode mode){
    if(!tts->contains_truth_value(truth)){
        run_single_truth(truth, true, 200);
        run_single_truth(truth, false, 200);
        update_truth_to_ts(truth_to_ts, ts_epsilon, idata);
        if(stop_execution) return;
        mode = M_ADAPTIVE;
    }
    double ts_value = truth_to_ts[truth];
    const int max_kk = idata == 0 ? 250 : 100;
    for(int kk=0; kk<max_kk; ++kk){
        if(stop_execution) return;
        truth_ts_values::cls_info res = tts->get_cls(ts_value, truth);
        double cls_uncertainty = sqrt(pow(res.cls_uncertainty_n0, 2) + pow(res.cls_uncertainty_n, 2));
        debug_out << "run_single_truth_adaptive (iteration " << kk << "; truth = " << truth << "; idata = " << idata << "; ts = " << ts_value << "):\n";
        debug_out << " clsb = " << res.clsb << " +- " << res.clsb_uncertainty << "\n";
        debug_out << " clb  = " << res.clb  << " +- " << res.clb_uncertainty << "\n";
        debug_out << " cls  = " << res.cls  << " +- " << cls_uncertainty << "\n";
        if(res.clb + res.clb_uncertainty < clb_cutoff){
            debug_out << "run_single_truth_adaptive: clb is below cutoff!\n";
            if(idata!=0){
                throw OutlierException("run_single_truth_adaptive: clb too low");
            }
        }
        if(mode == M_ADAPTIVE){
            if(fabs(res.cls - (1 - cl)) / cls_uncertainty > 3) return;
            if(cls_uncertainty < tol_cls) return;
            double target_cls_uncertainty = max(tol_cls, fabs(res.cls - (1 - cl))/3) * 0.7;
            if(res.cls_uncertainty_n0 > target_cls_uncertainty){
                size_t n0 = tts->get_n0(truth);
                size_t n0_new = pow(res.cls_uncertainty_n0 / target_cls_uncertainty, 2) * n0;
                theta_assert(n0_new >= n0);
                size_t n_toys = max<size_t>(50, min<size_t>(n0_new - n0, 500));
                debug_out << "run_single_truth_adaptive, mode 0: making " << n_toys << " more toys for background only\n";
                flush(debug_out);
                run_single_truth(truth, true, n_toys);
                if(stop_execution) return;
            }
            if(res.cls_uncertainty_n > target_cls_uncertainty){
                size_t n = tts->get_n(truth);
                size_t n_new = pow(res.cls_uncertainty_n / target_cls_uncertainty, 2) * n;
                theta_assert(n_new >= n);
                size_t n_toys = max<size_t>(50, min<size_t>(n_new - n, 500));
                debug_out << "run_single_truth_adaptive, mode 0: making " << n_toys << " more toys for signal + background\n";
                flush(debug_out);
                run_single_truth(truth, false, n_toys);
                if(stop_execution) return;
            }
        }
        else{
            //make some more toys for s+b or b-only, depending on which helps more, given the current uncertainty estimate on CLs...
            if(res.cls_uncertainty_n0 > res.cls_uncertainty_n){
                size_t n = min<size_t>(tts->get_n0(truth), 2000);
                debug_out << "run_single_truth_adaptive, mode 1: making " << n << " more toys for background only\n";
                flush(debug_out);
                run_single_truth(truth, true, n);
                if(stop_execution) return;
            }
            else{
                size_t n = min<size_t>(tts->get_n(truth), 2000);
                debug_out << "run_single_truth_adaptive, mode 1: making " << n << " more toys for signal + background\n";
                flush(debug_out);
                run_single_truth(truth, false, n);
                if(stop_execution) return;
            }
            return;
        }
    }
    throw OutlierException("run_single_truth_adaptive: too many iterations to reach CLs accuracy");
}


cls_limits::cls_limits(const Configuration & cfg): vm(cfg.pm->get<VarIdManager>()), truth_parameter(vm->get_par_id(cfg.setting["truth_parameter"])),
  runid(1), n_toys(0), n_toy_errors(0), n_toys_total(-1), pid_limit(vm->create_par_id("__limit")), pid_lambda(vm->create_par_id("__lambda")), 
  expected_bands(1000), limit_hint(NAN, NAN), reltol_limit(0.05), tol_cls(0.02), clb_cutoff(0.01), cl(0.95), beta_signal_expected(0.0){
    Setting s = cfg.setting;

    if(s.exists("debuglog")){
        string fname = s["debuglog"];
        debug_out.reset(new ofstream(fname.c_str()));
    }

    string s_mode = "set_limits";
    if(s.exists("mode")){
        s_mode = static_cast<string>(s["mode"]);
    }
    if(s_mode == "set_limits") mode = m_set_limits;
    else if(s_mode == "generate_grid") mode = m_generate_grid;
    else throw ConfigurationException("unknown mode '" + s_mode +"'");
    
    if(s.exists("beta_signal_expected")){
        beta_signal_expected = s["beta_signal_expected"];
    }

    truth_max = std::numeric_limits<double>::infinity();
    
    //1. setup common stuff:
    db = PluginManager<Database>::build(Configuration(cfg, s["output_database"]));

    std::auto_ptr<Table> logtable_underlying = db->create_table("log");
    logtable.reset(new LogTable(logtable_underlying));
    logtable->set_loglevel(LogTable::warning);
    
    std::auto_ptr<Table> rndinfo_table_underlying = db->create_table("rndinfo");
    rndinfo_table.reset(new RndInfoTable(rndinfo_table_underlying));
    cfg.pm->set("default", rndinfo_table);
    
    std::auto_ptr<Table> products_table_underlying = db->create_table("products");
    products_table.reset(new ProductsTable(products_table_underlying));
    string colname;
    if(mode == m_set_limits){
        colname = static_cast<string>(s["ts_column"]);
        sdc.reset(new SaveDoubleColumn(colname));
        std::vector<boost::shared_ptr<ProductsSink> > sinks;
        sinks.push_back(products_table);
        sinks.push_back(sdc);
        cfg.pm->set<ProductsSink>("default", boost::shared_ptr<ProductsSink>(new MultiplexingProductsSink(sinks)));
    }
    else{
        cfg.pm->set<ProductsSink>("default", products_table);
    }
    
    boost::shared_ptr<int> ptr_runid(new int(runid));
    cfg.pm->set("runid", ptr_runid);

    model = PluginManager<Model>::build(Configuration(cfg, s["model"]));

    producer = PluginManager<Producer>::build(Configuration(cfg, s["producer"]));
    pp_producer = dynamic_cast<ParameterDependentProducer*>(producer.get()); // either 0 or identical to producer ...
    if(pp_producer){
        ParIds pars = pp_producer->get_parameters();
        if(pars.size() == 0){
           // not parameter dependent after all:
           pp_producer = 0;
        }
        else{
            if(pars.size() > 1 or not pars.contains(truth_parameter)){
                 throw ConfigurationException("cls_limits only works for Producers which depend on the truth parameter only (or have no parameter dependence at all)");
            }
        }
    }

    // A.
    if(mode == m_set_limits){
        cls_limits_table = db->create_table("cls_limits");
        cls_limits__index = cls_limits_table->add_column("index", typeInt);
        cls_limits__limit = cls_limits_table->add_column("limit", typeDouble);
        cls_limits__limit_uncertainty = cls_limits_table->add_column("limit_uncertainty", typeDouble);
        minimizer = PluginManager<Minimizer>::build(Configuration(cfg, s["minimizer"]));
        tts.reset(new truth_ts_values());
        if(s.exists("truth_max")){
            truth_max = s["truth_max"];
        }
        if(s.exists("expected_bands")){
            expected_bands = s["expected_bands"];
        }
        
        // for data_source and data_source_expected: redirect the products of these to
        //   nowhere:
        boost::shared_ptr<ProductsSink> sink = cfg.pm->get<ProductsSink>();
        cfg.pm->set("default", boost::shared_ptr<ProductsSink>(new BlackholeProductsSink()));
        if(s.exists("data_source")){
            data_source = PluginManager<DataSource>::build(Configuration(cfg, s["data_source"]));
            data_source->fill(data_source_data);
        }
        if(s.exists("data_source_expected")){
            data_source_expected = PluginManager<DataSource>::build(Configuration(cfg, s["data_source_expected"]));
        }
        cfg.pm->set("default", sink);
        
        if(s.exists("reltol_limit")) reltol_limit = s["reltol_limit"];
        if(s.exists("limit_hint")){
            limit_hint.first = s["limit_hint"][0];
            limit_hint.second = s["limit_hint"][1];
        }
        if(s.exists("cl")){
            cl = s["cl"];
            if(cl >= 0.999 || cl <= 0) throw ConfigurationException("invalid value for cl. Valid range is (0, 0.999)");
        }
        if(s.exists("tol_cls")){
            tol_cls = s["tol_cls"];
        }
        if(s.exists("clb_cutoff")) {
            clb_cutoff = s["clb_cutoff"];
        }
        if(s.exists("reuse_toys")){
            Setting s2 = s["reuse_toys"];
            input_database = PluginManager<DatabaseInput>::build(Configuration(cfg, s2["input_database"]));
            input_ts_colname = colname;
            if(s2.exists("truth_column")){
                input_truth_colname = static_cast<string>(s2["truth_column"]);
            }
            else{
                input_truth_colname = "source__truth";
            }
            if(s2.exists("poi_column")){
                input_poi_colname = static_cast<string>(s2["poi_column"]);
            }
            else{
                input_poi_colname =  producer->get_name() + "__poi";
            }
        }
    }
    //B. mode = generate_grid
    else{
        theta_assert(mode == m_generate_grid);
        truth_range.first = s["truth_range"][0];
        truth_range.second = s["truth_range"][1];
        if(truth_range.first >= truth_range.second) throw ConfigurationException("invalid truth range");
        n_truth = s["n_truth"];
        if(n_truth < 2) throw ConfigurationException("n_truth < 2 not allowed");
        n_sb_toys_per_truth = s["n_sb_toys_per_truth"];
        n_b_toys_per_truth = s["n_b_toys_per_truth"];
        if(n_sb_toys_per_truth <= 0 or n_b_toys_per_truth <= 0) throw ConfigurationException("toy_per_truth <= 0 not allowed");
        n_toys_total = n_truth * (n_sb_toys_per_truth + n_b_toys_per_truth);
    }
    
    data_filler::t_toy_mode toy_mode = data_filler::ttm_prior;
    if(cfg.setting.exists("nuisancevalues-for-toys")){
        string s = cfg.setting["nuisancevalues-for-toys"];
        if(s=="prior"){}
        else if(s=="datafit"){
            toy_mode = data_filler::ttm_datafit;
        }
        else{
            throw ConfigurationException("invalid setting nuisancevalues-for-toys=\"" + s + "\"");
        }
    }
    
    source.reset(new data_filler(cfg, model, truth_parameter, toy_mode, data_source_data, minimizer, debug_out));
}


bool cls_is_significantly_larger(const cls_vs_truth_data & data, size_t i, double target_cls){
    if(i==0) return true;
    const vector<double> & cls_values = data.cls_values();
    const vector<double> & cls_uncertainties_n = data.cls_uncertainties_n();
    const vector<double> & cls_uncertainties_n0 = data.cls_uncertainties_n0();
    return cls_values[i] > target_cls && fabs(cls_values[i] - target_cls) / sqrt(pow(cls_uncertainties_n[i], 2) + pow(cls_uncertainties_n0[i], 2)) > 3;
}

bool cls_is_significantly_smaller(const cls_vs_truth_data & data, size_t i, double target_cls, double tol_cls){
    const vector<double> & cls_values = data.cls_values();
    const vector<double> & cls_uncertainties_n = data.cls_uncertainties_n();
    const vector<double> & cls_uncertainties_n0 = data.cls_uncertainties_n0();
    // cls == 0 is a special case as the error estimate for CLs is not very good.
    if(cls_values[i]==0 && sqrt(pow(cls_uncertainties_n[i], 2) + pow(cls_uncertainties_n0[i], 2)) < tol_cls) return true;
    return cls_values[i] < target_cls && fabs(cls_values[i] - target_cls) / sqrt(pow(cls_uncertainties_n[i], 2) + pow(cls_uncertainties_n0[i], 2)) > 3;
}


// find an interval seed within [i_low_min, i_high_max] of cls vs. truth in data.
// An interval seed has
// * a significantly smaller CLs value for large truth
// * a significantly larger CLs value for small truth
// if this can be fulfilled within [i_low_min, i_high_max].
pair<size_t, size_t> find_seed(const cls_vs_truth_data & data, double cl, double tol_cls, size_t i_low_min, size_t i_high_max, std::auto_ptr<ostream> & debug_out){
    theta_assert(i_high_max < data.cls_values().size());    
    // make a status vector: 1 = larger, -1 = smaller, 0 = none of these two.
    vector<int> status(data.cls_values().size());
    for(size_t i=0; i<status.size(); ++i){
        if(cls_is_significantly_larger(data, i, 1-cl)){
            status[i] = 1;
        }
        else if(cls_is_significantly_smaller(data, i, 1-cl, tol_cls)){
            status[i] = -1;
        }
        else{
            status[i] = 0;
        }
    }
    // check and insert points, s.t. algorithm terminates:
    if(status[i_low_min]!=1){
        debug_out << "WARNING: lower interval seeding index " << i_low_min << " is a point not significantly larger than target CLs!\n";
        status[i_low_min] = 1;
    }
    if(status[i_high_max]!=-1){
        debug_out << "WARNING: upper interval seeding index " << i_high_max << " is a point not significantly smaller than target CLs!\n";
        status[i_high_max] = -1;
    }
        
    // for a typical CLs versus truth curve, we should have something like
    // 1, 0, 1, 0, 0, -1, -1, 0, -1
    // i.e. first the 1s, then the -1s, with zeros in between.
    //
    // So the shortest seed in a usual case is given by [last 1, first -1].
    //
    // In case there are outliers, like in
    // 1, 1, -1, 1, 1, 1, 0, 0, -1, -1, -1
    // this fails. In such a case, we choose to ignore the single -1 between the 1s. This is prefered
    // over ignoring the three 1s after this -1, as this changes only one point.
    //
    // So the algorithm is: make the pair (last 1, first -1). Ideally this is an interval (=second number greater than first).
    // If not, look at the points within the interval and let the majority decide which one is the "outlier" (the 1 or the -1, or both),
    // and repeat the search, ignoring the outlier(s).
    size_t last_1 = i_high_max;
    size_t first_m1 = i_low_min;
    while(status[last_1]!=1) --last_1;
    while(status[first_m1]!=-1) ++ first_m1;
    // while there is an outlier: remove it ...
    while(last_1 > first_m1){
        size_t n1 = 0, nm1 = 0;
        for(size_t i = first_m1+1; i < last_1; ++i){
            if(status[i]==-1) ++nm1;
            else if(status[i]==1) ++n1;
        }
        // note that for n1 == nm1, both are outliers ...
        if(nm1 >= n1){ // 1 is the outlier:
            --last_1;
            while(status[last_1]!=1) --last_1;
        }
        if(n1 >= nm1){ // -1 is the outlier:
            ++first_m1;
            while(status[first_m1]!=-1) ++ first_m1;
        }
    }
    return make_pair(last_1, first_m1);
}

// update the truth_to_ts map to contain ts values for all truth_values in truth_to_ts using current_data.
// throws an OutlierException if the ts value is too extreme to be considered for CLs limits.
// All calling functions should catch this exception
void cls_limits::update_truth_to_ts(map<double, double> & truth_to_ts, double ts_epsilon, int idata){
    set<double> truth_values = tts->truth_values();
    boost::optional<double> constant_ts; //filled only in case of constant, non-poi dependent ts, i.e., pp_producer == 0
    // in case of data and constant ts, check whether we have calculated the ts value already:
    if(pp_producer==0 && truth_to_ts.size() > 0){
        constant_ts = truth_to_ts.begin()->second;
    }
    for(set<double>::const_iterator it = truth_values.begin(); it!=truth_values.end(); ++it){
        // in case we already know the (constant) ts, it's easy:
        if(constant_ts){
            truth_to_ts[*it] = *constant_ts;
            continue;
        }
        // if we have it already, do not calculate it again:
        if(truth_to_ts.find(*it)!=truth_to_ts.end()){
             continue;
        }
        // Now, calculate ts value, save it in truth_to_ts and save it in the table:
        if(pp_producer){
            ParValues vals;
            vals.set(truth_parameter, *it);
            pp_producer->set_parameter_values(vals);
        }
        try {
            producer->produce(current_data, *model);
        } catch (Exception & ex) {
            ex.message += " (while running ts_producer for current_data)";
            throw;
        }
        double ts = sdc->consum_value() + ts_epsilon;
        // for truth!=0, and if this is not real data, check if ts is an outlier and if it is, reject this toy by throwing an OutlierException:
        if(*it!=0){
            if(tts->is_outlier(*it, ts, ts_epsilon)){
                debug_out << "update_truth_to_ts: ts for truth = " << *it << " is " << ts << " which was identified as outlier.\n";
                if(idata!=0){
                    throw OutlierException("update_truth_to_ts: test statistic outlier");
                }
            }
        }
        truth_to_ts[*it] = ts;
        //call source->fill to set products column
        Data dummy_data;
        source->set_truth_value(*it);
        source->fill(dummy_data);
        products_table->add_row(0, idata);
        // if ts value is constant, fill it to constant_ts, for the next iteration:
        if(pp_producer==0){
            constant_ts = ts;
        }
    }
    correct_truth_to_ts(truth_to_ts);
}

void cls_limits::read_reuse_toys(){
    vector<string> colnames;
    colnames.push_back("runid");
    colnames.push_back(input_truth_colname);
    colnames.push_back(input_ts_colname);
    if(pp_producer){
        colnames.push_back(input_poi_colname);
    }
    std::auto_ptr<DatabaseInput::ResultIterator> it = input_database->query("products", colnames);
    set<double> truth_values;
    for(; it->has_data(); ++(*it)){
        int runid = it->get_int(0);
        if(runid<=0) continue;
        double truth = it->get_double(1);
        double ts = it->get_double(2);
        double poi = pp_producer ? it->get_double(3) : -1.0;
        // if the truth value was > 0, it's clear it was an s+b toy.
        if(truth > 0){
            truth_values.insert(truth);
            tts->add_point_sb(truth, ts);
        }
        else{
            // if truth is == 0.0, it's a background-only toy, but the situation is not so clear to which
            // truth value to add the point, unless we have a poi value :
            if(pp_producer){
                if(poi > 0)
                   tts->add_point_b(poi, ts);
                // so poi <= 0 are skipped
            }
            else{
                input_bonly_ts_pool.push_back(ts);
            }
        }
    }
    
    // in case of no parameter dependece, use at least a few b-only toys for each truth value (cls_limits requires to have
    // some b-only toys for each truth value!)
    if(pp_producer == 0){
        int n_truth = truth_values.size();
        int n_bonly_per_truth = input_bonly_ts_pool.size() / n_truth;
        debug_out << "reuse_toys: adding 200 b-only toys for each of the " << truth_values.size() << " truth values. Have " << input_bonly_ts_pool.size() << " toys, so ";
        if(n_bonly_per_truth < 200){
            debug_out << "they will suffice.\n";
        }
        else{
            debug_out << "will have to make some additional ones.\n";
        }
        flush(debug_out);
        for(set<double>::const_iterator it=truth_values.begin(); it!=truth_values.end(); ++it){
            // we have the code for pulling values from the pool already as long as available, no need to recode ...
            run_single_truth(*it, true, 200);
        }
    }
    tts->check_consistency();
}

void cls_limits::run(){
    switch(mode){
       case m_set_limits: run_set_limits(); break;
       case m_generate_grid: run_generate_grid(); break;
    }
}

void cls_limits::run_generate_grid(){
    double truth_step = (truth_range.second - truth_range.first) / (n_truth - 1);
    for(int i=0; i<n_truth; ++i){
        double truth = truth_range.first + i*truth_step;
        run_single_truth(truth, true, n_b_toys_per_truth);
        run_single_truth(truth, false, n_sb_toys_per_truth);
    }
}


void cls_limits::correct_truth_to_ts(std::map<double, double> & truth_to_ts){
    theta_assert(truth_to_ts.size() > 0);
    double ts_first = truth_to_ts.begin()->second;
    double ts_last = truth_to_ts.rbegin()->second;
    // sign is -1 if falling, +1 if raising (and kind of arbitray if constant):
    double sign = (ts_first - ts_last > 0) ? -1 : 1;
    for(map<double, double>::iterator it = truth_to_ts.begin(); it!= truth_to_ts.end(); ++it){
        map<double, double>::iterator it_next = it;
        ++it_next;
        if(it_next == truth_to_ts.end()) break;
        double diff = (it_next->second - it->second) * sign;
        if(diff < 0){
            debug_out << "correct_truth_to_ts: ts as a function of truth is not monotonic; (truth, ts) = (" << it->first << ", "
                      << it->second << "), (" << it_next->first << ", " << it_next->second << ")\n";
            // correct it_next value by setting it to the same value as it:
            it_next->second = it->second;
        }
    }
}

void cls_limits::run_set_limits(){
    debug_out << "tol_cls = " << tol_cls << ".\n";
    // 0. determine signal width
    double signal_width = limit_hint.second - limit_hint.first;
    if(!std::isfinite(limit_hint.first) || !std::isfinite(limit_hint.second)){
        ParValues widths = asimov_likelihood_widths(*model, boost::shared_ptr<Distribution>());
        signal_width = widths.get(truth_parameter);
        debug_out << "signal_width = " << signal_width << "\n";
        if(signal_width <= 0.0){
            debug_out << "Warning: signal_width <=0.0. Setting signal_width to 1.0.\n";
            signal_width = 1.0;
        }
        flush(debug_out);
    }

    //1. determine ts_epsilon. This is important for degenerate test statistic (such as likelihood ratios where
    // the backgruond only prior fixes the signal and is free in the s+b prior), because mathematically degenerate
    // values will not be necessarily exactly degenerate numerically. A possible solution to that is
    // to add small values to all 'observed' test statistic values, i.e., make them more signal-like (although
    // only by a very small amount).
    // truth0 is a guess for a 'high' truth value, which should be in about the right region of the limit.
    double truth0 = 2 * signal_width;
    if(std::isfinite(limit_hint.second)) truth0 = max(truth0, limit_hint.second);
    run_single_truth(truth0, false, 200);
    run_single_truth(truth0, true, 200);
    const double ts_epsilon = fabs(tts->get_ts_b_quantile(truth0, 0.68) - tts->get_ts_b_quantile(truth0, 0.16)) * 1e-3;
    debug_out << "ts_epsilon is " << ts_epsilon << "\n";

    // read in data from the input database, if set:
    if(input_database.get()){
        read_reuse_toys();
    }

    fitexp_parameters pars(*minimizer, pid_limit, pid_lambda);
    size_t n_results = 0;
    const size_t N_maxit = 200;
    double clb_cutoff0 = clb_cutoff;
    // idata == 0 means from data_source; 1..expected_baend are for the bands.
    // Note that the case idata=0 is moved to the end to have a more robust result for data
    for(int idata=1; ; ++idata){
        if(idata==expected_bands+1){
            idata=0;
        }
        if(idata==0 && data_source.get()==0) break;
        if(idata==0){
            clb_cutoff = 0.0;
        }
        else{
            clb_cutoff = clb_cutoff0;
        }
        debug_out << "starting idata " << idata  << "\n";
        flush(debug_out);
        try{
            if(idata==0){
                current_data = data_source_data;
            }
            else{
                // make a toy for the expected limit:
                if(data_source_expected.get()){
                    data_source_expected->fill(current_data);
                }
                else{
                    source->set_truth_value(beta_signal_expected);
                    source->fill(current_data);
                }
            }
            map<double, double> truth_to_ts;
            // run an update, to get the ts values at the truth points from previous idata iterations:
            update_truth_to_ts(truth_to_ts, ts_epsilon, idata);
            //3.a. find a seed:
            //3.a.i. make some adaptive toys at the highest current truth value to make it significantly smaller:
            cls_vs_truth_data data = tts->get_cls_vs_truth(truth_to_ts);
            run_single_truth_adaptive(truth_to_ts, ts_epsilon, data.truth_values().back(), idata, M_ADAPTIVE);
            if(stop_execution) return;
            data = tts->get_cls_vs_truth(truth_to_ts);
            while(not cls_is_significantly_smaller(data, data.cls_values().size()-1, 1 - cl, tol_cls)){
                 double next_value = data.truth_values().back()*2.0;
                 if(next_value > truth_max){
                     debug_out << "to high truth value: " << next_value << ". Marking as outlier\n";
                     throw OutlierException("run_set_limits: too high truth value");
                 }
                 debug_out << "making toys at high truth values to find upper limit on limit; next is truth=" << next_value << "\n";
                 flush(debug_out);
                 run_single_truth_adaptive(truth_to_ts, ts_epsilon, next_value, idata, M_ADAPTIVE);
                 if(stop_execution) return;
                 data = tts->get_cls_vs_truth(truth_to_ts);
            }
            //3.b. iterate
            fitexp_result latest_res;
            size_t i_low = 0, i_high = data.truth_values().size()-1;
            double truth_low;
            double truth_high;
            for(size_t i=0; i<=N_maxit; ++i){
                // 0. debugging stuff:
                if(i==N_maxit){
                    debug_out << idata << "." << i << ": maximum number of iterations reached. Marking as outlier.\n";
                    throw OutlierException("run_set_limits: too many iterations necessary");
                }
                debug_out << idata << "." << i << "\n";
                if(debug_out.get()){
                    debug_out << "cls vs truth:\n";
                    data.write_txt(*debug_out);
                    flush(debug_out);
                    debug_out << "truth_to_ts:\n" << truth_to_ts << "\n";
                }
                // 1. find fit range within current interval (in the first iteration, this will be in whole range):
                pair<size_t, size_t> seed = find_seed(data, cl, tol_cls, i_low, i_high, debug_out);
                i_low = seed.first;
                i_high = seed.second;
                debug_out << "proposed fit interval: " << data.truth_values()[i_low] << "--" << data.truth_values()[i_high] << "\n";
                // 1.b. make sure that fit range includes latest fit (if it converged ...):
                if(!std::isnan(latest_res.limit) && (latest_res.limit < data.truth_values()[i_low] || latest_res.limit > data.truth_values()[i_high])){
                    debug_out << "WARNING: latest fitted limit (" << latest_res.limit << ") not contained in proposed fit interval, making fit interval larger\n";
                    while(latest_res.limit < data.truth_values()[i_low] && i_low > 0){
                        --i_low;
                        while(!cls_is_significantly_larger(data, i_low, 1-cl) && i_low > 0) --i_low;
                    }
                    while(latest_res.limit > data.truth_values()[i_high] && i_high < data.truth_values().size()-1){
                        ++i_high;
                        while(!cls_is_significantly_smaller(data, i_high, 1-cl, tol_cls) && i_high < data.truth_values().size()-1) ++i_high;
                    }
                }
                truth_low = data.truth_values()[i_low];
                truth_high = data.truth_values()[i_high];
                // 2. make the fit
                latest_res = fitexp(data, 1 - cl, pars, truth_low, truth_high, debug_out);
                if(std::isnan(latest_res.limit)){
                    debug_out << "exp fit did not work; fill in some random point in fitted interval, with large error\n";
                    // u is a uniform random number between 0 and 1. Details don't play a role here, so just hard code a linear
                    // congruent generator using i+17 as seed:
                    double u = (((i + 17) * 1103515245 + 12345) % 4294967296) * 1.0 / 4294967296;
                    latest_res.limit = truth_low + (0.25 + 0.5 * u) * (truth_high - truth_low);
                    latest_res.limit_error = 0.5 * (truth_high - truth_low);
                }
                debug_out << idata << "." << i << ": fitted limit on range " << truth_low << "--" << truth_high << " using " << (i_high - i_low + 1) << " points:\n";
                debug_out << "       limit = " << latest_res.limit << " +- " << latest_res.limit_error << "\n";
                //accept result if truth=0 not included and error is small enough:
                if(latest_res.limit_error / latest_res.limit < reltol_limit && latest_res.limit >= truth_low && latest_res.limit <= truth_high && (i_high - i_low + 1) > 2){
                    debug_out << "       limit accepted.\n";
                    break;
                }
                else{
                    debug_out << "       limit not accepted.\n";
                }
                flush(debug_out);
                
                // 3. make new toys
                // 3.a. If the CLs range includes 0, make a new truth point close to 0:
                if(i_low == 0){
                    run_single_truth_adaptive(truth_to_ts, ts_epsilon, 0.5 * data.truth_values()[1], idata, M_ADAPTIVE);
                    if(stop_execution) return;
                }
                // 3.b. add a new point randomly in the interval:
                // add a random point in the fit interval:
                const double u = (((i + 17) * 1103515245 + 12345) % 4294967296) * 1.0 / 4294967296;
                const double next_truth0 = truth_low + u * (truth_high - truth_low);
                double next_truth = next_truth0;
                debug_out << "next truth for toys: " << next_truth << "\n";
                flush(debug_out);
                // "round" to a neighboring value if within the target reltol_limit:
                if((i_high - i_low + 1) > 2){
                    double min_diff = numeric_limits<double>::infinity();
                    for(size_t m=0; m < data.truth_values().size(); ++m){
                        double t = data.truth_values()[m];
                        if(fabs(t - next_truth0) < min_diff && fabs(t - next_truth0) / max(t, next_truth0) < reltol_limit){
                            next_truth = t;
                            min_diff = fabs(t - next_truth0);
                        }
                    }
                    if(!std::isinf(min_diff)){
                        debug_out << "rounded next truth value to existing point " << next_truth << "\n";
                        flush(debug_out);
                    }
                }
                run_single_truth_adaptive(truth_to_ts, ts_epsilon, next_truth, idata, M_MORE);
                if(stop_execution) return;
                
                // 4. update data and i_low, i_high, as new indices might have appeared:
                data = tts->get_cls_vs_truth(truth_to_ts);
                i_low = find(data.truth_values().begin(), data.truth_values().end(), truth_low) - data.truth_values().begin();
                i_high = find(data.truth_values().begin(), data.truth_values().end(), truth_high) - data.truth_values().begin();
                
                //5. pruning phase:
                // make interval [truth_low, truth_high] smaller if (border is far away (>2sigma) from limit or border is 0) AND the next point can serve as new interval border
                if((fabs(truth_low - latest_res.limit) / latest_res.limit_error > 2 or i_low==0) && i_low+1 < i_high){
                     debug_out << "lower border " << truth_low << " far away from current limit --> test if can be removed ...\n";
                     //re-run next candidate point:
                     run_single_truth_adaptive(truth_to_ts, ts_epsilon, data.truth_values()[i_low+1], idata, M_ADAPTIVE);
                     if(stop_execution) return;
                     if(cls_is_significantly_larger(data, i_low+1, 1 - cl)){
                         ++i_low;
                         truth_low = data.truth_values()[i_low];
                         debug_out << "--> yes, replaced by " << truth_low << ".\n";
                     }
                     else{
                         debug_out << "--> no, next point not valid.\n";
                     }
                }
                if(fabs(truth_high - latest_res.limit) / latest_res.limit_error > 2 && i_low < i_high-1){
                     debug_out << "upper border " << truth_high << " far away from current limit --> test if can be removed ...\n";
                     //re-run next candidate point:
                     run_single_truth_adaptive(truth_to_ts, ts_epsilon, data.truth_values()[i_high-1], idata, M_ADAPTIVE);
                     if(stop_execution) return;
                     if(cls_is_significantly_smaller(data, i_high-1, 1 - cl, tol_cls)){
                         --i_high;
                         truth_high = data.truth_values()[i_high];
                         debug_out << "--> yes, replaced by " << truth_high << ".\n";
                     }
                     else{
                         debug_out << "--> no, next point not valid.\n";
                     }
                }
            }
            ++n_results;
            Row r;
            r.set_column(cls_limits__index, idata);
            r.set_column(cls_limits__limit, latest_res.limit);
            r.set_column(cls_limits__limit_uncertainty, latest_res.limit_error);
            cls_limits_table->add_row(r);
            debug_out << "final limit for idata = " << idata << ": " << latest_res.limit << " +- " << latest_res.limit_error << "\n";
            flush(debug_out);
        }
        catch(OutlierException & ex){
            debug_out << "idata " << idata  << " identified as outlier (reason: " << ex.message << "); skipping\n";
        }
        catch(MinimizationException & ex){
            debug_out << "idata " << idata << " had minimizationexception; skipping.\n";
            
        }
        if(idata==0) break;
    }
    // print the total number of toys per point (hint for setting up grid ...)
    if(debug_out.get()){
        debug_out << "total toys made: " << n_toys << "; errorneous: " << n_toy_errors << "\n";
        debug_out << "number of toys per truth: truth n_b_toys n_sb_toys\n";
        set<double> truth_values = tts->truth_values();
        for(set<double>::const_iterator it=truth_values.begin(); it!=truth_values.end(); ++it){
            double t = *it;
            if(t==0) continue;
            debug_out << t << " " << tts->get_n0(t) << " " << tts->get_n(t) << "\n";
        }
        debug_out << "number of successful limits: " << n_results << " / " << (expected_bands + (data_source.get() ? 1 : 0)) << "\n";
    }
}

REGISTER_PLUGIN(cls_limits)

