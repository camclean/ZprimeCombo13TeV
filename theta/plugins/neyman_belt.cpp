#include "interface/main.hpp"
#include "interface/database.hpp"
#include "interface/plugin.hpp"
#include "interface/model.hpp"
#include "interface/distribution.hpp"

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <set>
#include <algorithm>

#include <boost/bimap/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>

using namespace std;
using namespace boost::bimaps;
using namespace theta;

class SaveDoubleProducts: public ProductsSink{
private:
    map<string, double> double_data;
    map<Column, string> column_names;
    int next_icol;
public:
    SaveDoubleProducts(): next_icol(0){}
    
    virtual Column declare_column_impl(const std::string & product_name, const data_type & type){
        Column result(next_icol++);
        column_names[result] = product_name;
        return result;
    }
    
    virtual void set_product(const Column & c, double d){
        double_data[column_names[c]] = d;
    }
    virtual void set_product(const Column & c, int i){
    }
    virtual void set_product(const Column & c, const std::string & s){
    }
    virtual void set_product(const Column & c, const Histogram1D & h){
    }
    
    double get_product(const std::string & name){
        map<string, double>::const_iterator it = double_data.find(name);
        if(it==double_data.end()) throw invalid_argument("unknown product '" + name + "'");
        return it->second;
    }
};

// information for a single toy experiment entering the neyman construction:
// the value of "truth", test statistic "ts" and "ordering":
struct tto{
   double truth, ts, ordering;

   tto(double truth_, double ts_, double ordering_=0): truth(truth_), ts(ts_), ordering(ordering_){}

   struct truth_ts_ordering{
      bool operator()(const tto & left, const tto & right) const{
          return (left.truth < right.truth) || (left.truth==right.truth && left.ts < right.ts);
      }
   };

   struct truth_o_ordering{
      bool operator()(const tto & left, const tto & right) const{
          return (left.truth < right.truth) || (left.truth==right.truth && left.ordering < right.ordering);
      }
   };
};

// information of a toy ensemble required for neyman construction. Basically
// a set of tto objects. However, provides views on these ttos both sorted
// by ts value and sorted by odering value, for a given truth value.
class tto_ensemble{
public:
    enum sortflag{ sorted_by_ts, sorted_by_ordering };

    typedef vector<tto>::const_iterator const_tto_iterator;
    typedef std::pair<const_tto_iterator, const_tto_iterator> tto_range;

    void reserve(size_t n){
        ttos__truth_ts_sorted.reserve(n);
        ttos__truth_o_sorted.reserve(n);
    }
    void add(const tto & t);
    size_t size(){
       return ttos__truth_ts_sorted.size();
    }
    tto_range get_ttos(double truth_value, const sortflag s);
    set<double> get_truth_values();

    void swap(tto_ensemble & rhs){
        std::swap(sorted, rhs.sorted);
        std::swap(ttos__truth_ts_sorted, rhs.ttos__truth_ts_sorted);
        std::swap(ttos__truth_o_sorted, rhs.ttos__truth_o_sorted);
        std::swap(truth_partitions, rhs.truth_partitions);
    }
    
private:
   void sort();
   bool sorted;

   // sort() sorts / fills all of the following objects:
   vector<tto> ttos__truth_ts_sorted;
   vector<tto> ttos__truth_o_sorted;
   map<double, size_t> truth_partitions;
};

namespace std{
   template<>
   void swap(tto_ensemble & left, tto_ensemble & right){
       left.swap(right);
   }
}


void tto_ensemble::sort(){
    if(sorted) return;
    truth_partitions.clear();
    if(ttos__truth_ts_sorted.size()==0) return;
    theta_assert(ttos__truth_ts_sorted.size() == ttos__truth_o_sorted.size());
    std::sort(ttos__truth_ts_sorted.begin(), ttos__truth_ts_sorted.end(), tto::truth_ts_ordering());
    std::sort(ttos__truth_o_sorted.begin(), ttos__truth_o_sorted.end(), tto::truth_o_ordering());
    truth_partitions.clear();
    size_t index = 0;
    double last_truth = ttos__truth_ts_sorted[0].truth - 1;
    for(vector<tto>::const_iterator it = ttos__truth_ts_sorted.begin(); it!=ttos__truth_ts_sorted.end(); ++it, ++index){
        const double truth = it->truth;
        theta_assert(ttos__truth_o_sorted[index].truth == truth);
        if(truth != last_truth){
            truth_partitions[truth] = index;
            last_truth = truth;
        }
    }
    sorted = true;
}

void tto_ensemble::add(const tto & t){
    sorted = false;
    ttos__truth_ts_sorted.push_back(t);
    ttos__truth_o_sorted.push_back(t);
}

set<double> tto_ensemble::get_truth_values(){
    sort();
    set<double> result;
    for(map<double, size_t>::const_iterator it=truth_partitions.begin(); it!=truth_partitions.end(); ++it){
        result.insert(it->first);
    }
    return result;
}


tto_ensemble::tto_range tto_ensemble::get_ttos(double truth_value, const sortflag s){
    sort();
    // get indices:
    map<double, size_t>::const_iterator it1 = truth_partitions.lower_bound(truth_value);
    map<double, size_t>::const_iterator it2 = truth_partitions.upper_bound(truth_value);
    if(it1==truth_partitions.end() || it1->first != truth_value){
        throw invalid_argument("tto_ensemble::get_ttos: unknown truth_value");
    }
    size_t index1 = it1->second;
    size_t index2 = (it2==truth_partitions.end()) ? ttos__truth_o_sorted.size() : it2->second;

    //return iterators to the appropriate vector:
    switch(s){
        case sorted_by_ts:
            return tto_range(ttos__truth_ts_sorted.begin() + index1, ttos__truth_ts_sorted.begin() + index2);
        case sorted_by_ordering:
            return tto_range(ttos__truth_o_sorted.begin() + index1, ttos__truth_o_sorted.begin() + index2);
    }
    throw invalid_argument("invalid sortflag in tto_ensemble::get_ttos");
}


/** \brief Neyman belt construction
 *
 * This Neyman construction for confidence intervals works by defining a <em>confidence belt</em>
 * on the (truth, test statistic) plane, given the result of toy experiments as input. For each toy,
 * the true value of the parameter of interest, "truth", and the value of the test statistic, "ts"
 * have to be known. The result of the toy experiments used here as input is the output of a previous theta run.
 *
 * Important: the values of the "truth" column in the input must be discrete. It is recommended to
 * use \link equidistant_deltas \endlink for the parameter of interest to make the input.
 *
 * The constructed belt will have (at least) the desired coverage on the ensemble used as input.
 *
 * \code
 * main = {
 *   type = "neyman_belt";
 *   cls = [0.68, 0.95];
 *   ordering_rule = "fclike";
 *   // fclike_options is only required if ordering_rule is "fclike".
 *   fclike_options = {
 *      model = "@some_model";
 *      ts_producer = "@mle";    // use the same producer as for the ts_column in toy_database below
 *      truth = "beta_signal";   // parameter name of the "truth"
 *      truth_n = 200;
 *      truth_range = (0.0, 20.0);
 *   };
 *   force_increasing_belt = true; // optional; default is false
 *   // input:
 *   toy_database = {
 *      type = "sqlite_database_in";
 *      filename = "toys.db";
 *   };
 *   truth_column = "source__beta_signal";
 *   ts_column = "mle__beta_signal";
 *   // output:
 *   output_database = {
 *      type = "sqlite_database";
 *      filename = "belt.db";
 *   };
 * };
 * \endcode
 *
 * \c type must always be "neyman_belt" to select this plugin
 *
 * \c cls are the confidence levels for which to construct the belts
 *
 * \c ordering_rule is the ordering rule to use for the belt construction. It controls which
 *   toys are included in the belt first. Valid values are 'central', 'central_shortest', 'lower', 'upper' and 'fclike'. If choosing
 *   'fclike', the additional setting group 'fclike_options' has to be given. It has to specify the model, the test statistic producer,
 *   the "truth" parameter. "truth_n" and "truth_range" are the number of points and range to use for the "truth" parameter
 *   to calculate asimov data for the test statistic calculation. Default is truth_n = 200 and a truth_range which
 *   is [truth_min, 2*truth_max - truth_min] where truth_min and truth_max are taken from the input toy_database.
 *
 * \c force_increasing_belt controls whether to the belt borders in "ts" as function of "truth" are forced to be monotonically
 *    increasing. This only makes sense if the test statistic increases as function of truth.
 *
 * \c toy_database is the input database in which the information about the toys to use for the construction is stored.
 *  \c truth_column and \c ts_column are the column names in this database for the value for "truth" and "ts", respectively.
 *
 * \c output_database configures the database in which the resulting belt is written.
 * 
 * Typical choices for the test statistic include the maximum likelihood estimate for the parameter of interest (via the mle plugin)
 * or a likelihood ratio (via the deltanll_hypotest plugin). However, it is also possible to use something
 * very different such as the quantiles of the Bayesian posterior. The latter has an interesting use case
 * as the Neyman construction using Bayesian quatiles as test statistic can be interpreted as coverage test
 * for the Bayesian method, at least for one-sided intervals.
 *
 * First, an ordering value is calculated for each input pair (truth, ts) from the toy_database according to the \c ordering_rule setting.
 * The belt construction works on an ensemble of toys where each toy is represented
 * by the three-tuple (truth, ts, ordering).
 *
 * The belt construction finds an interval in "ts" for each value of "truth". Apart from constraints
 * if \c force_increasing_belt is true, the interval construction for different "truth" values are independent.
 * The interval construction for a given "truth" value works in two steps:
 * (i) seed finding and (ii) interval growing. Seed finding is done by picking the
 * toy with smallest associated ordering value which respects constraints from the increasing belt
 * (if \c force_increasing_belt is true). Interval growing adds neighboring toys in ts, preferring those
 * with lower ordering value. This construction ensures that the resulting ts set belonging to the belt for a particular
 * truth value is always an interval in ts.
 * 
 * For ordering_rule 'lower' and 'upper', the resulting interval extends to +-infinity in one side.
 *
 * The result is written in a table named "belt" with (3*n + 1) columns if n is the number of confidence levels configured. The column
 * names are "truth" and for each confidence level "ts_lower<cl>", "ts_upper<cl>" and "coverage<cl>" where "<cl>" is the confidence level
 * as string using 5 digits, if placing the decimal point after the first one. For example, a confidence level 0.68 will create the columns
 * "ts_lower06800", "ts_upper06800", and "coverage06800". For each row, the ts_lower and ts_upper columns specify the "ts" interval
 * included in the belt. The coverage column contains the coverage of this interval which might be higher than
 * the configured one, especially if the test statistic is discrete (this coverage is calculated for a slightly increased
 * ts interval).
 * There is one row per input truth value.
 *
 * The reported progress is pretty arbitrary; no progress is reported for reading the input database;
 * progress is reported for ordering value calculation only if ordering_rule is "fclike". The progress for interval
 * finding is reported as per (cl, truth) value.
 */
class neyman_belt: public theta::Main{
public:
    neyman_belt(const Configuration & cfg);
    virtual void run();
    virtual ~neyman_belt(){}
private:
    vector<double> cls;
    string ordering_rule;
    bool force_increasing_belt;
    
    boost::shared_ptr<DatabaseInput> toy_database;
    string truth_column, ts_column;
    
    boost::shared_ptr<Database> output_database;
    
    // fclike options:
    boost::shared_ptr<Model> model;
    std::auto_ptr<ParId> poi;
    boost::shared_ptr<Producer> ts_producer;
    boost::shared_ptr<SaveDoubleProducts> save_double_products;
    std::string ts_name;
    size_t truth_n;
    std::pair<double, double> truth_range;
    
    int progress_done, progress_total, progress_errors;
    
    void progress(){
        if(progress_listener) progress_listener->progress(++progress_done, progress_total, progress_errors);
    }
    
    void add_ordering_fclike(tto_ensemble & inp);
};

// construct a belt interval; the input ranges are for a single truth value.
// The returned interval is an interval in the test statistic.
// ts_interval0_min and ts_interval1_min are the minimum values
//   for the lower and upper interval ends in the ts value. This is
//   set to the result of the previous interval to ensure that the resulting
//   belt yields actual confidence intervals (not just regions). Can be set
//   to -infinity to not impose any restrictions.
//   Note that these restrictions are not kept if coverage would be lost.
struct t_interval_coverage{
    pair<double, double> interval;
    double coverage;
};

t_interval_coverage construct_interval(const tto_ensemble::tto_range & tto_range_o_sorted,
    const tto_ensemble::tto_range & tto_range_ts_sorted,
    double cl, double ts_interval0_min, double ts_interval1_min){
    const size_t n_total = distance(tto_range_o_sorted.first, tto_range_o_sorted.second);
    theta_assert(distance(tto_range_o_sorted.first, tto_range_o_sorted.second) == distance(tto_range_ts_sorted.first, tto_range_ts_sorted.second));
    const size_t n_target = static_cast<size_t>(cl * n_total);
    //take the ts value corresponding to the lowest ordering value which respects the restrictions as starting point for the interval:
    pair<double, double> result;
    tto_ensemble::const_tto_iterator it_seed = tto_range_o_sorted.first;
    while(it_seed->ts < ts_interval0_min && it_seed != tto_range_o_sorted.second){
        ++it_seed;
    }
    if(it_seed == tto_range_o_sorted.second){
        it_seed = tto_range_o_sorted.first;
        result.first = result.second = it_seed->ts;
    }
    else{
       result.first = result.second = it_seed->ts;
       result.first = max(result.first, ts_interval0_min);
       result.second = max(result.second, ts_interval1_min);
    }
    theta_assert(result.second >= result.first);

    // from here on, only use the "ts sorted" input, not the "o sorted" one.
    // add ajacent toys directly left of or right of the result interval, preferring the value according to the ordering.
    // Repeat until coverage is reached.
    // Interval contains the values in [it_min, it_max), i.e., excludes where it_max points to; in particular,
    //   it_max might point beyond some input container.
    tto_ensemble::const_tto_iterator it_min = lower_bound(tto_range_ts_sorted.first,
                            tto_range_ts_sorted.second, tto(tto_range_o_sorted.first->truth, result.first), tto::truth_ts_ordering());
    tto_ensemble::const_tto_iterator it_max = upper_bound(tto_range_ts_sorted.first,
                            tto_range_ts_sorted.second, tto(tto_range_o_sorted.first->truth, result.second), tto::truth_ts_ordering());
    size_t n = distance(it_min, it_max);
    theta_assert(n > 0);
    const double nan = numeric_limits<double>::quiet_NaN();
    while(n < n_target){
        double order_value_left = nan, order_value_right = nan;
        double ts_value_left = nan, ts_value_right = nan;
        if(it_min != tto_range_ts_sorted.first){
            tto_ensemble::const_tto_iterator it_left = it_min;
            --it_left;
            order_value_left = it_left->ordering;
            ts_value_left = it_left->ts;
        }
        if(it_max != tto_range_ts_sorted.second){
            order_value_right = it_max->ordering;
            ts_value_right = it_max->ts;
        }
        bool left = false, right = false;
        if(std::isnan(order_value_right)){
           if(std::isnan(order_value_left)){
               throw Exception("no valid range found");//can only happen if input contained nans ...
           }
           else{
              --it_min; ++n; result.first = ts_value_left; left = true;
           }
        }
        else{
           if(std::isnan(order_value_left)){
              ++it_max; ++n; result.second = ts_value_right; right = true;
           }
           else{
               if(order_value_left < order_value_right && ts_value_left >= ts_interval0_min){--it_min; ++n; result.first = ts_value_left; left = true;}
               else { ++it_max; ++n; result.second = ts_value_right; right = true; }
           }
        }
        //A.3.a. Add points with same / very similar ts value of the one just added
        double ilen = result.second - result.first;
        if(left){
            while(it_min != tto_range_ts_sorted.first){
                tto_ensemble::const_tto_iterator it_left = it_min;
                --it_left;
                if(fabs(it_left->ts - ts_value_left) <= ilen * 0.001){
                    --it_min; ++n;
                    result.first = it_min->ts;
                }
                else break;
            }
        }
        if(right){
            while(it_max != tto_range_ts_sorted.second && fabs(it_max->ts - ts_value_right) <= ilen * 0.001){
                result.second = it_max->ts;
                ++it_max; ++n;
            }
        }
    }
    t_interval_coverage res;
    res.interval = result;
    res.coverage = n * 1.0 / n_total;
    return res;
}



namespace{

// fills result; all ordering values are set to 0.0
void get_input(tto_ensemble & result, DatabaseInput & db, const string & truth_column, const string & ts_column){
   vector<string> colnames;
   colnames.push_back(truth_column);
   colnames.push_back(ts_column);
   std::auto_ptr<DatabaseInput::ResultIterator> it = db.query("products", colnames);
   size_t i = 0;
   while(it->has_data()){
       ++i;
       double truth = it->get_double(0);
       double ts = it->get_double(1);
       result.add(tto(truth, ts));
       ++(*it);
   }
}


void add_ordering_lu(tto_ensemble & inp, const std::string & ordering_rule){
    double o_diff = 1;
    if(ordering_rule == "upper") o_diff = -1;
    set<double> truth_values = inp.get_truth_values();
    tto_ensemble new_ensemble;
    for(set<double>::const_iterator it=truth_values.begin(); it!=truth_values.end(); ++it){
        double o = 0;
        tto_ensemble::tto_range r = inp.get_ttos(*it, tto_ensemble::sorted_by_ts);
        for(tto_ensemble::const_tto_iterator it2=r.first; it2!=r.second; ++it2){
            new_ensemble.add(tto(*it, it2->ts, o));
            o += o_diff;
        }
    }
    swap(inp, new_ensemble);
}

void add_ordering_central(tto_ensemble & inp){
    set<double> truth_values = inp.get_truth_values();
    tto_ensemble new_ensemble;
    for(set<double>::const_iterator it=truth_values.begin(); it!=truth_values.end(); ++it){
        const double truth = *it;
        double o = 0;
        tto_ensemble::tto_range r = inp.get_ttos(truth, tto_ensemble::sorted_by_ts);
        tto_ensemble::const_tto_iterator it_left = r.first;
        advance(it_left, distance(r.first, r.second) / 2);
        tto_ensemble::const_tto_iterator it_right = it_left;
        while(it_left!=r.first || it_right!=r.second){
            if(it_left!=r.first){
                --it_left;
                new_ensemble.add(tto(truth, it_left->ts, o++));
            }
            if(it_right!=r.second){
                new_ensemble.add(tto(truth, it_right->ts, o++));
                ++it_right;
            }
        }
    }
    swap(inp, new_ensemble);
}


void add_ordering_central_shortest(tto_ensemble & inp){
    set<double> truth_values = inp.get_truth_values();
    tto_ensemble new_ensemble;
    for(set<double>::const_iterator it=truth_values.begin(); it!=truth_values.end(); ++it){
        const double truth = *it;
        double o = 0;
        tto_ensemble::tto_range r = inp.get_ttos(truth, tto_ensemble::sorted_by_ts);
        tto_ensemble::const_tto_iterator it_left = r.first;
        advance(it_left, distance(r.first, r.second) / 2);
        tto_ensemble::const_tto_iterator it_right = it_left;
        const double ts0 = it_left->ts;
        while(it_left!=r.first || it_right!=r.second){
           if(it_left!=r.first && it_right!=r.second){
                tto_ensemble::const_tto_iterator it_left_next = it_left;
                --it_left_next;
                double next_ts_left = it_left_next->ts;
                double next_ts_right = it_right->ts;
                if(fabs(next_ts_left - ts0) < fabs(next_ts_right - ts0)){
                    --it_left;
                    new_ensemble.add(tto(truth, it_left->ts, o++));
                }
                else{
                     new_ensemble.add(tto(truth, it_right->ts, o++));
                    ++it_right;
                }
            }
            else if(it_left!=r.first){
                --it_left;
                new_ensemble.add(tto(truth, it_left->ts, o++));
            }
            else {
                new_ensemble.add(tto(truth, it_right->ts, o++));
                ++it_right;
            }
        }
    }
    swap(inp, new_ensemble);
}



}


void neyman_belt::add_ordering_fclike(tto_ensemble & ttos){
    const Distribution & dist = model->get_parameter_distribution();
    Ranges ranges(dist);
    ParValues mode;
    dist.mode(mode);

    //get the set of true values of the poi:
    set<double> truth_values = ttos.get_truth_values();
    map<double, map<double, double> > truth_to__ts_to_ordering;
    Data data;
    double stepsize = (truth_range.second - truth_range.first) / (truth_n - 1);
    tto_ensemble ensemble_for_interpolation;
    for(size_t i=0; i<truth_n; ++i){
        //use "i" to generate asimov data, calculate the ts and nll_best:
        double poi_value1 = truth_range.first + i * stepsize;
        mode.set(*poi, poi_value1);
        model->get_prediction(data, mode);
        std::auto_ptr<NLLikelihood> nll = model->get_nllikelihood(data);
        double nll_best = (*nll)(mode);
        if(!isfinite(nll_best)) continue;
        bool exception = false;
        try{
            ts_producer->produce(data, *model);
        }
        catch(Exception & ex){
            exception = true;
        }
        if(exception) continue;
        double ts = save_double_products->get_product(ts_name);
        if(!isfinite(ts)) continue;
        //it2 is the "truth" for which to calculate the lr for ordering:
        for(set<double>::const_iterator it2=truth_values.begin(); it2!=truth_values.end(); ++it2){
            double truth = *it2;
            mode.set(*poi, truth);
            double nll_2 = (*nll)(mode);
            if(!isfinite(nll_2))continue;
            double ordering = nll_2 - nll_best;
            ensemble_for_interpolation.add(tto(truth, ts, ordering));
        }
    }
    if(ensemble_for_interpolation.size() < truth_n * truth_values.size() / 2){
        throw invalid_argument("add_ordering_fclike: too many failures in ts calculation (>50%)");
    }
    //now go through all truth values again and use truth_to__ts_to_ordering
    // to interpolate / extrapolate ...
    tto_ensemble new_ensemble;
    new_ensemble.reserve(ttos.size());
    for(set<double>::const_iterator truth_it=truth_values.begin(); truth_it!=truth_values.end(); ++truth_it){
        const double truth = *truth_it;
        tto_ensemble::tto_range r_interpolation = ensemble_for_interpolation.get_ttos(truth, tto_ensemble::sorted_by_ts);
        tto_ensemble::tto_range tto_range = ttos.get_ttos(truth, tto_ensemble::sorted_by_ts);
        theta_assert(distance(tto_range.first, tto_range.second) >= 2);
        tto_ensemble::const_tto_iterator it1 = r_interpolation.first;
        tto_ensemble::const_tto_iterator it2 = r_interpolation.first;
        ++it2;
        
        for(tto_ensemble::const_tto_iterator i=tto_range.first; i!=tto_range.second; ++i){
            const double ts = i->ts;
            //find ordering value from ts_to_ordering by linear interpolation:
            while(ts > it2->ts && it2!=r_interpolation.second){
                it1 = it2;
                ++it2;
            }
            const double ts_low = it1->ts;
            const double ts_high = it2->ts;
            const double ordering_low = it1->ordering;
            const double ordering_high = it2->ordering;
            double ordering = ordering_low + (ts - ts_low) * (ordering_high - ordering_low) / (ts_high - ts_low);
            new_ensemble.add(tto(truth, ts, ordering));
        }
        progress();
    }
    new_ensemble.swap(ttos);
}


neyman_belt::neyman_belt(const Configuration & cfg): force_increasing_belt(false), progress_errors(0){
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    progress_done = 0;
    size_t n = cfg.setting["cls"].size();
    if(n==0) throw ConfigurationException("'cls' must be a non-empty list");
    for(size_t i=0; i<n; ++i){
        cls.push_back(cfg.setting["cls"][i]);
    }
    ordering_rule = static_cast<string>(cfg.setting["ordering_rule"]);
    if(ordering_rule!="lower" && ordering_rule != "upper" && ordering_rule != "central"
         && ordering_rule != "central_shortest" && ordering_rule != "fclike"){
          throw ConfigurationException("unknown ordering_rule '" + ordering_rule + "'");
    }
    if(cfg.setting.exists("force_increasing_belt")){
        force_increasing_belt = cfg.setting["force_increasing_belt"];
    }
    
    toy_database = PluginManager<DatabaseInput>::build(Configuration(cfg, cfg.setting["toy_database"]));
    truth_column = static_cast<string>(cfg.setting["truth_column"]);
    ts_column = static_cast<string>(cfg.setting["ts_column"]);
    
    output_database = PluginManager<Database>::build(Configuration(cfg, cfg.setting["output_database"]));
    if(ordering_rule == "fclike"){
        save_double_products.reset(new SaveDoubleProducts());
        cfg.pm->set<ProductsSink>("default", save_double_products);
        model = PluginManager<Model>::build(Configuration(cfg, cfg.setting["fclike_options"]["model"]));
        poi.reset(new ParId(vm->get_par_id(cfg.setting["fclike_options"]["truth"])));
        ts_producer = PluginManager<Producer>::build(Configuration(cfg, cfg.setting["fclike_options"]["ts_producer"]));
        //ts_name is the "module only" part of ts_column:
        size_t p = ts_column.find("__");
        if(p==string::npos){
            throw ConfigurationException("ts_column must contain '__'");
        }
        ts_name = ts_column.substr(p+2);
        if(cfg.setting["fclike_options"].exists("truth_n")){
            truth_n = static_cast<unsigned int>(cfg.setting["fclike_options"]["truth_n"]);
        }
        else{
            truth_n = 200;
        }
        if(cfg.setting["fclike_options"].exists("truth_range")){
            truth_range.first = cfg.setting["fclike_options"]["truth_range"][0];
            truth_range.second = cfg.setting["fclike_options"]["truth_range"][1];
        }
        else{
            truth_range.first = truth_range.second = NAN;
        }
    }
}

void neyman_belt::run(){
    tto_ensemble ttos;
    get_input(ttos, *toy_database, truth_column, ts_column);
    set<double> truth_values = ttos.get_truth_values();
    truth_range.first = *truth_values.begin();
    truth_range.second = *truth_values.rbegin();
    progress_total = truth_values.size() * cls.size();
    if(ordering_rule=="lower" || ordering_rule=="upper"){
        add_ordering_lu(ttos, ordering_rule);
    }
    else if(ordering_rule=="central"){
        add_ordering_central(ttos);
    }
    else if(ordering_rule=="central_shortest"){
        add_ordering_central_shortest(ttos);
    }
    else{
       theta_assert(ordering_rule=="fclike" && model.get()!=0 && poi.get()!=0 && ts_producer!=0);
       progress_total += truth_n;
       add_ordering_fclike(ttos);
       truth_values = ttos.get_truth_values();
    }
    
    //find intervals:
    std::auto_ptr<Table> belt_table = output_database->create_table("belt");
    Column c_truth = belt_table->add_column("truth", typeDouble);
    std::vector<Column> c_lower, c_upper, c_coverage;
    for(size_t i=0; i<cls.size(); ++i){
        stringstream suffix;
        suffix << setw(5) << setfill('0') << static_cast<int>(cls[i] * 10000 + 0.5);
        c_lower.push_back(belt_table->add_column("ts_lower" + suffix.str(), typeDouble));
        c_upper.push_back(belt_table->add_column("ts_upper" + suffix.str(), typeDouble));
        c_coverage.push_back(belt_table->add_column("coverage_upper" + suffix.str(), typeDouble));
    }
    const double inf = numeric_limits<double>::infinity();
    for(size_t i=0; i<cls.size(); ++i){
        pair<double, double> previous_interval;
        previous_interval.first = previous_interval.second = -inf;
        const double cl = cls[i];
        for(set<double>::const_iterator truth_it=truth_values.begin(); truth_it!=truth_values.end(); ++truth_it){
            //relax the contraints of increasing intervals by one permille to circumvent rounding issues
            // arising from discrete ts values:
            const double truth = *truth_it;
            double l = 0.0;
            if(!std::isinf(previous_interval.first) && !std::isinf(previous_interval.second)){
                l = 0.001 * (previous_interval.second - previous_interval.first);
            }
            t_interval_coverage interval = construct_interval(ttos.get_ttos(truth, tto_ensemble::sorted_by_ordering),
                ttos.get_ttos(truth, tto_ensemble::sorted_by_ts), cl, previous_interval.first - l, previous_interval.second - l);
            if(ordering_rule=="lower"){
                interval.interval.second = inf;
            }
            else if(ordering_rule=="upper"){
                interval.interval.first = -inf;
            }
            Row row;
            row.set_column(c_truth, truth);
            row.set_column(c_lower[i], interval.interval.first);
            row.set_column(c_upper[i], interval.interval.second);
            row.set_column(c_coverage[i], interval.coverage);
            belt_table->add_row(row);
            if(force_increasing_belt){
               previous_interval.first = max(interval.interval.first, previous_interval.first);
               previous_interval.second = max(interval.interval.second, previous_interval.second);
            }
           //if force_increasing_belt is false, previous_interval is never updated and stays at -infinity, i.e., no restriction.
           progress();
        }
    }
}

REGISTER_PLUGIN(neyman_belt)

