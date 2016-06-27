#ifndef PLUGINS_CLS_LIMITS
#define PLUGINS_CLS_LIMITS

#include "interface/decls.hpp"
#include "interface/main.hpp"
#include "interface/phys.hpp"
#include "interface/database.hpp"
#include "interface/data.hpp"
#include "interface/histogram.hpp"

#include <map>
#include <memory>
#include <boost/shared_ptr.hpp>

class SaveDoubleColumn;
class data_filler;
class truth_ts_values;

/*
 *   reuse_toys = {
 *       input_database = {
 *          type = "sqlite_database_in";
 *          filename = "toys0.db";
 *       };
 *       truth_column = "source__beta_signal"; // defaults to "source__truth" which is always cls_limits compatible
 *       poi_column = "lr__poi"; // defaults to (producer name) + "__poi". Only used if current producer depends on poi
 *   };

 * (TODO: this is not implemented yet!) The optional setting \c reuse_toys is a list of theta::DatabaseInput
 * specifications. If present, it is assumed that these contain toys
 * from previous runs with the same model and ts definition (in particular, with the same truth_parameter and ts_column name ...).
 * Reusing toys can speed up the calculation substiantially.
*/



/** \brief Calculate CLs limits
 *
 *
 * \code
 * main = {
 *   type = "cls_limits";
 *
 *   // the first few settings are (almost) the same as usual in type = "default" (theta::Run):
 *   producer = "@lr";
 *   model = "@model";
 *   output_database = {
 *      type = "sqlite_database";
 *      filename = "toys.db";
 *   };
 *   
 *   truth_parameter = "beta_signal";
 *   debuglog = "debug.txt"; // optional, default is no debug output
 *   minimizer = { ... };
 * 
 *   nuisancevalues-for-toys = "datafit"; // optional, default is "prior"
 *
 *   // A. in limit setting mode (default):
 *   mode = "set_limits"; // optional; default is "set_limits"
 *   data_source = { ... }; //optional; for the observed data
 *   ts_column = "lr__nll_diff";
 *   cl = 0.90; // optional; default is 0.95
 *   expected_bands = 500; //optional; default is 1000
 *   beta_signal_expected = 1.0; // optional; default is 0.0
 *   data_source_expected = { ... }; // optional. DataSource specification for expected limits.
 *
 *   reltol_limit = 0.001; // optional; default is 0.05
 *   tol_cls = 0.015; //optional; default; is 0.02
 *   limit_hint = (200.0, 240.0); //optional; default is finding out on its own
 *
 *   // B. for grid generation mode (TODO: not implemented yet(!))
 *   mode = "generate_grid";
 *   truth_range = [0.0, 300.0];
 *   n_truth = 21;
 *   n_sb_toys_per_truth = 2000; // number of toys for s+b hypothesis per truth value
 *   n_b_toys_per_truth = 2000;  // number of b-only toys per truth value
 * };
 * \endcode
 * 
 * For the construction of CLs limits, toys are drawn from the \c model with different values for the signal parameter
 * given in \c truth_parameter (all other parameters are drawn randomly from the model's parameter distribution). For each toy, the \c producer is run and the
 * value from the given \c ts_column is used as test statistic for the construction of CLs limits.
 *
 * \c data_source and \c expected_bands control for which data / test statistic the limits are calculated.
 * \c data_source is a DataSource specification used as the <em>observed</em> data; provide the actual data here (do <em>not</em> use a random
 * data source here, this is not meaningful). \c expected_bands specifies the number of background-only toys --- which set the value
 * of \c beta_signal to \c beta_signal_expected -- to dice to determine the expected limits.
 *
 * The calculation of CLs limits in general requires making a large number of toys at different truth values, to scan for the 
 * desired CLs value 1 - cl. How many toys are drawn and where is done automatically and is driven by \c reltol_limit:
 * this is the relative 1sigma uncertainty tolerance for the limit. \c tol_cls is an absolute accuracy for the CLs value. It is used
 * to control how many toys are done: if the CLs uncertainty from the finite number of toys is below \c tol_cls, no more toys are made. Setting this to a smaller value
 * does not affect the accuracy of the result directly, but it can increase the robustness of the method.
 * \c nuisancevalues-for-toys specifies how the parameter values for the nuisance parameters are chosen to perform toy experiments. In the
 * default setting, "prior", the parameter-distributution from the model is used to sample values for those parameters. The alternative is "datafit"
 * in which case the nuisance parameter values at a given value for the truth parameter are always the same, namely the values from a maximum likelihood fit to data.
 *
 * If \c limit_hint is given, it should be an interval where the limits are probably contained.
 * It does not have to be exact, but wil be used for a starting point to make some toys which can make the
 * convergence considerably faster.
 *
 * The products produced by each toy is used immidiately for the construction of the CLs limits.
 * They are also written to the \c output_database. The tables are the same as for theta::Run. In addition, the table
 * "cls_limits" is created which contains the columns "index", "limit", and "limit_uncertainty".
 * The "index" column specifies on which data the limit was computed: 0 for actual data (from data_source), all other values for
 * background-only toys. In total, the "cls_limits" table should contain \c expected_bands rows and \c expected_rows + 1, if data_source
 * is non-empty. However, the limit computation might not converge for all toys and some limits might be missing. Also note that
 * a limit of infinity is reported for toys with very low CLb value, in order to save computing time.
 * 
 * \c debuglog is the filename of the debug log. If not filename is given, not debug output is written.
 *
 * The indicated progress and errors refer to the number of toys produced. It is hard to tell how many toys are
 * necessary; apart from the obvious things like \c reltol_limit, the number of required toys 
 * also depends on the CLb value: a low CLb value requires much more toys to reach the desired accuracy.
 * Especially when calculating the expected limit band, low CLb values occur naturally.
 * In such a case, it is not uncommon that between 50k to 200k toys are necessary, or even more.
 * For a typical case of calculating the observed limit, i.e., if CLb is around 0.5,
 * 5-10k toys might be sufficient. This is also true in cases where the test statistic is degenerate as
 * in such a case, very low CLb values do not occur.
 *
 * The indicated error counts the number of toy experiments with an exception. This indicates a problem running the producer
 * on the toy data and in most cases means that the minimization the producer performs was not successful. A non-zero value is
 * usually not of concern, as long as it is low (say < 5%).
 */
class cls_limits: public theta::Main{
public:
    cls_limits(const theta::Configuration & cfg);
    virtual void run();
    
private:
    
    void run_single_truth(double truth, bool bkg_only, int n_event);
    void run_set_limits();
    void run_generate_grid();
    
    // run toys at the given truth value (and at truth=0!) until either
    // * the uncertainty on the CLs value is below tol_cls
    //   or
    // * the CLs value is found to be more than 3 sigma away from the target CLs 1-cl.    
    enum e_mode{
        M_ADAPTIVE, M_MORE
    };
    void run_single_truth_adaptive(std::map<double, double> & truth_to_ts, double ts_epsilon, double truth, int idata, e_mode mode);
    void update_truth_to_ts(std::map<double, double> & truth_to_ts, double ts_epsilon, int idata);

    // make truth_to_ts monotonic
    void correct_truth_to_ts(std::map<double, double> & truth_to_ts);

    void read_reuse_toys();

    enum t_mode {
        m_set_limits, m_generate_grid
    };
    
    boost::shared_ptr<theta::VarIdManager> vm;
    boost::shared_ptr<theta::Model> model;
    //note: the tables hold a shared_ptr to this database to ensure proper destruction order
    boost::shared_ptr<theta::Database> db;

    std::auto_ptr<theta::LogTable> logtable;
    boost::shared_ptr<theta::RndInfoTable> rndinfo_table;

    //the producer to be run on the pseudo data which provides the test statistic:
    std::auto_ptr<theta::Producer> producer;
    theta::ParameterDependentProducer * pp_producer; // points to the ParameterDependentProducer in producer.get(), if type matched. 0 otherwise.
    boost::shared_ptr<theta::ProductsTable> products_table;
    boost::shared_ptr<SaveDoubleColumn> sdc;

    std::auto_ptr<data_filler> source; // source for s+b and b-only toys to make limits
    std::auto_ptr<theta::DataSource> data_source_expected; // source for expected limits. Can be 0; in this case, source is used.
    
    theta::ParId truth_parameter;
    int runid;
    int n_toys, n_toy_errors, n_toys_total;
    std::auto_ptr<std::ostream> debug_out;

    t_mode mode;
    boost::shared_ptr<theta::Minimizer> minimizer;

    // A. for mode = "set_limits":
    std::auto_ptr<theta::Table> cls_limits_table;
    theta::Column cls_limits__index, cls_limits__limit, cls_limits__limit_uncertainty;
    
    theta::ParId pid_limit, pid_lambda;
    
    theta::Data current_data;
    theta::Data data_source_data;
    std::auto_ptr<truth_ts_values> tts;
    
    int expected_bands;
    std::auto_ptr<theta::DataSource> data_source; // source for observed data
    
    std::pair<double, double> limit_hint;
    double reltol_limit, tol_cls;
    double clb_cutoff;
    double cl;
    double truth_max; // maximum truth value ever tried in toys. If larger values seem necessary, the toy is declared as outlier. TODO: should get rid of underlying numerical issue instead ...
    double beta_signal_expected;

    std::auto_ptr<theta::DatabaseInput> input_database;
    std::string input_truth_colname, input_poi_colname, input_ts_colname;
    std::vector<double> input_bonly_ts_pool;

    // B. for mode = "generate_grid":
    std::pair<double, double> truth_range;
    int n_truth, n_sb_toys_per_truth, n_b_toys_per_truth;
};


#endif

