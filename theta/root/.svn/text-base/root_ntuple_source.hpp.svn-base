#include "interface/plugin.hpp"
#include "interface/phys.hpp"
#include "interface/random-utils.hpp"

/** \brief A DataSource selecting a random subset of a ROOT ntuple
 *
 * This DataSource returns pseudo-data templates created from a subset
 * of data in a root ntuple. This is useful if studying correlations between
 * different analyses which share events.
 *
 * To study the correlation between analyses in pseudo experiments,
 * the procedure is the following:
 * <ol>
 *   <li>Create a root tree which contains the data of all events you want to use for modelling as float(!).
 *      For each observable o1 and o2, make a TBranch of type float and optionally a TBRanch of type float which contains an event weight.</li>
 *   <li>Configure a model in theta which describes both, o1 and o2.</li>
 *   <li>Use the root tree from step 1 and this DataSource plugin and run the producers you want to study with the model from step 2.</li>
 * </ol>
 *
 * For an example, see
 * <ol>
 *  <li>\c root/create_testtree.cxx which creates rootfile with the necessary structure containing the observable values for both analyses</li>
 *  <li>\c examples/roottreesource.cfg which contain a model configuration of a simple counting experiment with sigbnal only for the observables o1 and o2 which
 *    uses the tree by create_testtree in step 1.</li>
 * </ol>
 * Have a look at the fitted signal content after executing the example config in step 2 (type "products->Draw(mle__s1:mle__s2)" in the root prompt):
 * the measured signal in o1 and o2 will have a correlation, as expected.
 *
 * For more realistic studies, you have to draw many different background processes from trees. While this is possible by filling one large tree with the correct
 * event weight, it is recommended to use one \c root_ntuple_source per background and use the \link add_sources\endlink plugin for main.data_source. This facilitates
 * studies such as increasing one background which can then be done by adjusting the 'factor' setting (see below) in the according root_ntuple_source.
 *
 * Configuration is done via a setting like
 * \code
 *  data_source = {
 *    type = "root_ntuple_source";
 *    name = "source";
 *    filename = "events.root";
 *    observables = {
 *        o1 = {
 *            branchname = "o1";
 *            nbins = 30;
 *            range = [0.0, 30.0];
 *        };
 *        o2 = {
 *             branchname = "my_o2";
 *             nbins = 100;
 *             range = [-10.0, 42.3];
 *        };
 *    };
 *    //optional settings:
 *    treename = "my_events"; // default is "" which uses the top-level TTree in the file
 *    mean_nevents = 147.7;   // default: sum of relweights from the relweight_branch
 *    factor = 0.9;           // default: 1.0
 *    relweight_branchname = "weight"; //default is to use relweight = 1.0 for all events
 *    rnd_gen = {
 *       seed = 123;
 *    };
 *  };
 * \endcode
 *
 * \c type is always "root_ntuple_source" to select this plugin
 *
 * \c name is required for any DataSource; it is not used in this case.
 *
 * \c filename spefifies the path (absolute or relative to the directory from which theta is called) of the root file. 
 *
 * \c observables contains a settings block for each observable the pseudo-data should be created for. For each observable, you can specify
 *   the branchname, the number of bins and the range used to construct the pseudo-data template.
 *
 * \c treename is the name of the TTree in the root file given by \c filename. If it is the empty string, the top-level TTree in the file will be used.
 *  If there is more than one top-level TTree in the file, an exception will be thrown.
 *
 * \c mean_nevents is the mean number of events the DataSource will return. The number of events actually used will be distributed (approximately)
 *   according to a Poisson distribution with this mean. If not given, this is the sum of weights according to relweight_branchname, calculated before
 *   applying 'factor'.
 *
 * \c factor is an additional over factor for the weights, in addition to mean_nevents.
 *
 * \c relweight_branchname is the branch name which contains the relative weight of the event.
 *
 * \c rnd_gen is optional. It specifies details about the random number generator to use. Useful in sitatuations where you want to
 *    reproduce the exact same pseudo dataset in which case you can force the usage of a certain seed.
 *
 * Given the above configuration, the module loads all observable and weight values from all events in the tree into memory,
 * and calculates the total weight W_tot by adding the weights from the relweight branch. If no relweight_branchname is given,
 * 1.0 is used for each weight, hende W_tot is the number of events in the file in this case. If mean_nevents is not specified, it
 * is set to W_tot.
 *
 * To generate a pseudo dataset, each event i  is included with probability  w_i / W_tot * mean_nevents * factor where w_i is the event weight according
 * to the weight branch and 'mean_nevents' and 'factor' as given in the configuration.
 * From these included events, the pseudo data histograms are filled only if the observable values are within the histogram range.
 *
 * The above procedure will not work correctly if the above probability is larger than 1 for any event. Therefore,
 * if any event weight w_i is such that w_i / W_tot * mean_nevents * factor > 1, this is reported as error during the time of construction.
 *
 * For the common case that the relative weights are the same for all events, the number of events follow a binomial distribution around mean_nevents. This
 * agrees with a Poisson around mean_nevents in the limit of infinite (unweights) events.
 */
class root_ntuple_source: public theta::DataSource, public theta::RandomConsumer {
public:
    /// \brief Constructor used by the plugin system
    root_ntuple_source(const theta::Configuration & cfg);
    
    virtual void fill(theta::Data & data_out);
    
    virtual std::auto_ptr<theta::DataSource> clone(const theta::PropertyMap & pm) const{
        throw std::invalid_argument("clone not implemented for root_ntuple_source");
    }
    
private:
   std::vector<theta::ObsId> observables;
   std::vector<size_t> observable_nbins;
   std::vector<std::pair<double,double> > observable_ranges;
   //data from one file is saved in t_data:
   struct t_data{
       //probabilities[i] contains the probability to select event i, assuming that f_u is 1, i.e.,
       // probablities[i] = w_i / W_tot_nominal * mean_nevents    where w_i is the original event weight from the file,
       // W_tot_nominal is the sum of all event weights in the *nominal* sample and mean_nevents is the number of
       // events to select.
       std::vector<double> probabilities;
       // observable k for event i is saved in   (n_obs * i + k)
       std::vector<double> observable_values;
       double total_weight;
       void apply_factor(double factor){
           for(size_t i=0; i<probabilities.size(); ++i){
               double tmp = (probabilities[i] *= factor);
               if(tmp > 1) throw theta::ConfigurationException("probability > 1");
           }
           total_weight *= factor;
       }
       void fill(const std::string & filename, const std::string & treename,
               const std::string & relweight_branchname, const std::vector<std::string> & observable_branchnames);
   };
   t_data nominal_data;
};

