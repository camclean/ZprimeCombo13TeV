#include "root/root_ntuple_source.hpp"
#include "interface/random.hpp"
#include "interface/data.hpp"
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TKey.h"

using namespace theta;
using namespace std;

//relweight_branchname can be empty; in this case, 1 is used as weight.
//will fill result.probabilities with weights(!) and result.observable_values with the values
// for the observables.
void root_ntuple_source::t_data::fill(const string & filename, const string & treename,
    const string & relweight_branchname, const vector<string> & observable_branchnames){
    TFile f(filename.c_str(), "read");
    if(f.IsZombie()){
       throw ConfigurationException("could not open root file '" + filename + "'");
    }
    TTree * tree = 0;
    if(treename!=""){
        tree = dynamic_cast<TTree*>(f.Get(treename.c_str()));
    }
    else{
       //list file contents:
       TList * l = f.GetListOfKeys();
       for(int i=0; i<l->GetEntries();++i){
          TKey * ti = dynamic_cast<TKey*>(l->At(i));
          if(ti->GetClassName()!=string("TTree")) continue;
          TTree * tt_tmp = dynamic_cast<TTree*>(f.Get(ti->GetName()));
          assert(tt_tmp != 0);
          if(tree!=0) throw ConfigurationException("more than one TTree found in '" + filename + "' but no treename given!");
          tree = tt_tmp;
       }
    }
    if(tree==0) throw ConfigurationException("TTree '" + treename + "' in file '" + filename + "' not found");
    size_t n_events = tree->GetEntries();
    size_t n_obs = observable_branchnames.size();
    vector<float> tree_in(n_obs);
    
    //TODO: report error if type of branch is not float!!!
    TBranch * b;
    for(size_t i=0; i < n_obs; ++i){
        b = 0;
        tree->SetBranchAddress(observable_branchnames[i].c_str(), &tree_in[i], &b);
        if(b==0) throw ConfigurationException("tree branch '" + observable_branchnames[i] + "' in tree '" + filename + ":" + treename + "' not found");
    }
    float relweight = 1.0;
    if(relweight_branchname!=""){
        b = 0;
        tree->SetBranchAddress(relweight_branchname.c_str(), &relweight, &b);
        if(b==0) throw ConfigurationException("tree branch '" + relweight_branchname + "' in tree '" + filename + ":" + treename+ "' not found");
    }
    probabilities.resize(n_events);
    observable_values.resize(n_events * n_obs);
    total_weight = 0.0;
    for(size_t i=0; i<n_events; ++i){
        tree->GetEntry(i);
        probabilities[i] = relweight;
        total_weight += relweight;
        if(relweight < 0) throw ConfigurationException("got negative weight");
        copy(tree_in.begin(), tree_in.end(), observable_values.begin() + n_obs * i);
    }
    f.Close();
}


root_ntuple_source::root_ntuple_source(const Configuration & cfg): DataSource(cfg), RandomConsumer(cfg, get_name()){
    string filename = utils::replace_theta_dir(cfg.setting["filename"]);
    string treename;
    string relweight_branchname;
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    double cfg_factor = 1.0;
    if(cfg.setting.exists("relweight_branchname")){
       relweight_branchname = static_cast<string>(cfg.setting["relweight_branchname"]);
    }
    if(cfg.setting.exists("treename")){
       treename = static_cast<string>(cfg.setting["treename"]);
    }
    if(cfg.setting.exists("factor")){
       cfg_factor = cfg.setting["factor"];
    }
    double mean_nevents = -1;
    if(cfg.setting.exists("mean_nevents")){
        mean_nevents = cfg.setting["mean_nevents"];
    }
    
    Setting so = cfg.setting["observables"];
    size_t n_obs = so.size();
    vector<string> observable_branchnames;
    for(size_t i=0; i < n_obs; ++i){
        string obs_name = so[i].get_name();
        ObsId oid = vm->get_obs_id(obs_name);
        observables.push_back(oid);
        // set the histogram in data:
        int nbins = so[i]["nbins"];
        double xmin = so[i]["range"][0];
        double xmax = so[i]["range"][1];
        if(xmin >= xmax || nbins <=0) throw ConfigurationException("invalid range / binning for observable " + obs_name);
        observable_ranges.push_back(make_pair(xmin, xmax));
        observable_nbins.push_back(nbins);
        string obs_branchname = so[i]["branchname"];
        observable_branchnames.push_back(obs_branchname);
    }
    nominal_data.fill(filename, treename, relweight_branchname, observable_branchnames);
    // correct the probabilities by multiplying with   mean_nevents / total_weight.
    if(mean_nevents < 0) mean_nevents = nominal_data.total_weight;
    double factor = cfg_factor * mean_nevents / nominal_data.total_weight;
    nominal_data.apply_factor(factor);
}

void root_ntuple_source::fill(Data & data_out){
    data_out.reset();
    size_t n_obs = observables.size();
    for(size_t i=0; i<n_obs; ++i){
        data_out[observables[i]] = Histogram1D(observable_nbins[i], observable_ranges[i].first, observable_ranges[i].second);
    }
    //loop over all nominal events:
    size_t n_events = nominal_data.probabilities.size();
    for(size_t i=0; i<n_events; ++i){
        if(rnd_gen->uniform() >= nominal_data.probabilities[i]) continue;
        for(size_t j=0; j<n_obs; ++j){
            data_out[observables[j]].fill(nominal_data.observable_values[n_obs * i + j], 1.0);
        }
    }
}

REGISTER_PLUGIN(root_ntuple_source)
