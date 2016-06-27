#include "interface/random-utils.hpp"
#include "interface/random.hpp"
#include "interface/histogram.hpp"
#include "interface/histogram-with-uncertainties.hpp"
#include "interface/database.hpp"
#include "interface/plugin.hpp"

#include <boost/date_time/local_time/local_time.hpp>
#include <unistd.h>

using namespace theta;

RandomConsumer::RandomConsumer(const theta::Configuration & cfg, const std::string & name): seed(-1) {
   std::auto_ptr<RandomSource> rnd_source;
   std::string source_type = "taus";
   if(cfg.setting.exists("rnd_gen")){
       Setting s = cfg.setting["rnd_gen"];
       if(s.exists("source_type")){
           source_type = static_cast<std::string>(s["source_type"]);
       }
       if(s.exists("seed")){
          seed = s["seed"];
       }
   }
   if(source_type=="taus"){
      rnd_source.reset(new RandomSourceTaus());
   }
   else if(source_type == "mt"){
      rnd_source.reset(new RandomSourceMersenneTwister());
   }
   else{
      throw ConfigurationException("unknown source_type given for rnd_gen (valid values are 'taus' and 'mt')");
   }
   if(seed == -1){
       using namespace boost::posix_time;
       using namespace boost::gregorian;
       ptime t(microsec_clock::universal_time());
       time_duration td = t - ptime(date(1970, 1, 1));
       seed = td.total_microseconds();
       // to avoid clashes with other RandomConsumers initialized in the same process in the same
       // clock resolution interval, use a static counter:
       static int rnd_id = 0;
       seed = seed * 33 + rnd_id;
       ++rnd_id;
       // to avoid clashes in case of batch system usage with jobs starting in the same clock resolution
       // interval, also use the hostname for the seed.
       //Note that we should use a length of "HOST_NAME_MAX + 1" here, however, HOST_NAME_MAX is not defined on
       // all POSIX systems, so just use the value HOST_NAME_MAX=255 ...
       char hname[256];
       gethostname(hname, 256);
       // In case the hostname does not fit into hname, the name is truncated but no error is returned.
       // On the other hand, using some random bytes here for seeding is also Ok, just make sure we have a null byte
       // somewhere:
       hname[255] = '\0';
       int c;
       for(size_t i=0; (c = hname[i]); ++i){
           seed = seed * 33 + c;
       }
       // and finally, also include the process id, if theta starts on the same host at the same time:
       seed = seed * 33 + (int)getpid();
   }
   rnd_gen.reset(new Random(rnd_source));
   int runid = *(cfg.pm->get<int>("runid"));
   // note: for runid = 1, do not change the seed.
   // For runid > 1, change the seed by a large value to prevent multithreading jobs from using consecutive seeds.
   // This allows e.g. to have many .cfg files for "background-only" toys with consecutive seeds
   // in a multithreaded setup without seed collissions.
   seed += (runid - 1) * 14657;
   rnd_gen->set_seed(seed);
   if(cfg.pm->exists<RndInfoTable>()){
       cfg.pm->get<RndInfoTable>()->append(runid, name, seed);
   }
}

RandomConsumer::~RandomConsumer(){}

void theta::randomize_poisson(DoubleVector & d, Random & rnd){
    const size_t n = d.size();
    double * data = d.get_data();
    for(size_t i=0; i<n; ++i){
        double mu = data[i];
        if(mu > 0.){
            data[i] = rnd.poisson(mu);
        }
    }
}

Histogram1D theta::randomize_gauss(const Histogram1DWithUncertainties & histo, Random & rnd){
    const size_t n = histo.get_nbins();
    Histogram1D result(n, histo.get_xmin(), histo.get_xmax());
    for(size_t i=0; i<n; ++i){
        const double val = histo.get_value(i);
        const double unc = histo.get_uncertainty(i);
        double new_value;
        do{
            new_value = rnd.gauss(unc) + val;
        }while(val >=0.0 and new_value < 0.0);
        result.set(i, new_value);
    }
    return result;
}


