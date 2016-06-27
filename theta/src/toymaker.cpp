#include "interface/toymaker.hpp"
#include "interface/main.hpp"
#include "interface/phys.hpp"
#include "interface/data.hpp"
#include "interface/model.hpp"

#include "interface/atomic.hpp"

#include <sstream>

using namespace std;
using namespace theta;

ToyMaker::ToyMaker(std::auto_ptr<DataSource> & data_source_, std::auto_ptr<Model> & model_, boost::ptr_vector<Producer> & producers_, const boost::shared_ptr<BufferingProductsSink> & bs_):
    data_source(data_source_), model(model_), bs(bs_){
    producers.transfer(producers.begin(), producers_.begin(), producers_.end(), producers_);
}

void ToyMaker::run(int n_event, atomic_int * toy_count, atomic_int * toy_error_count){
    bs->clear();
    log_messages.clear();
    Data data;
    for (int eventid = 1; eventid <= n_event; eventid++) {
        if(stop_execution)break;
        try{
            data_source->fill(data);
        }
        catch(theta::Exception & ex){
           ex.message += " (in ToyMaker::run while throwing toy data)";
           throw;
        }
        bool error = false;
        for (size_t j = 0; j < producers.size(); j++) {
            try {
                producers[j].produce(data, *model);
            } catch (Exception & ex) {
                error = true;
                std::stringstream ss;
                ss << "Producer '" << producers[j].get_name() << "' failed: " << ex.message << ".";
                log_messages.push_back(LogMessage(eventid, LogTable::error, ss.str()));
                break;
            }
            catch(std::logic_error & f){
                stringstream ss;
                ss << "Producer '" << producers[j].get_name() << "': " << f.what();
                throw logic_error(ss.str());
            }
        }
        if(toy_count) atomic_inc(toy_count);
        if(!error){
            bs->add_row(eventid);
        }
        else if(toy_error_count) atomic_inc(toy_error_count);
    }
}

