#ifndef TOYMAKER_HPP
#define TOYMAKER_HPP

#include "interface/decls.hpp"
#include "interface/data_type.hpp"
#include "interface/producer.hpp"
#include "interface/database.hpp"

#include <vector>
#include <map>


namespace theta {

// ProductsSink which buffers everything in memory.
class BufferingProductsSink: public ProductsSink{
public:
    BufferingProductsSink(): next_id(0){
        c_eventid = declare_column("eventid", typeInt);
    }

    const std::vector<Row> & get_rows() const{
        return rows;
    }

    void clear(){
        rows.clear();
    }    

    virtual void set_product(const Column & c, double d){
        current_row.set_column(c, d);
    }

    virtual void set_product(const Column & c, int i){
        current_row.set_column(c, i);
    }

    virtual void set_product(const Column & c, const std::string & s){
        current_row.set_column(c, s);
    }   

    virtual void set_product(const Column & c, const Histogram1D & h){
        current_row.set_column(c, h);
    }

    void add_row(int eventid){
        current_row.set_column(c_eventid, eventid);
        rows.push_back(current_row);
        current_row.clear();
    }

    // declares all own columns in ps
    void declare_columns(ProductsSink & ps) const {
        for(std::map<std::string, std::pair<Column, data_type> >::const_iterator it = name_to_column_type.begin(); it!=name_to_column_type.end(); ++it){
            if(it->first=="eventid") continue;
            ps.declare_column(it->first, it->second.second);
        }
    }

    // add_row_cmd(int eventid) must be valid. ps is assumed to have the same columns
    // than this. (This can be ensured by calling declare_columns(ps) ).
    template<typename T>
    void replay_to(ProductsSink & ps, const T & add_row_cmd) const;

private:
    virtual Column declare_column_impl(const std::string & full_column_name, const data_type & type){
        return Column(next_id++);
    }
    int next_id;
    std::vector<Row> rows;
    Row current_row;
    Column c_eventid;
};

// Toy maker is used to separate the (computational intensive) task of making
// toys and running some producers on them from the rest.
// This separation is to enable multi-threading or distributed computing.
//
//
// To use the ToyMaker,
// 1. create a new Configuration instance, where the ProductsSink is set to a BufferingProductsSink, the runid is different from the one
//    before (for uniqueness in the RndInfoTable). The other things like VarIdManager can be the same.
// 2. with the configuration from 1., and the PluginManager, create the instances for DataSource, Model and Producers to run.
// 3. pass the arguments to the constructor of the ToyMaker, which takes the ownership of these objects.
// 4. call ToyMaker::run(); this can be done in a different thread
// 5. the result is in the BufferingProductsSink. To get the log messages, call get_log_messages. Note that the message buffer
//   is cleared at each invocation of run() and the BufferingProductsSink are cleared at the invocation of run, so you should
//   transfer all data you want to keep before the next call.
//
// ToyMaker can be re-used, i.e., execute steps 1., 2., 3. once at the beginning (from the main thread), then 4., 5., repeatedly.
//
// Note that the bufferingProductsSink and the internal log buffer are cleared at each invocation of run.
class ToyMaker{
public:
    struct LogMessage{
        int eventid;
        LogTable::e_severity severity;
        std::string message;
        LogMessage(int eventid_, LogTable::e_severity severity_, const std::string & message_): eventid(eventid_), severity(severity_), message(message_){}
    };

    ToyMaker(std::auto_ptr<DataSource> & data_source, std::auto_ptr<Model> & model, boost::ptr_vector<Producer> & producers, const boost::shared_ptr<BufferingProductsSink> & bs);

    const std::vector<LogMessage> & get_log_messages() const{
        return log_messages;
    }

    // toy_count and toy_error_count can be 0. If they are not 0,
    // they are incremented to indicate progress.
    void run(int n_events, atomic_int * toy_count, atomic_int * toy_error_count);
private:
    std::auto_ptr<DataSource> data_source;
    std::auto_ptr<Model> model;
    boost::ptr_vector<Producer> producers;

    std::vector<LogMessage> log_messages;
    boost::shared_ptr<BufferingProductsSink> bs;
};





template<typename T>
void BufferingProductsSink::replay_to(ProductsSink & ps, const T & add_row_cmd) const{
    const std::map<std::string, std::pair<Column, data_type> > & ps_name_to_cols =  ps.get_name_to_column_type();
    // now, build a map 'our Column' -> 'ps Column', based on the name. Also check that the types match.
    std::map<Column, Column> ours_to_theirs;
    std::map<Column, data_type> ours_to_dt;
    for(std::map<std::string, std::pair<Column, data_type> >::const_iterator it = name_to_column_type.begin(); it!=name_to_column_type.end(); ++it){
        if(it->first == "eventid") continue;
        std::map<std::string, std::pair<Column, data_type> >::const_iterator ps_it = ps_name_to_cols.find(it->first);
        if(ps_it == ps_name_to_cols.end()) throw std::invalid_argument("BufferingProductsSink::replay_to: column '" + it->first + "' not found in target");
        if(it->second.second != ps_it->second.second) throw std::invalid_argument("BufferingProductsSink::replay_to: column type for column '" + it->first + "' does not match");
        ours_to_theirs[it->second.first] = ps_it->second.first;
        ours_to_dt[it->second.first] = it->second.second;
    }
    for(std::vector<Row>::const_iterator it_row = rows.begin(); it_row!=rows.end(); ++it_row){
        const Row & row = *it_row;
        std::map<Column, data_type>::const_iterator dt_it = ours_to_dt.begin();
        for(std::map<Column, Column>::const_iterator c_it=ours_to_theirs.begin(); c_it != ours_to_theirs.end(); ++c_it, ++dt_it){
            theta_assert(dt_it != ours_to_dt.end() and c_it->first == dt_it->first);
            switch(dt_it->second){
                case typeInt: ps.set_product(c_it->second, row.get_column_int(c_it->first)); break;
                case typeDouble: ps.set_product(c_it->second, row.get_column_double(c_it->first)); break;
                case typeHisto: ps.set_product(c_it->second, row.get_column_histogram(c_it->first)); break;
                case typeString: ps.set_product(c_it->second, row.get_column_string(c_it->first)); break;
            }
        }
        int eventid = row.get_column_int(c_eventid);
        add_row_cmd(eventid);
    }
}

}

#endif


