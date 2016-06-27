#include <string>
#include <sstream>

#include <boost/date_time/local_time/local_time.hpp>

#include "interface/database.hpp"
#include "interface/histogram.hpp"
#include "interface/plugin.tcc"

using namespace std;
using namespace theta;

REGISTER_PLUGIN_BASETYPE(theta::Database);
REGISTER_PLUGIN_BASETYPE(theta::DatabaseInput);

DatabaseInput::~DatabaseInput(){}

DatabaseInput::ResultIterator::~ResultIterator(){}

Database::~Database(){}

void Database::check_name(const string & column_name) {
    if (column_name.size() == 0)
        throw DatabaseException("Database::check_name: name was empty.");
    if (not ((column_name[0] >= 'A' && column_name[0] <= 'Z') or (column_name[0] >= 'a' && column_name[0] <= 'z'))) {
        stringstream ss;
        ss << "Database::check_name: '" << column_name << "' does not start with a letter as is required.";
        throw DatabaseException(ss.str());
    }
    for (size_t i = 1; i < column_name.size(); i++) {
        if (column_name[i] >= 'A' && column_name[i] <= 'Z')
            continue;
        if (column_name[i] >= 'a' && column_name[i] <= 'z')
            continue;
        if (column_name[i] >= '0' && column_name[i] <= '9')
            continue;
        if (column_name[i] == '_')
            continue;
        stringstream ss;
        ss << "Database::check_name: name '" << column_name << "' invalid at letter " << i;
        throw DatabaseException(ss.str());
    }
}

Table::~Table(){}

Table::Table(const boost::shared_ptr<Database> & db_): db(db_){}

/* Row */
void Row::set_column(const Column & col, double d){
    doubles[col] = d;
}
void Row::set_column(const Column & col, int i){
    ints[col] = i;
}
void Row::set_column(const Column & col, const std::string & s){
    strings[col] = s;
}
void Row::set_column(const Column & col, const Histogram1D & h){
    histos[col] = h;
}


double Row::get_column_double(const Column & col) const{
    std::map<Column, double>::const_iterator it = doubles.find(col);
    if(it==doubles.end()) throw DatabaseException("Row: column not set");
    return it->second;
}
int Row::get_column_int(const Column & col) const{
    std::map<Column, int>::const_iterator it = ints.find(col);
    if(it==ints.end()) throw DatabaseException("Row: column not set");
    return it->second;
}
const std::string & Row::get_column_string(const Column & col) const{
    std::map<Column, std::string>::const_iterator it = strings.find(col);
    if(it==strings.end()) throw DatabaseException("Row: column not set");
    return it->second;
}
const Histogram1D & Row::get_column_histogram(const Column & col) const{
    std::map<Column, Histogram1D>::const_iterator it = histos.find(col);
    if(it==histos.end()) throw DatabaseException("Row: column not set");
    return it->second;
}

void Row::clear(){
    histos.clear();
    ints.clear();
    doubles.clear();
    strings.clear();
}

/* ProductsTable */
ProductsTable::ProductsTable(std::auto_ptr<Table> & table_): table(table_){
    c_runid = table->add_column("runid", typeInt);
    c_eventid = table->add_column("eventid", typeInt);
}

Column ProductsTable::declare_column_impl(const std::string & colname, const data_type & type){
    return table->add_column(colname, type);
}

void ProductsTable::add_row(int runid, int eventid){
    current_row.set_column(c_runid, runid);
    current_row.set_column(c_eventid, eventid);
    table->add_row(current_row);
}

void ProductsTable::set_product(const Column & c, double d){
    current_row.set_column(c, d);
}

void ProductsTable::set_product(const Column & c, int i){
    current_row.set_column(c, i);
}

void ProductsTable::set_product(const Column & c, const std::string & s){
    current_row.set_column(c, s);
}

void ProductsTable::set_product(const Column & c, const theta::Histogram1D & histo){
    current_row.set_column(c, histo);
}



/* LogTable */
LogTable::LogTable(std::auto_ptr<Table> & table_): level(info), table(table_){
   for(int i=0; i<4; ++i){
       n_messages[i]=0;
   }
   c_runid = table->add_column("runid", typeInt);
   c_eventid = table->add_column("eventid", typeInt);
   c_severity = table->add_column("severity", typeInt);
   c_message = table->add_column("message", typeString);
   c_time = table->add_column("time", typeDouble);
}

const int* LogTable::get_n_messages() const{
    return n_messages;
}
    

void LogTable::set_loglevel(e_severity s){
    level = s;
}

LogTable::e_severity LogTable::get_loglevel() const{
    return level;
}

void LogTable::really_append(int runid, int eventid, e_severity s, const string & message) {
    n_messages[s]++;
    Row row;
    row.set_column(c_runid, runid);
    row.set_column(c_eventid, eventid);
    row.set_column(c_severity, s);
    row.set_column(c_message, message);
    using namespace boost::posix_time;
    using namespace boost::gregorian;
    ptime t(microsec_clock::universal_time());
    time_duration td = t - ptime(date(1970, 1, 1));
    double time = td.total_microseconds() / 1000000.0;
    row.set_column(c_time, time);
    table->add_row(row);
}

//RndInfoTable
RndInfoTable::RndInfoTable(std::auto_ptr<Table> & table_): table(table_){
    c_runid = table->add_column("runid", typeInt);
    c_name = table->add_column("name", typeString);
    c_seed = table->add_column("seed", typeInt);
}

void RndInfoTable::append(int runid, const string & name, int seed){
    Row row;
    row.set_column(c_runid, runid);
    row.set_column(c_name, name);
    row.set_column(c_seed, seed);
    table->add_row(row);
}

