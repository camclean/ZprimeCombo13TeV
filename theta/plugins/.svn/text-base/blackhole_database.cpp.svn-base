#include "plugins/blackhole_database.hpp"
#include "interface/plugin.hpp"

using namespace std;
using namespace theta;

blackhole_database::blackhole_database(const Configuration & cfg){}

blackhole_database::~blackhole_database() {}

std::auto_ptr<Table> blackhole_database::create_table(const string & table_name){
    check_name(table_name);
    return std::auto_ptr<Table>(new blackhole_table(boost::dynamic_pointer_cast<blackhole_database>(shared_from_this())));
}

blackhole_database::blackhole_table::~blackhole_table(){}

Column blackhole_database::blackhole_table::add_column(const std::string & name, const data_type & type){
    return Column();
}

void blackhole_database::blackhole_table::add_row(const Row & row){}

REGISTER_PLUGIN(blackhole_database)
