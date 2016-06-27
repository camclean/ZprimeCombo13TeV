#include "plugins/multiplexing_database.hpp"
#include "interface/plugin.hpp"

using namespace std;
using namespace theta;

multiplexing_database::multiplexing_database(const Configuration & cfg){
   const size_t n = cfg.setting["databases"].size();
   for(size_t i=0; i<n; ++i){
       boost::shared_ptr<Database> db(PluginManager<Database>::build(Configuration(cfg, cfg.setting["databases"][i])).release());
       underlying_dbs.push_back(db);
   }
}

multiplexing_database::~multiplexing_database() {}

std::auto_ptr<Table> multiplexing_database::create_table(const string & table_name){
    check_name(table_name);
    boost::ptr_vector<theta::Table> underlying_tables;
    for(size_t i=0; i<underlying_dbs.size(); ++i){
        underlying_tables.push_back(underlying_dbs[i]->create_table(table_name));
    }
    return std::auto_ptr<Table>(new multiplexing_table(boost::dynamic_pointer_cast<multiplexing_database>(shared_from_this()), underlying_tables));
}

multiplexing_database::multiplexing_table::~multiplexing_table(){}

Column multiplexing_database::multiplexing_table::add_column(const std::string & name, const data_type & type){
    Column result = Column(next_icol++);
    column_types[result] = type;
    for(size_t i=0; i<underlying_tables.size(); ++i){
        external_to_underlying_cols[result].push_back(underlying_tables[i].add_column(name, type));
    }
    return result;
}

void multiplexing_database::multiplexing_table::add_row(const Row & row){
    for(size_t i=0; i<underlying_tables.size(); ++i){
        Row underlying_row;
        // convert columns:
        for(std::map<theta::Column, std::vector<theta::Column> >::const_iterator it=external_to_underlying_cols.begin(); it!=external_to_underlying_cols.end(); ++it){
            const theta::Column & col = it->first;
            const theta::Column & underlying_col = it->second[i];
            switch(column_types[it->first]){
                case typeDouble: underlying_row.set_column(underlying_col, row.get_column_double(col)); break;
                case typeInt: underlying_row.set_column(underlying_col, row.get_column_int(col)); break;
                case typeString: underlying_row.set_column(underlying_col, row.get_column_string(col)); break;
                case typeHisto: underlying_row.set_column(underlying_col, row.get_column_histogram(col)); break;
                default: theta_assert(false);
            }
        }
        underlying_tables[i].add_row(underlying_row);
    }
}

multiplexing_database::multiplexing_table::multiplexing_table(const boost::shared_ptr<Database> & db, boost::ptr_vector<theta::Table> & underlying_tables_): Table(db), next_icol(0){
    underlying_tables.transfer(underlying_tables.end(), underlying_tables_.begin(), underlying_tables_.end(), underlying_tables_);
}

REGISTER_PLUGIN(multiplexing_database)


