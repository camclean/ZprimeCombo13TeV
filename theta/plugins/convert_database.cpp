#include "plugins/convert_database.hpp"

#include "interface/database.hpp"
#include "interface/plugin.hpp"

using namespace theta;
using namespace std;

void convert_database::run(){
    vector<string> tables = db_in->get_all_tables();
    for(size_t i=0; i<tables.size(); ++i){
        // 1. create table+column structure in output database:
        std::auto_ptr<Table> table_out = db_out->create_table(tables[i]);
        std::vector<std::pair<std::string, theta::data_type> > columns_in = db_in->get_all_columns(tables[i]);
        vector<Column> columns_out(columns_in.size());
        vector<string> colnames;
        const size_t ncols = columns_in.size();
        for(size_t icol=0; icol < columns_in.size(); ++icol){
            colnames.push_back(columns_in[icol].first);
            columns_out[icol] = table_out->add_column(columns_in[icol].first, columns_in[icol].second);
        }
        // 2. transfer data:
        std::auto_ptr<DatabaseInput::ResultIterator> it = db_in->query(tables[i], colnames);
        while(it->has_data()){
            Row row_out;
            for(size_t icol=0; icol < ncols; ++icol){
                switch(columns_in[icol].second){
                    case typeDouble: row_out.set_column(columns_out[icol], it->get_double(icol)); break;
                    case typeInt: row_out.set_column(columns_out[icol], it->get_int(icol)); break;
                    case typeHisto: row_out.set_column(columns_out[icol], it->get_histogram(icol)); break;
                    case typeString: row_out.set_column(columns_out[icol], it->get_string(icol)); break;
                }
            }
            table_out->add_row(row_out);
            ++(*it);
        }
    }
}

convert_database::convert_database(const theta::Configuration & cfg){
    db_out = PluginManager<Database>::build(Configuration(cfg, cfg.setting["output_database"]));
    db_in = PluginManager<DatabaseInput>::build(Configuration(cfg, cfg.setting["input_database"]));
}

convert_database::~convert_database(){}


REGISTER_PLUGIN(convert_database)

