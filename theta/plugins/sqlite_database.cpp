#include "plugins/sqlite_database.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"
#include "interface/redirect_stdio.hpp"

#include <sstream>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/scoped_array.hpp>

using namespace std;
using namespace theta;

sqlite_database::sqlite_database(const Configuration & cfg) :
    db(0), transaction_active(false), save_all_products(true){
    std::string filename = cfg.setting["filename"];
    if(cfg.setting.exists("products_data")){
        if(cfg.setting["products_data"].get_type()==Setting::TypeString){
            string s = cfg.setting["products_data"];
            if(s=="*")save_all_products = true;
            else throw ConfigurationException("products_data setting is a string but not '*'");
        }
        else{
            save_all_products = false;
            size_t n = cfg.setting["products_data"].size();
            for(size_t i=0; i<n; ++i){
                string column_name = cfg.setting["products_data"][i];
                products_data.insert(column_name);
                if(column_name=="*"){
                    save_all_products = true;
                    products_data.clear();
                    break;
                }
            }
            //if anything is written at all, also write runid and eventid:
            if(products_data.size()){
                products_data.insert("runid");
                products_data.insert("eventid");
            }
        }
    }
    if (boost::filesystem::exists(filename)) {
        boost::filesystem::remove(filename);
    }
    int res = sqlite3_open(filename.c_str(), &db);
    if (res != SQLITE_OK) {
        stringstream ss;
        ss << "sqlite_database constructor failed (filename=\"" << filename << "\") with SQLITE code " << res;
        error(ss.str());//throws.
    }
    maxcol_per_table = sqlite3_limit(db, SQLITE_LIMIT_COLUMN, -1); // default: 2000
    maxcol_per_table = min(maxcol_per_table, sqlite3_limit(db, SQLITE_LIMIT_VARIABLE_NUMBER, -1)); // default: 999
    beginTransaction();
}

void sqlite_database::close() {
    if (!db)
        return;
    endTransaction();
    //finalize all statements associated with the connection:
    sqlite3_stmt *pStmt;
    while ((pStmt = sqlite3_next_stmt(db, 0)) != 0) {
        sqlite3_finalize(pStmt);
    }
    int res = sqlite3_close(db);
    //set db to zero even in case of failure as it is unlikely
    // that the error is recoverable ...
    db = 0;
    if (res != 0) {
        error("sqlite_database::close()");
    }
}

sqlite_database::~sqlite_database() {
    //Close, but do not throw on failure, just print it:
    try {
        close();
    } catch (Exception & e) {
        theta::err << "Exception in sqlite_database destructor: " << e.message << endl << "Ingoring." << endl;
    }
}

void sqlite_database::beginTransaction() {
    if(!transaction_active){
        exec("BEGIN;");
        transaction_active = true;
    }
}

void sqlite_database::endTransaction() {
    if (transaction_active)
        exec("END;");
    transaction_active = false;
}

void sqlite_database::exec(const string & query) {
    char * err = 0;
    if(db==0){
        throw DatabaseException("called sqlite_database::exec although database is closed");
    }
    sqlite3_exec(db, query.c_str(), 0, 0, &err);
    if (err != 0) {
        stringstream ss;
        ss << "sqlite_database::exec(\"" << query << "\") returned error: " << err;
        sqlite3_free(err);
        //database errors should not happen at all. So throw an exception which usually propagates all the way to main.
        throw DatabaseException(ss.str());
    }
}

sqlite3_stmt* sqlite_database::prepare(const string & sql) {
    sqlite3_stmt * statement = 0;
    int ret = sqlite3_prepare_v2(db, sql.c_str(), sql.size() + 1, &statement, 0);
    if (ret != SQLITE_OK) {
        if (statement != 0) {
            sqlite3_finalize(statement);
        }
        try {
            error( __FUNCTION__);
        } catch (DatabaseException & ex) {
            stringstream ss;
            ss << ex.message << " (return code " << ret << "); SQL was " << sql;
            throw DatabaseException(ss.str());
        }
    }
    return statement;
}

void sqlite_database::error(const string & functionName) {
    stringstream ss;
    ss << "Error in function " << functionName << ": sqlite said '" << sqlite3_errmsg(db) << "'";
    throw DatabaseException(ss.str());
}

std::auto_ptr<Table> sqlite_database::create_table(const string & table_name){
    check_name(table_name);
    sqlite_table * result = new sqlite_table(table_name, boost::dynamic_pointer_cast<sqlite_database>(shared_from_this()));
    if(table_name == "products"){
        result->save_all_columns = save_all_products;
        result->save_columns = products_data;
    }
    return std::auto_ptr<Table>(result);
}


sqlite_database::sqlite_table::sqlite_table(const string & name_, const boost::shared_ptr<sqlite_database> & db_) : Table(db_),
    name(name_), table_created(false), next_insert_index(1), db(db_), save_all_columns(true) {
}

Column sqlite_database::sqlite_table::add_column(const std::string & name, const data_type & type){
    if(table_created) throw invalid_argument("sqlite_table::add_column called after table already created (via call to add_row).");
    if(!save_all_columns && save_columns.find(name) == save_columns.end()) return Column(-1);
    Column result(next_insert_index++);
    column_infos[result] = column_info(name, type);
    return result;
}

void sqlite_database::sqlite_table::create_table(){
    map<Column, column_info>::const_iterator column_it = column_infos.begin();
    int n_tables = column_infos.size() / db->maxcol_per_table;
    if(column_infos.size() % db->maxcol_per_table > 0) ++n_tables;
    for(int itable=0; itable < n_tables; ++itable){
        stringstream ss_create_table, ss_insert;
        ss_create_table << "CREATE TABLE '" << name;
        ss_insert << "INSERT INTO '" << name;
        if(itable > 0){
            ss_create_table << "__" << (itable - 1);
            ss_insert << "__" << (itable - 1);
        }
        ss_create_table << "' (";
        ss_insert << "' (";
        int icol = 0;
        for(; icol < db->maxcol_per_table && column_it != column_infos.end(); ++icol, ++column_it){
            if(icol > 0){
                ss_create_table << ", ";
                ss_insert << ", ";
            }
            ss_create_table << "'" << column_it->second.name << "' ";
            ss_insert << "'" << column_it->second.name << "' ";
            switch(column_it->second.type){
                case typeDouble: ss_create_table << "DOUBLE"; break;
                case typeInt: ss_create_table << "INTEGER(4)"; break;
                case typeString: ss_create_table << "TEXT"; break;
                case typeHisto: ss_create_table << "BLOB"; break;
                default:  throw invalid_argument("sqlite_database::create_table: invalid type parameter given.");
            }
        }
        theta_assert(icol > 0);
        ss_create_table << ");";
        ss_insert << ") VALUES (?";
        for(int i=1; i<icol; ++i){
            ss_insert << ",?";
        }
        ss_insert << ");";
        db->exec(ss_create_table.str());
        insert_statements.push_back(db->prepare(ss_insert.str()));
    }
    theta_assert(column_it == column_infos.end());
    table_created = true;    
}

//create the table if it is empty to ensure that all tables have been created
// even if there are no entries
sqlite_database::sqlite_table::~sqlite_table(){
    // we do not want exceptions to be generated in the destructor, so
    // silently ignore an error here because an error here only means that an error already happened
    // somewhere else. If not, all we loose is an empty table anyway ...
    try {
        if(not table_created){
            create_table();
        }
    }
    catch(...){}
}


void sqlite_database::sqlite_table::add_row(const Row & row){
    if(not table_created) create_table();
    const size_t n_tables = insert_statements.size();
    theta_assert(n_tables > 0);
    for(size_t itable=0; itable < n_tables; ++itable){
        map<Column, column_info>::const_iterator column_it = column_infos.begin();
        int icol = 0; // 0-based index of the column in the current table
        try{// try for catching missing data in the row
            for(; icol < db->maxcol_per_table && column_it != column_infos.end(); ++icol, ++column_it){
                if(column_it->second.type == typeDouble){
                    sqlite3_bind_double(insert_statements[itable], icol+1, row.get_column_double(column_it->first));
                }
                else if(column_it->second.type == typeInt){
                    sqlite3_bind_int(insert_statements[itable], icol+1, row.get_column_int(column_it->first));
                }
                else if(column_it->second.type == typeString){
                    const string & s = row.get_column_string(column_it->first);
                    sqlite3_bind_text(insert_statements[itable], icol+1, s.c_str(), s.size(), SQLITE_TRANSIENT);
                }
                else if(column_it->second.type == typeHisto){
                    const Histogram1D & h = row.get_column_histogram(column_it->first);
                    boost::scoped_array<double> blob_data(new double[h.get_nbins()+4]);
                    blob_data[0] = h.get_xmin();
                    blob_data[1] = h.get_xmax();
                    //set underflow to 0.0:
                    blob_data[2] = 0.0;
                    std::copy(h.get_data(), h.get_data() + h.size(), &blob_data[3]);
                    //set overflow to 0.0:
                    blob_data[h.get_nbins()+3] = 0.0;
                    size_t nbytes = sizeof(double) * (h.get_nbins() + 4);
                    sqlite3_bind_blob(insert_statements[itable], icol+1, &blob_data[0], nbytes, SQLITE_TRANSIENT);
                }
            }
        }
        catch(DatabaseException & ex){ 
            stringstream ss;
            ss << "error '" << ex.message << "' in sqlite_table::add_row";
            if(column_it != column_infos.end()){
                ss << " icol = " << icol << " column name: " << column_it->second.name;
            }
            ex.message = ss.str();
            throw ex;
        }
        int res = sqlite3_step(insert_statements[itable]);
        sqlite3_reset(insert_statements[itable]);
        sqlite3_clear_bindings(insert_statements[itable]);
        if (res != 101) {
            db->error(__FUNCTION__);
        }
    }
}

REGISTER_PLUGIN(sqlite_database)
