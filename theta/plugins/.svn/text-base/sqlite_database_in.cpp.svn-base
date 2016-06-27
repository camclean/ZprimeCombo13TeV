#include "plugins/sqlite_database_in.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"

#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>

using namespace theta;
using namespace std;

namespace{
    void sqlite_error(const string & prefix, sqlite3 * db){
       stringstream error_ss;
       error_ss << prefix << sqlite3_errmsg(db);
       throw DatabaseException(error_ss.str());
    }

    string sqlite_type_to_string(int sqlite_type){
        switch(sqlite_type){
            case SQLITE_FLOAT: return "SQLITE_FLOAT";
            case SQLITE_INTEGER: return "SQLITE_INTEGER";
            case SQLITE_TEXT: return "SQLITE_TEXT";
            case SQLITE_BLOB: return "SQLITE_BLOB";
            default: return "(unknown)";
        }
    }
}


sqlite_database_in::sqlite_database_in(const theta::Configuration & cfg){
   if(cfg.setting.exists("filename") && cfg.setting.exists("filenames")) throw ConfigurationException("both 'filename' and 'filenames' given");
   vector<string> filenames;
   if(cfg.setting.exists("filename")){
       string filename = utils::replace_theta_dir(cfg.setting["filename"]);
       filenames.push_back(filename);
       n_files = 1;
   }
   else{
       n_files = cfg.setting["filenames"].size();
       if(n_files == 0) throw ConfigurationException("'filenames' is empty");
       for(size_t i=0; i<n_files; ++i){
           string filename = utils::replace_theta_dir(cfg.setting["filenames"][i]);
           filenames.push_back(filename);
       }
   }

   for(size_t i=0; i<n_files; ++i){
       if(!boost::filesystem::exists(filenames[i])){
           throw ConfigurationException("file '" + filenames[i] + "' does not exist");
      }
      if(i==0){
         int res = sqlite3_open_v2(filenames[i].c_str(), &db, SQLITE_OPEN_READONLY, 0);
         if(res!=SQLITE_OK){
            throw DatabaseException("could not open file " + filenames[i]);
         }
      }
      else{
          stringstream ss;
          ss << "attach \"" << filenames[i] << "\" as \"file" << i << "\"";
          int res = sqlite3_exec(db, ss.str().c_str(), 0, 0, 0);
          if(res != SQLITE_OK){
              throw DatabaseException("could not attach file " + filenames[i]);
          }
      }
   }
}

sqlite_database_in::~sqlite_database_in(){
    sqlite3_close(db);
}

std::vector<std::string> sqlite_database_in::get_all_tables(){
    const char query[] = "select name from sqlite_master where type='table'";
    sqlite3_stmt * st = 0;
    int ret = sqlite3_prepare_v2(db, query, sizeof(query), &st, 0);
    if(ret!=SQLITE_OK){
        sqlite_error("get_all_tables: ", db); // throws
    }
    boost::shared_ptr<sqlite3_stmt> finalizer(st, sqlite3_finalize);
    vector<string> result;
    while(SQLITE_ROW == (ret = sqlite3_step(st))){
        const char * table_name = reinterpret_cast<const char*>(sqlite3_column_text(st, 0));
        theta_assert(table_name != 0);
        result.push_back(table_name);
    }
    if(ret!=SQLITE_DONE){
        sqlite_error("get_all_tables, after stepping through result: ", db);
    }
    return result;
}


std::vector<std::pair<std::string, theta::data_type> > sqlite_database_in::get_all_columns(const std::string & table_name){
    stringstream q;
    q << "pragma table_info(\"" << table_name << "\")";
    sqlite3_stmt * st;
    int ret = sqlite3_prepare_v2(db, q.str().c_str(), q.str().size() + 1, &st, 0);
    if(ret!=SQLITE_OK){
        sqlite_error("sqlite_database_in::get_all_columns: ", db); // throws
    }
    boost::shared_ptr<sqlite3_stmt> finalizer(st, sqlite3_finalize);
    std::vector<std::pair<std::string, theta::data_type> > result;
    while((ret=sqlite3_step(st))==SQLITE_ROW){
        string colname = reinterpret_cast<const char*>(sqlite3_column_text(st, 1));
        string s_coltype = reinterpret_cast<const char*>(sqlite3_column_text(st, 2));
        theta::data_type type = (data_type)(-1);
        if(s_coltype.find("INTEGER")==0){
            type = typeInt;
        }
        else if(s_coltype.find("TEXT")==0){
            type = typeString;
        }
        else if(s_coltype.find("BLOB")==0){
            type = typeHisto;
        }
        else if(s_coltype.find("DOUBLE")==0){
            type = typeDouble;
        }
        else sqlite_error("sqlite_database_in::get_all_columns: unknown column type " + s_coltype, db);
        result.push_back(make_pair(colname, type));
    }
    if(ret!=SQLITE_DONE){
        sqlite_error("get_all_columns, while stepping through result: ", db);
    }
    return result;
}

std::auto_ptr<DatabaseInput::ResultIterator> sqlite_database_in::query(const std::string & table_name, const std::vector<std::string> & column_names){
    // check that all columns exist:
    {
        std::vector<std::pair<std::string, theta::data_type> > cols = get_all_columns(table_name);
        std::set<std::string> colset;
        for(size_t i=0; i<cols.size(); ++i){
            colset.insert(cols[i].first);
        }
        for(size_t i=0; i<column_names.size(); ++i){
            if(colset.find(column_names[i]) == colset.end()){
                throw DatabaseException("sqlite_database_in: asked for column '" + column_names[i] + "', but it does not exist for table '" + table_name + "'");
            }
        }
        
    }
    
    stringstream q;
    q << "select ";
    for(size_t i=0; i<column_names.size(); ++i){
        if(column_names[i]==""){
           throw DatabaseException("sqlite_database_in::query: got empty column name. This is not allowed.");
        }
        if(i>0) q << ", ";
        q << "\"" << column_names[i] << "\"";
    }
    string sel_cols = q.str();
    q << " from \"" << table_name << "\"";
    
    for(size_t j=1; j<n_files; ++j){
        q << " union all " << sel_cols << " from file" << j << ".\"" << table_name << "\"";
    }
     
    sqlite3_stmt * statement = 0;
    int ret = sqlite3_prepare_v2(db, q.str().c_str(), q.str().size() + 1, &statement, 0);
    if(ret!=0){
       if(statement!=0){
           sqlite3_finalize(statement);
       }
       stringstream error_ss;
       error_ss << "Could not compile SQL statement " << q.str() << "; sqlite said " << sqlite3_errmsg(db);
       throw DatabaseException(error_ss.str());
   }
   return std::auto_ptr<ResultIterator>(new SqliteResultIterator(q.str(), statement, db));
}


void sqlite_database_in::SqliteResultIterator::operator++(){
    int res = sqlite3_step(statement);
    if(res==SQLITE_ROW) has_data_ = true;
    else{
        has_data_ = false;
        if(res!=SQLITE_DONE){
             sqlite_error("Error while stepping through result; sqlite said: ", db);
        }
    }
}


double sqlite_database_in::SqliteResultIterator::get_double(size_t icol){
    if(sqlite3_column_type(statement, icol)!=SQLITE_FLOAT){
        stringstream ss;
        ss << "sqlite_database_in: column type mismatch for column " << icol << " in query '" << sql_statement << "': asked for double (SQLITE_FLOAT), but type is "
           << sqlite_type_to_string(sqlite3_column_type(statement, icol));
        throw DatabaseException(ss.str());
    }
    return sqlite3_column_double(statement, icol);
}

int sqlite_database_in::SqliteResultIterator::get_int(size_t icol){
    if(sqlite3_column_type(statement, icol)!=SQLITE_INTEGER){
        throw DatabaseException("column type mismatch: asked for int (SQLITE_INTEGER) but sqlite type was " + sqlite_type_to_string(sqlite3_column_type(statement, icol)));
    }
    return sqlite3_column_int(statement, icol);
}

std::string sqlite_database_in::SqliteResultIterator::get_string(size_t icol){
    if(sqlite3_column_type(statement, icol)!=SQLITE_TEXT){
        throw DatabaseException("column type mismatch: asked for string (SQLITE_TEXT) but sqlite type was " + sqlite_type_to_string(sqlite3_column_type(statement, icol)));
    }
    const char * res = reinterpret_cast<const char*>(sqlite3_column_text(statement, icol)); // note: res memory is handeled by sqlite
    if(res==0) throw DatabaseException("error in get_string: NULL pointer");//usually means out of memory
    return string(res);
}

theta::Histogram1D sqlite_database_in::SqliteResultIterator::get_histogram(size_t icol){
    if(sqlite3_column_type(statement, icol)!=SQLITE_BLOB){
        throw DatabaseException("column type mismatch: asked for Histogram (SQLITE_BLOB) but sqlite type was " + sqlite_type_to_string(sqlite3_column_type(statement, icol)));
    }
    const double * data = static_cast<const double*>(sqlite3_column_blob(statement, icol));
    //the blob holds nbins + 4 doubles: xmin, xmax, underflow, histgramdata (nbin doubles), overflow
    int nbins = sqlite3_column_bytes(statement, icol) / sizeof(double) - 4;
    if(nbins <= 0) throw DatabaseException("illegal Histogram value: nbins <= 0");
    double xmin = data[0];
    double xmax = data[1];
    if(xmin >= xmax) throw DatabaseException("illegal Histogram value: xmin >= xmax");
    Histogram1D result(nbins, xmin, xmax);
    //data[2] is underflow, data[3] is the data for the first bin in the range
    for(int i=0; i<nbins; ++i){
         result.set(i, data[i+3]);
    }
    return result;
}

REGISTER_PLUGIN(sqlite_database_in)

