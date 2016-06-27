#include "sqlite3.h"

#include <string>
#include <sstream>
#include <iostream>
#include <stdlib.h>

#include "interface/exception.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

using namespace std;
using namespace theta;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

void sqlite3_exec(sqlite3 * db, const char * query){
   char * err = 0;
   sqlite3_exec(db, query, 0, 0, &err);
   if(err!=0){
      stringstream ss;
      ss << err;
      cerr << "sqlite3_exec returned error: " << ss.str() << endl;
      sqlite3_free(err);
      throw Exception(ss.str());
   }
}

sqlite3_stmt* sqlite3_prepare(sqlite3 * db, const char* sql){
    sqlite3_stmt * statement = 0;
    int ret = sqlite3_prepare(db, sql, -1, &statement, 0);
    if(ret!=0){
        if(statement!=0){
            sqlite3_finalize(statement);
        }
        stringstream error_ss;
        error_ss << "Could not compile SQL statement " << sql;
        throw Exception(error_ss.str());
    }
    return statement;
}

void sqlite_error(int code, const string & function){
    cerr << "SQL error " << code << " in " << function << ". Exiting." << endl;
    exit(1);
}

//returns a list of column names in the given table
vector<string> get_column_names(sqlite3 * db, const string & table, const string & database=""){
    vector<string> result;
    stringstream ss;
    if(database!="")
        ss << "PRAGMA " << database << ".";
    else
        ss << "PRAGMA ";
    ss << "table_info('" << table << "');";
    sqlite3_stmt * statement = sqlite3_prepare(db, ss.str().c_str());
    int ret;  
    while((ret=sqlite3_step(statement))==SQLITE_ROW){
        result.push_back(reinterpret_cast<const char*>(sqlite3_column_text(statement, 1)));
    }
    sqlite3_finalize(statement);
    if(ret!=SQLITE_DONE){
        sqlite_error(ret, __FUNCTION__); //exits
    }    
    return result;
}

//returns a list of sorted table names in the given database
vector<string> get_tables(sqlite3 * db, const string & database = ""){
    vector<string> result;
    string master_table = "sqlite_master";
    if(database!="") master_table = database + ".sqlite_master";
    stringstream ss;
    ss << "SELECT name FROM " << master_table << " WHERE type='table' ORDER BY name;";
    sqlite3_stmt * statement = sqlite3_prepare(db, ss.str().c_str());
    int ret;  
    while((ret=sqlite3_step(statement))==SQLITE_ROW){
        result.push_back(reinterpret_cast<const char*>(sqlite3_column_text(statement, 0)));
    }
    sqlite3_finalize(statement);
    if(ret!=SQLITE_DONE){
        sqlite_error(ret, __FUNCTION__); //exits
    }
    return result;
}

int select_int(sqlite3* db, const string & sql){
    sqlite3_stmt * statement = sqlite3_prepare(db, sql.c_str());
    int ret;
    if((ret = sqlite3_step(statement))!=SQLITE_ROW){
        sqlite_error(ret, __FUNCTION__);
    }
    int result = sqlite3_column_int(statement, 0);
    sqlite3_finalize(statement);
    return result;
}

//merge the sqlite databases in file1 and file2 into file1.
//The files must contain the same tables with the same schema.
//Otherwise, the command will fail with an exception.
//Note that the contents of file1is undefined if trying to merge with a file of different database schema:
// some tables might already be merged, others not.
//TODO: could be done better (=faster), if merging more than two
// files at once ...
// errorflag is set to a non-zero value in case of errors
void merge(const string & file1, const string & file2){
    sqlite3 * db=0;
    sqlite3_open(file1.c_str(), &db);
    //this saves some time ...:
    sqlite3_exec(db, "PRAGMA journal_mode=OFF;");
    sqlite3_exec(db, "PRAGMA cache_size=5000;");

    stringstream ss;
    ss << "attach \"" << file2 << "\" as o";
    sqlite3_exec(db, ss.str().c_str());

    //checks:    
    //1. compare tables in the files
    vector<string> tables = get_tables(db);
    vector<string> other_tables = get_tables(db, "o");
    if(tables != other_tables){
        cerr << "Error: tables are not identical in the files to merge ('" << file1 << "', '" << file2 << "'). Skipping." << endl;
        return;
    }
    
    //2. compare random number seed
    int count_rnd = select_int(db, "select count(*) from 'rndinfo' as r, o.'rndinfo' as s where r.seed = s.seed;");
    if(count_rnd!=0){
        cerr << "Error: the random seeds are identical in the files to merge ('" << file1 << "', '" << file2 << "')." << endl;
        throw Exception("random seeds are the same");
        return;
    }
    
    //Ok, we are through with the checks; now compute the runid offset.
    //find out the runids in the first and second table:
    int max_runid_file1, min_runid_file2;
    //use the 'rndinfo' table:
    max_runid_file1 = select_int(db, "select max(runid) from 'rndinfo';");
    min_runid_file2 = select_int(db, "select min(runid) from o.'rndinfo';");
    int offset = max_runid_file1 - min_runid_file2 + 1;

    //Go through the tables and merge the result:
    try{
        sqlite3_exec(db, "BEGIN");
        for(size_t itable=0; itable<tables.size(); itable++){
            vector<string> column_names1 = get_column_names(db, tables[itable]);
            vector<string> column_names2 = get_column_names(db, tables[itable], "o");
            if(column_names1!=column_names2){
                cerr << "Error: columns in table '"<< tables[itable] << "' not identical in the files to merge ('"
                     << file1 << "', '" << file2 << "'): " << endl << "   columns in file 1: ";
                for(size_t i=0; i<column_names1.size(); ++i){
                    cerr << column_names1[i] << " ";
                }
                cerr << endl << "   columns in file 2: ";
                for(size_t i=0; i<column_names2.size(); ++i){
                    cerr << column_names2[i] << " ";
                }
                cerr << endl;
                throw Exception("columns not identical");
            }
            if(column_names1[0]=="runid"){
                ss.str("");
                ss << "INSERT INTO '" << tables[itable] << "' SELECT runid + " << offset;
                //skip runid, which should be the first column in all tables
                for(size_t ic=1; ic<column_names1.size(); ++ic){
                    ss << ", \"" << column_names1[ic] << "\"";
                }
                ss << " FROM o.'" << tables[itable] << "';";
                sqlite3_exec(db, ss.str().c_str());
            }
            else{
                cerr << "Error: table '" << tables[itable] << "' does not contain runid! Exiting." << endl;
                throw Exception("no runid in a table");
            }
        }
        sqlite3_exec(db, "END");
    }
    catch(Exception & e){
        sqlite3_close(db);
        throw;
    }
    sqlite3_close(db);
}


namespace{
    //note: these functions are useful to have compatibility
    // with both V2 and V3 of boost::filesystem
    // as in V2, path::filename returns a string whereas in
    // V3, path::filename returns a path.
    std::string to_string(const std::string & s){
        return s;
    }
    std::string to_string(const fs::path & p){
        return p.string();
    }
}

/** searches path non-recursively for files matching pattern and fills the found files in \c files (prefixed with \c path).
 */
void find_files(const string & path, const string & pattern, vector<string> & files){
  if (!fs::exists(path)){
      stringstream ss;
      ss << "path '" << path << "' does not exist!";
      throw Exception(ss.str());
  }
  fs::directory_iterator end_itr;
  boost::regex pattern_regex(pattern);
  for (fs::directory_iterator itr(path); itr != end_itr; ++itr ) {
      if(boost::regex_search(to_string(itr->path().filename()), pattern_regex)) {
          stringstream ss;
          ss << path;
          if(path[path.size()-1]!='/') ss << "/";
          ss << to_string(itr->path().filename());
          files.push_back(ss.str());
    }
  }
}


int main(int argc, char** argv){
    po::options_description desc("Supported options");
    desc.add_options()("help", "show help message")
    ("outfile", po::value<string>(), "output file of merging")
    ("verbose,v", "verbose output (with progress)")
    ("in-dir", po::value<string>(), "input directory (all files matching *.db there will be merged)");

    po::options_description hidden("Hidden options");

    hidden.add_options()
    ("in-file", po::value< vector<string> >(), "in file");

    po::positional_options_description p;
    p.add("in-file", -1);

    po::options_description cmdline_options;
    cmdline_options.add(desc).add(hidden);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
    po::notify(vm);

    if(vm.count("help")){
        cout << desc << endl;
        return 0;
    }
    if(vm.count("outfile")==0){
        cerr << "please specify an output file with --outfile=..." << endl;
        return 1;
    }

    string outfile = vm["outfile"].as<string>();
    bool verbose = vm.count("verbose");

    string in_dir;
    if(vm.count("in-dir")){
        in_dir = vm["in-dir"].as<string>();
    }

    vector<string> input_files;
    if(vm.count("in-file")){
        input_files = vm["in-file"].as< vector<string> >();
    }

    //go through input dir and add *.db to input_files:
    if(in_dir!=""){
        try{
            find_files(in_dir, "\\.db$", input_files);
        }
        catch(Exception & e){
            cerr << "Error while adding files in input directory: " << e.message << endl;
            return 1;
        }
    }

    if(input_files.size() < 1){
        cerr << "no input files" << endl;
        return 1;
    }

    bool created_output = false;
    for(size_t i=0; i<input_files.size(); i++){
        if(not fs::is_regular_file(input_files[i])){
            cerr << "Input file '" << input_files[i] << "' not found (or not a file). Skipping." << endl;
            continue;
        }
        if(not created_output){
            //create output by copying:
            if(verbose){
                cout << "Copying '" << input_files[i] << "' to '" << outfile << "' ... "<< flush;
            }
            try{
                if(fs::exists(outfile)){
                    fs::remove(outfile);
                }
                fs::copy_file(input_files[i], outfile);
            }
            catch(fs::filesystem_error & ex){
                cout << "error while copying '" << input_files[i] << "' to '" << outfile << "': " << ex.what() << endl;
                return 2;
            }
            if(verbose){
                cout << "done." << endl;
            }
            created_output = true;
        }
        else{
            if(verbose){
                cout << "Merging '" << input_files[i] << "' to '"<< outfile <<"' ... " << flush;
            }
            try{
                merge(outfile, input_files[i]);
            }
            catch(Exception & e){
                cerr << "Error while merging '" << input_files[i] << "': " << e.message << endl;
                cerr << "Deleting output and exiting." << endl;
                fs::remove(outfile);
                exit(1);
            }
            if(verbose){
                cout << "done." << endl;
            }
        }
    }
    return 0;
}


