//creates Histograms from result database.

#include "sqlite3.h"

#include "TH1.h"
#include "TFile.h"

#include <string>
#include <sstream>
#include <iostream>

#include <cmath>

#include "libconfig/libconfig.h++"

using namespace std;
using namespace libconfig;

//these two routines for easier sqlite3 are taken from database/merge.cpp:
void sqlite3_exec(sqlite3 * db, const char * query){
   char * err = 0;
   sqlite3_exec(db, query, 0, 0, &err);
   if(err!=0){
      stringstream ss;
      ss << err;
      cerr << "sqlite3_exec returned error: " << ss.str() << endl;
      sqlite3_free(err);
      throw ss.str();
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
        error_ss << "Could not compile SQL statement " << sql << "; sqlite said " << sqlite3_errmsg(db);
        throw error_ss.str();
    }
    return statement;
}

sqlite3* sqlite3_open(const string & fname){
    sqlite3 * db=0;
    int res = sqlite3_open(fname.c_str(), &db);
    if(res!=SQLITE_OK){
        stringstream ss;
        ss << "sqlite3_open(" << fname << "): " << sqlite3_errmsg(db);
        throw ss.str();
    }
    return db;
}

void create_histo_from_doublequery(const string & infile, const string & query, TDirectory * cd, const string & name, int nbins, double xmin, double xmax){
    cd->cd();
    TH1D* histo = new TH1D(name.c_str(), name.c_str(), nbins, xmin, xmax);
    histo->SetDirectory(cd);
    sqlite3 * db = sqlite3_open(infile);
    sqlite3_stmt* st = sqlite3_prepare(db, query.c_str());
    int res;
    while(SQLITE_ROW == (res=sqlite3_step(st))){
        double result = sqlite3_column_double(st, 0);
        histo->Fill(result);
    }
    if(res!=SQLITE_DONE){
        sqlite3_finalize(st);
        stringstream ss;
        ss << "Error in create_histo: stepping through results returned " << res;
        throw ss.str();
    }
    sqlite3_finalize(st);//ignore error
    sqlite3_close(db);//ignore error
}

void create_histos_from_histoquery(const string & infile, const string & query, TDirectory * cd, const string & name){
    cd->cd();
    sqlite3 * db = sqlite3_open(infile);
    sqlite3_stmt* st = sqlite3_prepare(db, query.c_str());
    int res;
    int counter = 0;
    while(SQLITE_ROW == (res=sqlite3_step(st))){
        const double * data = static_cast<const double*>(sqlite3_column_blob(st, 0));
        
        //the blob holds nbins + 4 doubles: xmin, xmax, underflow, histgramdata (nbin doubles), overflow
        int nbins = sqlite3_column_bytes(st, 0) / sizeof(double) - 4;
        
        double xmin = data[0];
        double xmax = data[1];
        stringstream ss_name;
        ss_name << name << counter++;
        
        //copy the data to the histogram:
        TH1F* histo = new TH1F(ss_name.str().c_str(), ss_name.str().c_str(), nbins, xmin, xmax);
        histo->SetDirectory(cd);
        //histo memory is owned by cd.
        for(int i=0; i<=nbins+1; ++i){
            histo->SetBinContent(i, data[i+2]);
        }
    }
    if(res!=SQLITE_DONE){
        sqlite3_finalize(st);
        stringstream ss;
        ss << "Error in create_histo: stepping through results returned " << res;
        throw ss.str();
    }
    sqlite3_finalize(st);//ignore error
    sqlite3_close(db);//ignore error
}

int main(int argc, char** argv){
    if(argc!=2){
        cerr << "Usage: " << argv[0] << " <config file>" << endl;
        return 1;
    }
    try{
        Config cfg;
        cfg.readFile(argv[1]);
        Setting & root = cfg.getRoot();
        string outfile = root["outfile"];
        TFile file(outfile.c_str(), "recreate");
        int n = root.getLength();
        for(int i=0; i<n; ++i){
            if(root[i].getType() != Setting::TypeGroup) continue;
            string histo_name = root[i].getName();
            string infile = root[i]["file"];
            if(root[i].exists("histo_query")){
                string histo_query = root[i]["histo_query"];
                create_histos_from_histoquery(infile, histo_query, &file, histo_name);
            }
            else{
                int nbins = root[i]["nbins"];
                double xmin = root[i]["range"][0];
                double xmax = root[i]["range"][1];
                string query = root[i]["query"];
                create_histo_from_doublequery(infile, query, &file, histo_name, nbins, xmin, xmax);
            }
        }
        file.Write();
        file.Close();
    }
    catch(SettingException & ex){
        cerr << "Error processing configuration file: " << ex.what() << " at path " << ex.getPath() << endl;
        return 1;
    }
    catch(std::string & ex){
        cerr << "Error: " << endl << ex << endl;
        return 1;
    }
}

