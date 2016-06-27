#include "root/rootfile_database.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"

#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/scoped_array.hpp>

using namespace std;
using namespace theta;

rootfile_database::rootfile_database(const Configuration & cfg): file(0){
   std::string filename = utils::replace_theta_dir(cfg.setting["filename"]);
   file = new TFile(filename.c_str(), "recreate");
   if(not file->IsOpen()){
       delete file;
       throw ConfigurationException("error opening output root file '" + filename + "'");
   }
   
   if(cfg.setting.exists("products_data")){
      try{
         string s = cfg.setting["products_data"];
         if(s=="*")save_all_products = true;
         else throw ConfigurationException("products_data setting is a string but not '*'");
      }
      catch(libconfig::SettingTypeException & e){
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
   else{
      save_all_products = true;
   }
   
   if(cfg.setting.exists("products_histograms")){
      size_t n = cfg.setting["products_histograms"].size();
      for(size_t i=0; i<n; ++i){
          int nbins = cfg.setting["products_histograms"][i]["nbins"];
          double xmin = cfg.setting["products_histograms"][i]["range"][0];
          double xmax = cfg.setting["products_histograms"][i]["range"][1];
          string name = cfg.setting["products_histograms"][i]["name"];
          string colname = cfg.setting["products_histograms"][i]["column"];
          hist_infos.push_back(hist_info(Histogram1D(nbins, xmin, xmax), name, colname));
      }
   }
}

rootfile_database::~rootfile_database() {
    if(file){
       //write those root histos:
       if(hist_infos.size()){
            TDirectory * dir = file->mkdir("products_histograms");
            for(size_t i=0; i<hist_infos.size(); ++i){
                size_t nbins;
                TH1D * root_histo = new TH1D(hist_infos[i].name.c_str(), hist_infos[i].name.c_str(), nbins = hist_infos[i].h.get_nbins(),
                                                hist_infos[i].h.get_xmin(), hist_infos[i].h.get_xmax());
                for(size_t ibin=0; ibin<nbins; ++ibin){
                    root_histo->SetBinContent(ibin+1, hist_infos[i].h.get(ibin));
                }
                root_histo->SetDirectory(dir);
            }
       }
       file->Write();
       file->Close();
       delete file;
       file = 0;
    }
}

std::auto_ptr<Table> rootfile_database::create_table(const string & table_name){
    check_name(table_name);
    rootfile_table * result = new rootfile_table(table_name, boost::dynamic_pointer_cast<rootfile_database>(shared_from_this()));
    if(table_name == "products"){
        result->save_all_columns = save_all_products;
        result->save_columns = products_data;
    }
    return std::auto_ptr<Table>(result);
}


rootfile_database::rootfile_table::rootfile_table(const std::string & tablename, const boost::shared_ptr<rootfile_database> & db_):Table(db_),
   db(db_), save_all_columns(true), next_id(0){
       tree = new TTree(tablename.c_str(), tablename.c_str());
       tree->SetDirectory(db->file);
       products_table = tablename == "products";
}

rootfile_database::rootfile_table::~rootfile_table(){}

Column rootfile_database::rootfile_table::add_column(const std::string & name, const data_type & type){
    Column result(next_id++);
    column_datas[result].name = name;
    column_datas[result].type = type;
    bool make_branch = save_all_columns || save_columns.find(name)!=save_columns.end();
    switch(type){
        case theta::typeDouble:
            if(make_branch)
                tree->Branch(name.c_str(), &column_datas[result].data.d, "data/D");
            break;
        case theta::typeInt:
            if(make_branch)
                tree->Branch(name.c_str(), &column_datas[result].data.i, "data/I");
            break;
        case theta::typeString:
            column_datas[result].data.s = new TString();
            if(make_branch)
                tree->Branch(name.c_str(), "TString", &column_datas[result].data.s);
            break;
        case theta::typeHisto:
            column_datas[result].data.h = new TH1D(name.c_str(), name.c_str(), 1, 0, 1);
            if(make_branch)
                tree->Branch(name.c_str(), "TH1D", &column_datas[result].data.h);
            break;
    }
    return result;
}


void rootfile_database::rootfile_table::add_row(const Row & row){
    for(map<Column, column_data>::iterator it=column_datas.begin(); it!=column_datas.end(); ++it){
        bool is_histo = false;
        switch(it->second.type){
            case typeDouble:
                it->second.data.d = row.get_column_double(it->first);
                if(products_table){
                    for(size_t i=0; i< db->hist_infos.size(); ++i){
                        if(db->hist_infos[i].column_name == it->second.name){
                            db->hist_infos[i].h.fill(it->second.data.d, 1.0);
                        }
                    }
                }
                break;
            case typeInt:
                it->second.data.i = row.get_column_int(it->first);
                break;
            case typeString:
                *(it->second.data.s) = row.get_column_string(it->first).c_str();
                break;
            case typeHisto:
                is_histo = true;
        }
        if(is_histo){
            const Histogram1D & h = row.get_column_histogram(it->first);
            it->second.data.h->SetBins(h.get_nbins(), h.get_xmin(), h.get_xmax());
            for(size_t i=0; i<h.get_nbins(); ++i){
                it->second.data.h->SetBinContent(i+1, h.get(i));
            }
        }
    }
    tree->Fill();
}

REGISTER_PLUGIN(rootfile_database)
