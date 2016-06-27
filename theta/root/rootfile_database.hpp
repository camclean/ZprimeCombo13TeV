#ifndef PLUGIN_ROOTFILE_DATABASE_HPP
#define PLUGIN_ROOTFILE_DATABASE_HPP

#include "interface/decls.hpp"
#include "interface/database.hpp"
#include "interface/histogram.hpp"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include <string>
#include <set>

/** \brief Database which makes histograms of values instead of saving them all
 *
 * This is especially useful if producing large-scale test statistic which would take
 * up a lot of space if saved.
 *
 * Configured via a setting group like
 * \code
 * database = {
 *   type = "rootfile_database";
 *   filename = "root_out.root";
 *   products_histograms = ({  //optional. Default is to save no histograms.
 *       name = "d_nll";
 *       nbins = 100;
 *       range = [0.0, 100.0];
 *       column = "deltanll__nll_diff";
 *   }, {
 *      name = "nll_sb";
 *      nbins = 1000;
 *      range = [-8000.0, 8000.0];
 *      column = "deltanll__nll_sb";
 *   });
 *   products_data = ("deltanll__nll_sb"); //optional. Default is "*", i.e., to save all columns
 * };
 * \endcode
 *
 * \c type must always be "rootfile_database" in order to select this plugin
 *
 * \c filename is the name of the root file to be created by the plugin. It will be overwritten
 *  if it already exists.
 *
 * The settings \c product_data and \c products_histograms contains settings specifically for the \c products
 * database:
 *
 * \c product_histograms is a list of setting blocks. Each setting block contains:
 * <ul>
 *   <li>\c name the name of the histogram in the root file</li>
 *   <li>\c nbins the number of bins for the histogram</li>
 *   <li>\c range the range for the histogram</li>
 *   <li>\c column the name of the column to fill into the histogram. It must be a column of type double,
 *     otherwise, it will not be filled. As in other databases, this is the name of the producer,
 *     followed by "__" as separator, concatenated with the producer-defined column name.</li>
 * </ul>
 *  Note that overflow and underflow bin are filled for the histograms.
 *
 * \c product_data is a list of column names of the products table. Those columns will be saved as brnaches in a tree.
 *   It can be an empty list "()". It is also possible to use the special string "*" instead of the list (or as part of the list)
 *   to save <em>all</em> columns.
 *   Note that if you specify any column to write in this list, the "usual" columns runid and eventid will be written as well,
 *   even if you do not specify them.
 *
 * The root file contains a directory named "products_histograms" which contains the configured histograms with the columns
 * from the products table.
 *
 * Apart from these histograms, the root file contains one tree per table with one branch for each column in this table.
 * Note that the "products" table is treated specially and will only contain those columns configured in the \c products_data
 * setting.
 *
 * The mapping from theta column types and root/C++ types is straight-forward:
 * <ul>
 *  <li>Table::typeDouble and Table::typeInt are \c double and \c int, respectively</li>
 *  <li>Table::typeString is mapped to \c TString </li>
 *  <li>Table::typeHisto is mapped to \c TH1D </li>
 * </ul>
 */
class rootfile_database: public theta::Database {
public:
    
    /** \brief Constructor for the plugin system
     */
    rootfile_database(const theta::Configuration & cfg);
    
    virtual ~rootfile_database();
    
    /** \brief See documentation of Database::create_table
     */
    virtual std::auto_ptr<theta::Table> create_table(const std::string & table_name);
        
private:
    TFile * file;
    
    //which product columns to write as tree branch
    std::set<std::string> products_data;
    bool save_all_products;
    
    //histogram infos:
    struct hist_info{
        theta::Histogram1D h;
        std::string name;
        std::string column_name;
        hist_info(const theta::Histogram1D & h_, const std::string & name_, const std::string & colname_): h(h_), name(name_), column_name(colname_){}
    };
    std::vector<hist_info> hist_infos;
    
    
    friend class rootfile_table;
    //declare privately the table class:
    class rootfile_table: public theta::Table {
        friend class rootfile_database;
        class rootfile_column_int;
        
        boost::shared_ptr<rootfile_database> db;
        
        std::set<std::string> save_columns;
        bool save_all_columns;
        
        bool products_table;
        
        TTree * tree;
        int next_id;
        
        rootfile_table(const std::string & tablename, const boost::shared_ptr<rootfile_database> & db);
        virtual ~rootfile_table();
        
        theta::Column add_column(const std::string & name, const theta::data_type & type);
        virtual void add_row(const theta::Row & row);
        
        struct column_data{
            std::string name;
            theta::data_type type;
            union u_data{
                int i;
                double d;
                TString * s;
                TH1D * h;
            } data;
        };
        std::map<theta::Column, column_data> column_datas;
    };
};

#endif


