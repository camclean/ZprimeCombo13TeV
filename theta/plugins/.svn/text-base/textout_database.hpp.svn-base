#ifndef PLUGIN_TEXTOUT_DATABASE_HPP
#define PLUGIN_TEXTOUT_DATABASE_HPP

#include "interface/decls.hpp"
#include "interface/database.hpp"
#include <boost/enable_shared_from_this.hpp>

#include <sqlite3.h>
#include <memory>
#include <string>
#include <set>


/** \brief Database which writes all information in a human-readable text form to stdout
 * 
 * This is mainly useful for debugging or if you want to see the result immediately.
 *
 * Configured via a setting group like
 * \code
 * output_database = {
 *   type = "textout_database";
 *   products_data = ("deltanll__nll_sb"); // optional, default is '*'
 * }
 * \endcode
 *
 * \c type must always be "sqlite_database" in order to select this plugin
 *
 * \c products_data is a list of column names to write. The default is to save all columns which can
 *    be configured by setting it to "*".
 */
class textout_database: public theta::Database{
public:
    
    /** \brief Constructor for the plugin system
     *
     * See class documentation for a description of the parsed Configuration settings.
     */
    textout_database(const theta::Configuration & cfg);
    
    
    virtual ~textout_database();
    
    /** \brief See documentation of Database::create_table
     */
    virtual std::auto_ptr<theta::Table> create_table(const std::string & table_name);
        
private:
    
    bool save_all_products;
    std::set<std::string> products_data;
    
    class textout_table: public theta::Table {
    friend class textout_database;
        
        virtual ~textout_table(){}
        
        virtual theta::Column add_column(const std::string & name, const theta::data_type & type);
        virtual void add_row(const theta::Row & row);

    private:
        
        textout_table(const std::string & name_, const boost::shared_ptr<textout_database> & db_);
        
        int next_colid;
        size_t irow;
        
        std::string name;
        
        boost::shared_ptr<textout_database> db;

        bool save_all_columns;
        std::set<std::string> save_columns;
        
        struct column_info{
            std::string name;
            theta::data_type type;
            column_info(){}
            column_info(const std::string & name_, const theta::data_type & type_): name(name_), type(type_){}
        };
        std::map<theta::Column, column_info> column_infos;
    };
};

#endif
