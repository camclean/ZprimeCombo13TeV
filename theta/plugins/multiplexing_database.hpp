#ifndef PLUGIN_MULTIPLEXING_DATABASE_HPP
#define PLUGIN_MULTIPLEXING_DATABASE_HPP

#include "interface/decls.hpp"
#include "interface/database.hpp"

#include <memory>
#include <string>


/** \brief Database which saves information to several underlying databases
 *
 * This is useful to save information to more than one database at once, for example
 * in a sqlite and root file at the same time, or displaying the result on stdout and saving
 * it to a file, too.
 *
 * Configured via
 * \code
 * output_database = {
 *   type = "multiplexing_database";
 *   databases = ("@db1", "@db2");
 * };
 * db1 = { ... }; // some database configuration
 * db2 = { ... }; // some other database configuration
 * \endcode
 *
 * \c type must always be "multiplexing_database" in order to select this plugin
 *
 * \c databases is a list of databases which are used to write out results to
 */
class multiplexing_database: public theta::Database{
public:
    
    /** \brief Constructor for the plugin system
     */
    multiplexing_database(const theta::Configuration & cfg);
    
    virtual ~multiplexing_database();
    virtual std::auto_ptr<theta::Table> create_table(const std::string & table_name);
        
private:
    std::vector<boost::shared_ptr<theta::Database> > underlying_dbs;

    class multiplexing_table: public theta::Table {
        friend class multiplexing_database;
        std::map<theta::Column, std::vector<theta::Column> > external_to_underlying_cols;
        std::map<theta::Column, theta::data_type> column_types;
        boost::ptr_vector<theta::Table> underlying_tables;
        int next_icol;
        
        multiplexing_table(const boost::shared_ptr<Database> & db, boost::ptr_vector<theta::Table> & underlying_tables);

     public:
        virtual ~multiplexing_table();
        theta::Column add_column(const std::string & name, const theta::data_type & type);
        virtual void add_row(const theta::Row & row);
    };
};

#endif

