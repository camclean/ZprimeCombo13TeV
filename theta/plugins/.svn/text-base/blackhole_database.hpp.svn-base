#ifndef PLUGIN_BLACKHOLE_DATABASE_HPP
#define PLUGIN_BLACKHOLE_DATABASE_HPP

#include "interface/decls.hpp"
#include "interface/database.hpp"

#include <memory>
#include <string>


/** \brief Database which discards all information
 *
 * This is only useful for testing and benchmarking purposes.
 *
 * Configured via a setting group like
 * \code
 * output_database = {
 *   type = "blackhole_database";
 * };
 * \endcode
 *
 * \c type must always be "blackhole_database" in order to select this plugin
 */
class blackhole_database: public theta::Database{
public:
    
    /** \brief Constructor for the plugin system
     */
    blackhole_database(const theta::Configuration & cfg);
    
    virtual ~blackhole_database();
    virtual std::auto_ptr<theta::Table> create_table(const std::string & table_name);
        
private:
    //declare privately(!) the sqlite_table class:
    class blackhole_table: public theta::Table {
        friend class blackhole_database;
        blackhole_table(const boost::shared_ptr<Database> & db): Table(db){}
        virtual ~blackhole_table();
        theta::Column add_column(const std::string & name, const theta::data_type & type);
        virtual void add_row(const theta::Row & row);

    };
};

#endif
