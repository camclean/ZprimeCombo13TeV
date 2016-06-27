#ifndef PLUGIN_SQLITE_DATABASE_HPP
#define PLUGIN_SQLITE_DATABASE_HPP

#include "interface/decls.hpp"
#include "interface/database.hpp"
#include <boost/enable_shared_from_this.hpp>

#include <sqlite3.h>
#include <memory>
#include <string>
#include <sstream>
#include <set>


/** \brief Database which stores all information in a single sqlite3 database file
 *
 * Configured via a setting group like
 * \code
 * output_database = {
 *   type = "sqlite_database";
 *   filename = "abc.db";
 *   products_data = ("deltanll__nll_sb"); // optional, default is '*'
 * }
 * \endcode
 *
 * \c type must always be "sqlite_database" in order to select this plugin
 *
 * \c filename is the filename of the sqlite3 output file. It is a path relative to the path where theta is invoked
 *
 * \c products_data is a list of column names to save. The default is to save all columns which can
 *    be configured by setting it to "*". Note that the 'runid' and 'eventid' columns will always be saved.
 *
 * If the file already exists, it is overwritten silently.
 *
 * The types theta::typeDouble, theta::typeInt and theta::typeString are translated directly
 * to their SQL counterparts \c DOUBLE, \c INT(4) and \c TEXT, respectively. For theta::typeHisto,
 * an SQL BLOB is saved which contains the lower and upper border of the histogram and the raw histogram data,
 * including one underflow bin and one overflow bin, which are always set to zero.
 * 
 * Note that due to limits in sqlite3, more than 999 columns per table are not supported. To allow
 * more columns nevertheless, this plugin creates secondary (and tertiary, ...) tables in the sqlite3 database
 * which store additional columns with the name &lt;table name&gt;__&lt;n&gt; where n is an integer starting to count at 0.
 * (Note that these tables are only created as necessary, i.e. if all tables have less than 999 columns, each theta::Table
 * corresponds to one sqlite3 table of the same name). The connection between corresponding rows in those tables is done with the
 * special column "rowid" which is always implicitly created by sqlite3 for all tables. So in principle, one
 * could quey the complete column with a statement like "select * from products left join products__0 on rowid=products__0.rowid;", although
 * this query might hit other limits in sqlite3, and the safest choice is to query each table in a separate statement, using the "rowid"
 * from the first query result as SQL "where" condition.
 */
class sqlite_database: public theta::Database{
public:
    
    /** \brief Constructor for the plugin system
     *
     * See class documentation for a description of the Configuration settings.
     */
    sqlite_database(const theta::Configuration & cfg);
    
    
    virtual ~sqlite_database();
    
    /** \brief See documentation of Database::create_table
     */
    virtual std::auto_ptr<theta::Table> create_table(const std::string & table_name);
        
private:
    int maxcol_per_table; // maximum number of columns as supported by sqlite
    
    /** Execute the sql string \c query.
     */
    void exec(const std::string & query);
    
    /** Prepare the sql statement \c query.
     */
    sqlite3_stmt* prepare(const std::string & query);
    
    /** Start a Transaction. In case of an error, a DatabaseException is thrown.
     */
    void beginTransaction();
    
    /** End a transaction previously started with \c beginTransaction. Calling
     *  endTransaction without a call to beginTransaction is valid (it is a no-op).
     * */
    void endTransaction();

    /** Throw a DatabaseException using the last database error message.
     * 
     * This method is intended for Table objects which want to propagate
     * an error via Exceptions in a consistent manner.
     * 
     * \param functionName: name of the function (in a \c Table object) in which the database error ocurred.
     */
    void error(const std::string & functionName);
    
    void close();
    
    sqlite3* db;
    bool transaction_active;
    bool save_all_products;
    std::set<std::string> products_data;
    
    //declare privately(!) the sqlite_table class:
    class sqlite_table: public theta::Table {
    friend class sqlite_database;

        // destructor; creates the table if empty
        virtual ~sqlite_table();
        
        virtual theta::Column add_column(const std::string & name, const theta::data_type & type);

        virtual void add_row(const theta::Row & row);

    private:
        
        sqlite_table(const std::string & name_, const boost::shared_ptr<sqlite_database> & db_);
        
        std::string name;
        bool table_created;
        
        int next_insert_index; // next free insert_index to use by add_column to construct an sqlite_column, starting at 1.
        
        std::vector<sqlite3_stmt*> insert_statements; // these are owned by sqlite_database.
        boost::shared_ptr<sqlite_database> db;

        bool save_all_columns;
        std::set<std::string> save_columns;
        
        void create_table();
        
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
