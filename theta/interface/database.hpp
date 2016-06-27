#ifndef DATABASE_HPP
#define DATABASE_HPP

#include "interface/decls.hpp"
#include "interface/data_type.hpp"
#include "interface/histogram.hpp"
#include "interface/producer.hpp"


#include <string>
#include <map>
#include <memory>

#include <boost/utility.hpp>
#include <boost/enable_shared_from_this.hpp>

namespace theta {

/** \brief Abstract class for data input
 *
 * A database in theta is a collection of tables. Tables are identified by a unique name
 * within the database and can contain an arbitrary number of columns (identified by name) of type double, int, string and
 * Histogram.
 *
 * This is an abstract class. Concrete instances can be retrieved via the plugin system using
 * this class a template argument to PluginManager.
 */
class DatabaseInput: private boost::noncopyable{
public:

    /// Required for the plugin system
    typedef DatabaseInput base_type;
    
    /** \brief Iterator-like class to step through the result of a Database query
     *
     * A typical use of the iterator is:
     * \code
     * DatabaseInput & db = ...;
     * vector<string> cols;
     * ... // fill cols
     * std::auto_ptr<ResultIterator> it = db.query("some_table", cols);
     * while(it->has_data()){
     *    // use row data via it->get_* methods
     *    ++(*it);
     * }
     * \endcode
     * This works also correctly in case of zero rows.
     */
    class ResultIterator: private boost::noncopyable{
    public:
        /** \brief Retrieve next row in the table
         *
         * If there is no next row, a subsequent call to has_data() will return false.
         *
         * In case of an error, a DatabaseException is thrown.
         */
        virtual void operator++() = 0;
        
        /** \brief Returns true iff this iterator points to a valid column containing data
         *
         * Returns false if there are no more result columns
         */
        virtual bool has_data() = 0;
        
        ///@{
        /** \brief Retrieve the column values of the current row
         *
         * Before calling any of these, make sure has_data() returns true.
         *
         * The column index argument is zero-based and refers to the column_names vector
         * passed to DatabaseInput::query.
         *
         * In case of type mismatch, a DatabaseException is thrown.
         */
        virtual double get_double(size_t icol) = 0;
        virtual int get_int(size_t icol) = 0;
        virtual theta::Histogram1D get_histogram(size_t icol) = 0;
        virtual std::string get_string(size_t icol) = 0;
        ///@}
        
        virtual ~ResultIterator();
    };

    /** \brief create a row-iterator, retrieving some columns for all rows in a table
     *
     * The returned ResultIterator points to the first result row.
     * In case of an error, such as if the column is not found, a DatabaseException is thrown.
     */
    virtual std::auto_ptr<ResultIterator> query(const std::string & table_name, const std::vector<std::string> & column_names) = 0;

    /// Get all table names in the database
    virtual std::vector<std::string> get_all_tables() = 0;

    /// Get all column names and types for the given table
    virtual std::vector<std::pair<std::string, data_type> > get_all_columns(const std::string & table_name) = 0;

    virtual ~DatabaseInput();
};


/** \brief Abstract database class for data output
 *
 * A database in theta is a collection of tables. Tables are identified by a unique name
 * within the database and can contain an arbitrary number of columns (identified by a name string)
 * of type double, int, string and Histogram.
 *
 * This is an abstract class. Concrete instances can be retrieved via the plugin system using
 * this class a template argument to PluginManager.
 */
class Database: private boost::noncopyable, public boost::enable_shared_from_this<Database> {
public:
    
    /// Required for the plugin system
    typedef Database base_type;

    /** \brief Virtual Destructor
     *
     * Virtual, as polymorphic access to derived classes will happen.
     *
     * Should do any cleanup work (like closing the database file for file-based databases,
     * closing the network connection for network-based, etc.).
     */
    virtual ~Database();
    
    /** \brief Create a new table instance within this database
     *
     * Table names must start with a letter and otherwise consist only of letters,
     * digits and underscores. Throws an InvalidArgumentException if the name does not meet the requirements.
     *
     * The returned Table can be used only as long as this Database instance is not destroyed.
     * Using Tables created from a destroyed Database yields undefined behavior.
     */
    virtual std::auto_ptr<Table> create_table(const std::string & table_name) = 0;

protected:
    
    /** \brief Checks the table name requirements
     *
     * Derived classes can call this at the beginning of their implementation of create_table
     * in order to meet the specification: this method will throw an InvalidArgumentException
     * if \c name violates the specification.
     */
    void check_name(const std::string & table_name);
};


/** \brief Data for a single row in a Table
 * 
 * This class represents (the data of) a single row within a theta::Table. It exposes getter
 * and setter methods to get/set the values of a column identified by an instance
 * of theta::Column.
 */
class Row{
public:
    ///@{
    /// Setters for the row values, given a Column instance
    void set_column(const Column & col, double d);
    void set_column(const Column & col, int i);
    void set_column(const Column & col, const std::string & s);
    void set_column(const Column & col, const Histogram1D & h);
    ///@}
    
    ///@{ Getters for the row values, given a Column instance
    /** \brief Setters for the row values, given a Column instance
     * 
     * Will throw a DatabseException if the Column value asked for has not been
     * set previously by a (type-matching) setter method.
     */
    double get_column_double(const Column & col) const;
    int get_column_int(const Column & col) const;
    const std::string & get_column_string(const Column & col) const;
    const Histogram1D & get_column_histogram(const Column & col) const;
    ///@}

    void clear();
private:
    std::map<Column, double> doubles;
    std::map<Column, int> ints;
    std::map<Column, std::string> strings;
    std::map<Column, Histogram1D> histos;
};


/** \brief Abstract class for a table in a Database
 *
 * Tables are always constructed via a Database instance. Once created, it can be used by:
 * <ol>
 *   <li>Call add_column one or more times and save the returned Column handles</li>
 *   <li>To write out a row, create a Row instance and call Row::set_column,
 *      using the Column handles from the previous step. After the Row instance
 *      is fully populated (=has values for all Columns in this Table), call Table::add_row.</li>
 * </ol>
 *
 * You must not call add_column after the first call to add_row. Roes added via add_row must
 * provide data for all previously defined columns.
 */
/* Implementation notes: there used to be an old interface for Table consisting of
 * set_column(Column &, const T & value)    [T=double, int, string, Histogram]
 * and
 * add_row()
 * methods; a user of this class would first call set_column for all columns; a subsequent call to add_row would
 * add a row with the latest set values.
 * Instead, the new interface encapsulates the data state for a new row into a Row instance. This is done
 * for easier implementation of concurrency: different threads can populate their own Row instance and add it via a single call
 * to add_row, which is much easier to synchronize than a block of set_column calls, followed by a add_row.
 */
class Table: private boost::noncopyable {
public:

    /// destructor; creates the table if empty
    virtual ~Table();
    
    /** \brief Add a column to the table
     *
     * Create a new column of the given name and type and return the corresponding
     * Column handle which can be used for the methods in the Row class.
     *
     * This method should only be called at the beginning, after the construction of the object and
     * before using it to actually write any data.
     */
    virtual Column add_column(const std::string & name, const data_type & type) = 0;
    
    /** \brief Add a row to the table
     * 
     * row must provide a value for all Columns previously created with add_column.
     */
    virtual void add_row(const Row & row) = 0;
protected:
    //the tables always hold a shared ptr to the database to prevent
    // the database being destroyed earlier than all its tables (!)
    boost::shared_ptr<Database> db;
    Table(const boost::shared_ptr<Database> & db_);
private:
    Table(); //not implemented
};


/** \brief A Table to store products, i.e., per-event output
 *
 * Per \link theta::Run Run \endlink, there is exactly one ProductsTable. Products are produced by
 * <ul>
 *   <li>Producer instances to save the (per-event) result of their computation</li>
 *   <li>DataSource to save some per-event information about data production</li>
 * </ul>
 *
 * An ProductsTable will have an integer column named "runid" and an integer column named "eventid". Additionally,
 * any columns defined by the Producer or DataSource, with column names as described in ProductsTable::add_column.
 *
 * Clients use ProductsTable very similarly to a Table. Differences are (i) the
 * signature of the add_column method which takes an additional \c name argument and (ii)
 * the add_row column, which takes a Run instance as argument here.
 *
 * The actual write is done by the \link Run \endlink instance; the Producer / DataSource must not
 * call add_row directly.
 */
class ProductsTable: public ProductsSink {
public:
    virtual void set_product(const Column & c, double d);
    virtual void set_product(const Column & c, int i);
    virtual void set_product(const Column & c, const std::string & s);
    virtual void set_product(const Column & c, const theta::Histogram1D & histo);
    
    /** \brief Construct a new ProductsTable based on the given table
     *
     * Memory ownership of table will be transferred to this.
     */
    ProductsTable(std::auto_ptr<Table> & table);
    
    /** \brief Add a row to the table, given the current run
     *
     * This is called by theta::Run after a producer has executed.
     */
    void add_row(int runid, int eventid);

    virtual ~ProductsTable(){}
        
private:
    virtual Column declare_column_impl(const std::string & full_product_name, const data_type & type);

    std::auto_ptr<Table> table;
    Row current_row;
    Column c_runid, c_eventid;
};



/** \brief Table to store all logging information
 *
 * The corresponding table has following columns:
 * <ol>
 * <li>runid (typeInt)</li>
 * <li>eventid (typeInt)</li>
 * <li>severity (typeInt)</li>
 * <li>message (typeString)</li>
 * <li>time (typeDouble)</li>
 * </ol>
 * 
 * \c runid and \c eventid are the run and event id, respectively, the log entry
 * is associated to. \c eventid is set to 0, if no particular event (but the run as a whole)
 * is referred to.
 * 
 * \c severity is one level from LogTable::e_severity. See there for a description of the meaning
 * of these levels.
 * 
 * \c message is the human-readable log message.
 * 
 * \c time is the number of seconds since the unix epoch (1970-01-01 UTC),
 *    with sub-second accuracy.
 **/
class LogTable: private boost::noncopyable {
public:
    /** \brief Severity levels for log messages
     * 
     * <ol>
     * <li>\c error should be used if a serious condition is reported which will likely
     *     affect the whole result.</li>
     * <li>\c warning should be used if a problem is not as serious as \c error, but might affect
     *     the result in an undesired way. </li>
     * <li>\c info is used purely informational, without the implicit action request of
     *    \c error and \c warning. It covers such events like the start and end of a pseudo experiment. </li>
     * <li>\c debug is used like \c info, but for very detailed reporting which is usually not required.</li>
     * </ol>
     **/
    enum e_severity {
        error = 0, warning, info, debug
    };
    
    /** \brief Construct a new logTable based on the given Table
     *
     * The default loglevel is warning.
     *
     * Ownership of table will be transferred.
     */
    LogTable(std::auto_ptr<Table> & table);
    
    /** \brief Set the current log level.
     *
     * Future calls to append() will only write the message to the table if the messages
     * has a severity which is equal to or exceeds the level given here.
     *
     * Note that is it not possible to disable logging of error messages.
     */
    void set_loglevel(e_severity s);
    
    /** \brief Get the currently active log level
     */
    e_severity get_loglevel() const;
    
    /** \brief Append message to log table, if severity is larger than currently configured level
     * 
     *
     * \c s is the severity level of the log message
     *
     * \c message is the message, in human-readable english
     */
    void append(int runid, int eventid, e_severity s, const std::string & message){
        //define inline as hint to the compiler; keep this function as short as possible to
        // encourage the compiler actually inlining it and to have high performance benefits in
        // case the user disables logging.
        if(s <= level) really_append(runid, eventid, s, message);
    }
    
    /** \brief Returns the number of messages
     *
     * Use severity::e_severity converted to int to access the number of messages written to the log.
     * Messages suppressed by the set_loglevel() method <em>not</em> counted.
     */
    const int* get_n_messages() const;

private:

    /// really append the log message. Called from append() in case severity is large enough
    void really_append(int runid, int eventid, e_severity s, const std::string & message);
    e_severity level;
    int n_messages[4];
    Column c_runid, c_eventid, c_severity, c_message, c_time;
    std::auto_ptr<Table> table;
};


/** \brief Table to store information about the random number seeds.
 *
 * There is a one-to-one relationship between this class and \link theta::Run Run \endlink, i.e.,
 * to each RndInfoTable instance, there is one Run instance and vice versa.
 *
 * The corresponding table has following columns:
 * <ol>
 * <li>runid (typeInt): the run id the entry refers to.</li>
 * <li>seed (typeInt): the seed of the random number generator used in the run.</li>
 * </ol>
 */
class RndInfoTable: private boost::noncopyable {
public:
    /** \brief Construct with an underlying Table
     *
     * Ownership of table will be transferred.
     */
    RndInfoTable(std::auto_ptr<Table> & table);

    /** \brief append an entry to the RndInfoTable
     *
     * \c runid is the current runid
     * \c name is the name of the module of the seed, according to \link theta::ProductsTableWriter ProductsTableWriter \endlink . <br />
     * \c seed is the seed of the random number generator used for this module
     */
    void append(int runid, const std::string & name, int seed);
private:
    Column c_runid, c_name, c_seed;
    std::auto_ptr<Table> table;
};

}

#endif
