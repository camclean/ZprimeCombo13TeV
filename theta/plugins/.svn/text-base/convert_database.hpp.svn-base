#include "interface/decls.hpp"
#include "interface/main.hpp"

#include <boost/shared_ptr.hpp>

/** \brief Read in a data from a theta::DatabaseInput and write all data out via a theta::Database
 *
 * This is useful for merging databases or converting databases from one to another format.
 *
 * Configuration via
 * \code
 * main = {
 *   type = "convert_database";
 *   input_database = {     // some DatabaseInput specification, for example:
 *     type = "sqlite_database_in";
 *     filenames = ("a.db", "b.db");
 *   };
 *   output_database = {    // some Database specification, for example:
 *       type = "rootfile_database";
 *       filename = "ab.root";
 *   };
 * };
 * \endcode
 */
class convert_database: public theta::Main{
public:
    virtual void run();
    convert_database(const theta::Configuration & cfg);
    virtual ~convert_database();
private:
    boost::shared_ptr<theta::DatabaseInput> db_in;
    boost::shared_ptr<theta::Database> db_out;
};


