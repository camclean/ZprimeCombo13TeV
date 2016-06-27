#ifndef DATA_TYPE_HPP
#define DATA_TYPE_HPP

#include <ostream>

//note: data_type and Column are used by both, database and ProductSource. In order
// not to generate cyclic dependencies, provide the definition in this separate file.

namespace theta{
    enum data_type { typeDouble, typeInt, typeString, typeHisto };


/** \brief Class representing a single Column within a Table
 * 
 * Column is a thin wrapper around an \c int called id.
 * Instances should always be retrieved via Table::add_column.
 */
class Column{
public:
    /// Comparison operator
    bool operator<(const Column & rhs) const{
        return id < rhs.id;
    }

    bool operator==(const Column & rhs) const{
        return id == rhs.id;
    }
    
    /// Construct a Column instance, given the id
    explicit Column(int id_ = -1): id(id_){}
    
    /// Get the id as set at construction time
    int get_id() const{
        return id;
    }
private:
    int id;
};

inline std::ostream & operator<<(std::ostream & out, const Column & c){
    return out << "Column(" << c.get_id() << ")";
}



}




#endif

