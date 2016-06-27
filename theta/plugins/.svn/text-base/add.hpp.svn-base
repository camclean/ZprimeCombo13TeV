#include "interface/phys.hpp"

/** \brief A Function which adds different addends which can be parameters, literal values and other Functions
 *
 * It is the counterpart of the \link multilply multiply \endlink plugin.
 *
 * Example configuration:
 * \code
 * function1 = {
 *   type = "add";
 *   addends = ("p1", "p2", 1.7, 2.3, "@function2");
 * }
 * 
 * function2 = {  ... some function specification ...   };
 * \endcode
 *
 * \c addends is a list of addends. Each addend can be a literal floating point value will be included directly,
 *   a string which is interpreted as parameter name which will be included, or a setting group
 *   used to construct a Function object via the plugin system.
 */
class add: public theta::Function{
public:
    add(const theta::Configuration & cfg);
    virtual double operator()(const theta::ParValues & v) const;
private:
    std::vector<theta::ParId> v_pids;
    double literal_addend;
    boost::ptr_vector<theta::Function> functions;
};

