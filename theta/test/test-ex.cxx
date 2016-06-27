//this is a minimal FunctionPlugin used by test-exception.hpp
// to test inter-compilatio unit exception catching with shared objects.

#include "interface/phys.hpp"
#include "interface/plugin.hpp"
#include "interface/phys.hpp"

#include <string>

using namespace std;
using namespace theta;

class test_exception: public Function{
public:
    test_exception(const Configuration & cfg){}
    virtual double operator()(const ParValues & v) const{
       throw Exception("test-exception message 42");
    }
};

REGISTER_PLUGIN(test_exception)
