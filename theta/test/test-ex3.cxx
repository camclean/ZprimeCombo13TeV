//test exceptions during Factory::build

#include "interface/phys.hpp"
#include "interface/plugin.hpp"

#include <string>

using namespace std;
using namespace theta;

class test_ex_during_build: public Function{
public:
  test_ex_during_build(const Configuration & cfg){
       throw Exception("exception message 23");
  }
  virtual double operator()(const ParValues & v) const{       return 0.0;  }
};

REGISTER_PLUGIN(test_ex_during_build)
