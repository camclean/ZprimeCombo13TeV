#include "interface/utils.hpp"
#include "test/utils.hpp"

#include <boost/test/unit_test.hpp>

using namespace theta::utils;
using namespace std;

BOOST_AUTO_TEST_SUITE(utils_tests)


BOOST_AUTO_TEST_CASE(roots_quad_tests){
    double b = 2, c = -15;
    double x1, x2;
    
    // 2 solutions:
    int res = roots_quad(x1, x2, b, c);
    BOOST_CHECK(res == 2);
    BOOST_CHECK(x1 == -5);
    BOOST_CHECK(x2 == 3);
    
    // 1 solution:
    b = -2;
    c = 1;
    res = roots_quad(x1, x2, b, c);
    BOOST_CHECK(res == 1);
    BOOST_CHECK(x1 == 1);
    BOOST_CHECK(x2 == 1);
    
    // no solutions:
    b = 1;
    c = 2;
    res = roots_quad(x1, x2, b, c);
    BOOST_CHECK(res == 0);
    BOOST_CHECK(std::isnan(x1));
    BOOST_CHECK(std::isnan(x2));
    
    // b == 0.0:
    b = 0;
    c = -4;
    res = roots_quad(x1, x2, b, c);
    BOOST_CHECK(res == 2);
    BOOST_CHECK(close_to_relative(x1,-2.0));
    BOOST_CHECK(close_to_relative(x2,+2.0));
}


BOOST_AUTO_TEST_SUITE_END()
