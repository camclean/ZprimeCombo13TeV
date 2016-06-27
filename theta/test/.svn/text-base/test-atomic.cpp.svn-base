#include <boost/test/unit_test.hpp>

#include "interface/atomic.hpp"

BOOST_AUTO_TEST_SUITE(atomic)


BOOST_AUTO_TEST_CASE(basic){
    theta::atomic_int a;
    atomic_set(&a, 0);
    BOOST_CHECK(atomic_get(&a) == 0);
}


BOOST_AUTO_TEST_SUITE_END()
