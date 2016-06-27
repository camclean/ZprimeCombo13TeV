#include "interface/atomic.hpp"
#include <boost/static_assert.hpp>

// For the assembler version to work correctly, the atomic_int has to be
// aligned at an 8-byte-boundary and the address for atomic_int
// must coincide with the address for atomic_int::value. This macros
// check these two conditions:
BOOST_STATIC_ASSERT(sizeof(theta::atomic_int)==8);
BOOST_STATIC_ASSERT(__alignof__(theta::atomic_int)==8);

