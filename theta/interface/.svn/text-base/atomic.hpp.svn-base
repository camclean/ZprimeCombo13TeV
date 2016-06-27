#ifndef ATOMIC_HPP
#define ATOMIC_HPP

#include <stdint.h>
#include "interface/decls.hpp"


extern "C" {

///@{ \brief Functions to access and manipulate theta::atomic_int
void atomic_add(theta::atomic_int * a, uint64_t addend);
void atomic_inc(theta::atomic_int * a);
uint64_t atomic_get(const theta::atomic_int * a);
void atomic_set(theta::atomic_int * a, uint64_t value);
///@}

}

namespace theta{


/** \brief 64-bit unsigned integer value for atomic operations
 *
 * This class has only private members; use the functions atomic_* to
 * read and manipulate the internal state.
 */
// note: keep this class to the minimum, so it's a POD and has guaranteed alignment and size
// on as many platformas as possible ...
class atomic_int{
void friend ::atomic_add(theta::atomic_int * a, uint64_t addend);
void friend ::atomic_inc(theta::atomic_int * a);
uint64_t friend ::atomic_get(const theta::atomic_int * a);
void friend ::atomic_set(theta::atomic_int * a, uint64_t value);
private:
    uint64_t value;
};

}

#endif

