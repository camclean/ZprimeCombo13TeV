#include "interface/atomic.hpp"

#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>

// This is the "generic" implementation of the functions in atomic.hpp, a replacement for atomic.s
// for non-x86_64 platforms.
// This will not have a great performance, especially if using more than a handful of atomic_ints,
// as we use one global mutex for *all* atomic_ints.
// Note that we work with intermediate volatiles here because otherwise, link-time-optimization could optimize
// away non-volatile accesses to atomic_int::value.

namespace{
   boost::mutex m;
}

void atomic_add(theta::atomic_int * a, uint64_t addend){
    boost::unique_lock<boost::mutex>(m);
    a->value += addend;
}

void atomic_inc(theta::atomic_int * a){
    boost::unique_lock<boost::mutex>(m);
    ++a->value;
}

uint64_t atomic_get(const theta::atomic_int * a){
    boost::unique_lock<boost::mutex>(m);
    return a->value;
}

void atomic_set(theta::atomic_int * a, uint64_t newval){
    boost::unique_lock<boost::mutex>(m);
    a->value = newval;
}
