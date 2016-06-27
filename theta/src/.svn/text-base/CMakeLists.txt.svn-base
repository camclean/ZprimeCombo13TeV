file(GLOB libtheta_sources *.cpp)

if(generic)
   set(libtheta_sources ${libtheta_sources} atomic.cxx)
else()
   set(libtheta_sources ${libtheta_sources} atomic.s atomic_asserts.cxx)
endif()

add_library(libtheta SHARED ${libtheta_sources})
target_link_libraries(libtheta m dl config ${Boost_LIBRARIES})

if(amdlibm)
   target_link_libraries(libtheta ${amdlibm_lib})
endif()

set_property(TARGET libtheta PROPERTY OUTPUT_NAME theta)

