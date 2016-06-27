file(GLOB plugin_sources *.cpp)
add_library(plugin MODULE ${plugin_sources})
target_link_libraries(plugin libtheta sqlite3 ${Boost_LIBRARIES})
add_custom_command(TARGET plugin POST_BUILD COMMAND mv ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libplugin.so ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/core-plugins.so)
