# Defines targets for the C++ restraints implemented here. These CMake targets are used by the
# unit tests and by the Python module target defined in ../pythonmodule/CMakeLists.txt

# Create a shared object library for our restrained ensemble plugin.
add_library(mdstring_potential STATIC
            mdstring_potential.h
            mdstring_potential.cpp
            sessionresources.cpp)
set_target_properties(mdstring_potential PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(mdstring_potential PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:include>
                           )
# RPATH management
# don't skip the full RPATH for the build tree
set_target_properties(mdstring_potential PROPERTIES SKIP_BUILD_RPATH FALSE)

# (but later on when installing)
# If building with setuptools, CMake will not be performing the install
set_target_properties(mdstring_potential PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)

target_link_libraries(mdstring_potential PRIVATE Gromacs::gmxapi)
