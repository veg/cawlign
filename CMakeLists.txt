cmake_minimum_required(VERSION 3.10)

project(cawlign)

#-------------------------------------------------------------------------------
# default installation prefix
#-------------------------------------------------------------------------------
set(INSTALL_PREFIX /usr/local CACHE PATH "Installation prefix")
set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE INTERNAL "Installation prefix" FORCE)

set(CMAKE_CONFIGURATION_TYPES Release)


# SET VERSION FROM FILE
file (STRINGS "VERSION.txt" VERSION_NUMBER)
add_definitions (-DVERSION_NUMBER=\"${VERSION_NUMBER}\")

include_directories(
    src/
)

add_executable(
    cawlign
    src/alignment.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
    src/argparse.cpp
    src/cawlign.cpp
    src/configparser.cpp
    src/scoring.cpp
)

target_compile_options (cawlign PRIVATE -fsigned-char -O3 -std=c++14  -funroll-loops )
target_link_options (cawlign PRIVATE -O3 -march=native -std=c++14 -fsigned-char -funroll-loops)


add_executable(
    cawlign_debug EXCLUDE_FROM_ALL
    src/alignment.cpp
    src/stringBuffer.cc
    src/tn93_shared.cc
    src/argparse.cpp
    src/cawlign.cpp
    src/configparser.cpp
    src/scoring.cpp
    
)

target_compile_options (cawlign_debug PRIVATE -O0 -std=c++14  -g -fsanitize=address)
target_link_options (cawlign_debug PRIVATE -O0 -std=c++14  -g -fsanitize=address)


set_property(
    TARGET cawlign
    APPEND PROPERTY COMPILE_DEFINITIONS LIBRARY_PATH="${INSTALL_PREFIX}/share/cawlign"
)





#add_definitions(${DEFAULT_COMPILE_OPTIONS})

find_package(OpenMP)

if(${OPENMP_FOUND})
   add_definitions(${OpenMP_CXX_FLAGS})
   target_link_libraries(cawlign PRIVATE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
   target_link_libraries(cawlign_debug PRIVATE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
endif(${OPENMP_FOUND})


install(
    TARGETS cawlign
	RUNTIME DESTINATION bin
	OPTIONAL
)

install(DIRECTORY res/ DESTINATION share/cawlign)

