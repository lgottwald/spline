cmake_minimum_required(VERSION 2.8)

PROJECT(libspline)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

include(ExternalProject)

set(LOCAL_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/install )
set(LOCAL_PREFIX_PATH ${CMAKE_INSTALL_PREFIX} ${LOCAL_INSTALL_PREFIX} ${CMAKE_PREFIX_PATH})
string(REPLACE ";" "^^" LOCAL_PREFIX_PATH "${LOCAL_PREFIX_PATH}")

ExternalProject_Add(
  libspline
  LIST_SEPARATOR ^^
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/spline
  UPDATE_COMMAND ""
  INSTALL_COMMAND make install > ${LOCAL_INSTALL_PREFIX}/install_output.log
  INSTALL_DIR ${LOCAL_INSTALL_PREFIX}
  CMAKE_ARGS -DCMAKE_PREFIX_PATH=${LOCAL_PREFIX_PATH} -DCMAKE_INSTALL_PREFIX=${LOCAL_INSTALL_PREFIX} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
)

ExternalProject_Add_Step(libspline forceconfigure
   COMMAND ${CMAKE_COMMAND} -E echo "Force configure of libspline"
   DEPENDEES update
   DEPENDERS configure
   ALWAYS 1)


find_package( libsimd QUIET )
if( NOT libsimd_FOUND )
 ExternalProject_Add(
  libsimd
  LIST_SEPARATOR ^^
  GIT_REPOSITORY https://github.com/rgottwald/simd.git
  GIT_TAG origin/master
  INSTALL_DIR ${LOCAL_INSTALL_PREFIX}
  INSTALL_COMMAND make install > ${LOCAL_INSTALL_PREFIX}/install_output.log
  CMAKE_ARGS -DCMAKE_PREFIX_PATH=${LOCAL_PREFIX_PATH} -DCMAKE_INSTALL_PREFIX=${LOCAL_INSTALL_PREFIX} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
 )
 add_dependencies(libspline libsimd)
endif()


ExternalProject_Add(
  libspline-test
  LIST_SEPARATOR ^^
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test
  TEST_BEFORE_INSTALL 1
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  CMAKE_ARGS -DCMAKE_PREFIX_PATH=${LOCAL_PREFIX_PATH}
  DEPENDS libspline
)

find_package( cpplsq QUIET )
if( NOT cpplsq_FOUND )
 ExternalProject_Add(
  libcpplsq
  LIST_SEPARATOR ^^
  GIT_REPOSITORY https://github.com/rgottwald/cpplsq.git
  GIT_TAG origin/master
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  CMAKE_ARGS -DCMAKE_PREFIX_PATH=${LOCAL_PREFIX_PATH} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
 )
 add_dependencies( libspline libcpplsq )
endif()

INSTALL(
  DIRECTORY
    ${LOCAL_INSTALL_PREFIX}/include
    ${LOCAL_INSTALL_PREFIX}/share
  DESTINATION . )