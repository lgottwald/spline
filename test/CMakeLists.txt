cmake_minimum_required(VERSION 2.8)

include(ExternalProject)

ExternalProject_Add(libcatch
  DOWNLOAD_DIR "${CMAKE_CURRENT_BINARY_DIR}/catch"
  DOWNLOAD_COMMAND wget https://raw.githubusercontent.com/catchorg/Catch2/Catch1.x/single_include/catch.hpp  
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND ""
)

find_package( libspline REQUIRED)
FIND_PACKAGE( Boost COMPONENTS serialization system REQUIRED )

include_directories(
  ${Boost_INCLUDE_DIRS} 
  ${libspline_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR}
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic-errors -Wall -Wextra -Wno-unused-parameter")

enable_testing()
add_test( SplineTest spline_test )
add_executable( spline_test HornerSchemeTest.cpp KnotTransformTest.cpp PiecewisePolynomialTest.cpp PolynomialRootsTest.cpp Main.cpp )
target_link_libraries( spline_test ${Boost_LIBRARIES} )
add_dependencies(spline_test libcatch )
