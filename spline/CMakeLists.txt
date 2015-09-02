cmake_minimum_required(VERSION 2.8)

find_package( libsimd REQUIRED )
find_package( cpplsq REQUIRED )

include ( CMakePackageConfigHelpers )

set ( CMAKE_CONFIG_DEST "share/cmake/libspline" )
set ( libspline_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. )

configure_package_config_file (                    
        "${CMAKE_CURRENT_SOURCE_DIR}/libspline-config.cmake.in"
        "${CMAKE_BINARY_DIR}/libspline-config.cmake"  
        INSTALL_DESTINATION "${CMAKE_CONFIG_DEST}"
        PATH_VARS 
          libsimd_INCLUDE_DIRS
          cpplsq_INCLUDE_DIRS
          libspline_INCLUDE_DIR )

export ( PACKAGE libspline )

set ( libspline_INCLUDE_DIR "include" )

configure_package_config_file ( 
       "${CMAKE_CURRENT_SOURCE_DIR}/libspline-config.cmake.in" 
       "${CMAKE_BINARY_DIR}/export/libspline-config.cmake"
       INSTALL_DESTINATION "${CMAKE_CONFIG_DEST}"
       PATH_VARS 
          libsimd_INCLUDE_DIRS
          cpplsq_INCLUDE_DIRS
          libspline_INCLUDE_DIR )

FILE(GLOB header_files "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
INSTALL(FILES ${header_files} DESTINATION include/spline)
INSTALL(FILES ${CMAKE_BINARY_DIR}/export/libspline-config.cmake DESTINATION ${CMAKE_CONFIG_DEST})