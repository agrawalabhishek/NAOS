# Copyright (c) <year> <author> (<email>)
# Distributed under the MIT License.
# See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT

# Include script to build external libraries with CMake.
include(ExternalProject)

# -------------------------------

# Catch: https://github.com/philsquared/Catch

if(BUILD_TESTS)
  if(NOT BUILD_DEPENDENCIES)
    find_package(CATCH)
  endif(NOT BUILD_DEPENDENCIES)

  if(NOT CATCH_FOUND)
    message(STATUS "Catch will be downloaded when ${CMAKE_PROJECT_NAME} is built")
    ExternalProject_Add(catch-lib
      PREFIX ${EXTERNAL_PATH}/Catch
      #--Download step--------------
      URL https://github.com/philsquared/Catch/archive/master.zip
      TIMEOUT 30
      #--Update/Patch step----------
      UPDATE_COMMAND ""
      PATCH_COMMAND ""
      #--Configure step-------------
      CONFIGURE_COMMAND ""
      #--Build step-----------------
      BUILD_COMMAND ""
      #--Install step---------------
      INSTALL_COMMAND ""
      #--Output logging-------------
      LOG_DOWNLOAD ON
    )
    ExternalProject_Get_Property(catch-lib source_dir)
    set(CATCH_INCLUDE_DIRS ${source_dir}/include CACHE INTERNAL "Path to include folder for Catch")
  endif(NOT CATCH_FOUND)

  if(NOT APPLE)
    include_directories(SYSTEM AFTER "${CATCH_INCLUDE_DIRS}")
  else(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${CATCH_INCLUDE_DIRS}\"")
  endif(NOT APPLE)
endif(BUILD_TESTS)

# -------------------------------

if(NOT BUILD_DEPENDENCIES)
  find_package(GSL)
endif(NOT BUILD_DEPENDENCIES)

if(NOT GSL_FOUND)
  message(STATUS "GSL will be downloaded when ${CMAKE_PROJECT_NAME} is built")
  ExternalProject_Add(gsl-lib
    PREFIX ${EXTERNAL_PATH}/GSL
    #--Download step--------------
    URL https://github.com/ampl/gsl/archive/master.zip
    TIMEOUT 30
    #--Update/Patch step----------
    #--Configure step-------------
    #--Build step-----------------
    BUILD_IN_SOURCE 1
    #--Install step---------------
    INSTALL_COMMAND ""
    #--Output logging-------------
    LOG_DOWNLOAD ON
  )
  ExternalProject_Get_Property(gsl-lib source_dir)
  set(GSL_INCLUDE_DIRS ${source_dir} CACHE INTERNAL "Path to include folder for SGP4")
  set(GSL_LIBRARY_DIRS ${source_dir} CACHE INTERNAL "Path to library folder for SGP4")
  set(GSL_LIBRARIES "gsl" "gslcblas")
endif(NOT GSL_FOUND)

if(NOT APPLE)
  include_directories(SYSTEM AFTER "${GSL_INCLUDE_DIRS}")
else(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${GSL_INCLUDE_DIRS}\"")
endif(NOT APPLE)
link_directories(${GSL_LIBRARY_DIRS})

# -------------------------------

find_package(Threads)

if(NOT BUILD_DEPENDENCIES)
  find_package(sqlite3)
endif(NOT BUILD_DEPENDENCIES)

if(NOT SQLITE3_FOUND)
  message(STATUS "SQLite3 will be downloaded when ${CMAKE_PROJECT_NAME} is built")
  ExternalProject_Add(sqlite3-lib
    PREFIX ${EXTERNAL_PATH}/SQLite3
    #--Download step--------------
    URL https://github.com/kartikkumar/sqlite3-cmake/archive/master.zip
    TIMEOUT 30
    #--Update/Patch step----------
    #--Configure step-------------
    #--Build step-----------------
    BUILD_IN_SOURCE 1
    #--Install step---------------
    INSTALL_COMMAND ""
    #--Output logging-------------
    LOG_DOWNLOAD ON
  )
  ExternalProject_Get_Property(sqlite3-lib source_dir)
  set(SQLITE3_INCLUDE_DIR ${source_dir}/src CACHE INTERNAL "Path to include folder for SQLite3")
  set(SQLITE3_LIBRARY_DIR ${source_dir} CACHE INTERNAL "Path to library folder for SQLite3")
  set(SQLITE3_LIBRARY "sqlite3-static")
endif(NOT SQLITE3_FOUND)
link_directories(${SQLITE3_LIBRARY_DIR})

if(NOT APPLE)
  include_directories(SYSTEM AFTER "${SQLITE3_INCLUDE_DIR}")
else(NOT APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SQLITE3_INCLUDE_DIR}\"")
endif(NOT APPLE)

# -------------------------------

if(NOT BUILD_DEPENDENCIES)
  find_package(SQLiteCpp 1003001)
endif(NOT BUILD_DEPENDENCIES)

if(NOT SQLITECPP_FOUND)
  message(STATUS "SQLiteCpp will be downloaded when ${CMAKE_PROJECT_NAME} is built")
  ExternalProject_Add(sqlitecpp-lib
    PREFIX ${EXTERNAL_PATH}/SQLiteCpp
    #--Download step--------------
    URL https://github.com/SRombauts/SQLiteCpp/archive/master.zip
    TIMEOUT 30
    #--Update/Patch step----------
    #--Configure step-------------
    #--Build step-----------------
    BUILD_IN_SOURCE 1
    #--Install step---------------
    INSTALL_COMMAND ""
    #--Output logging-------------
    LOG_DOWNLOAD ON
  )
  ExternalProject_Get_Property(sqlitecpp-lib source_dir)
  set(SQLITECPP_INCLUDE_DIRS ${source_dir}/include
    CACHE INTERNAL "Path to include folder for SQLiteCpp")
  set(SQLITECPP_LIBRARY_DIR ${source_dir} CACHE INTERNAL "Path to library folder for SQLiteCpp")
  set(SQLITECPP_LIBRARY "SQLiteCpp")
endif(NOT SQLITECPP_FOUND)

if(NOT APPLE)
  include_directories(SYSTEM AFTER "${SQLITECPP_INCLUDE_DIRS}")
else(NOT APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SQLITECPP_INCLUDE_DIRS}\"")
endif(NOT APPLE)
link_directories(${SQLITECPP_LIBRARY_DIR})

# -------------------------------
