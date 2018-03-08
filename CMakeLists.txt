####################################################
#
# CMakeLists.txt for top-level EMDAG repository
#
#####################################################

CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

ADD_DEFINITIONS(
    -std=c++11 # Or -std=c++0x
)

###########################
#
# FIND DEPENDENCIES
#
##########################
# Needed for FindPyne.cmake
SET ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake)

INCLUDE_DIRECTORIES("/usr/include/eigen3/")

FIND_PACKAGE( MPI REQUIRED )

FIND_PACKAGE ( MOAB REQUIRED ) 
MESSAGE ( STATUS "Found MOAB Version: " ${MOAB_VERSION} )
  MESSAGE ( STATUS "MOAB_LIBRARIES is " ${MOAB_LIBRARIES})
  MESSAGE ( STATUS "MOAB_INCLUDE_DIRS is " ${MOAB_INCLUDE_DIRS})
INCLUDE_DIRECTORIES( ${MOAB_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})

FIND_PACKAGE ( DAGMC REQUIRED ) 
MESSAGE ( STATUS "Found DAGMC Version: " ${DAGMC_VERSION} )
MESSAGE ( STATUS "DAGMC_LIBRARIES is " ${DAGMC_LIBRARIES})
MESSAGE ( STATUS "DAGMC_LIBRARY_DIRS is " ${DAGMC_LIBRARY_DIRS})
MESSAGE ( STATUS "DAGMC_INCLUDE_DIRS is " ${DAGMC_INCLUDE_DIRS})
LINK_DIRECTORIES(${DAGMC_LIBRARY_DIRS})
INCLUDE_DIRECTORIES( ${DAGMC_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})

IF(SIMD_BVH_BUILD)
  ADD_DEFINITIONS(-DSIMD_BVH)
ENDIF(SIMD_BVH_BUILD)

SET ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH})

###################################
#
# Install dirs 
#
###################################

SET( INSTALL_LIB_DIR     lib     CACHE PATH "Installation directory for libraries" )
SET( INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header" )
SET( INSTALL_BIN_DIR     bin     CACHE PATH "Installation directory for executables" )

###################################
#
# Set default install location if not specified
#
###################################

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX ../)
endif()


ADD_SUBDIRECTORY(src)
