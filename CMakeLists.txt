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

FIND_PACKAGE( MPI REQUIRED )

FIND_PACKAGE ( MOAB REQUIRED ) 
MESSAGE ( STATUS "Found MOAB Version: " ${MOAB_VERSION} )
  MESSAGE ( STATUS "MOAB_LIBRARIES is " ${MOAB_LIBRARIES})
  MESSAGE ( STATUS "MOAB_INCLUDE_DIRS is " ${MOAB_INCLUDE_DIRS})
INCLUDE_DIRECTORIES( ${MOAB_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})

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

INCLUDE_DIRECTORIES( ${MOAB_INCLUDE_DIRS} )

ADD_EXECUTABLE(test_model test_model.cpp gen_mb_funcs.cpp ray_fire.cpp)
ADD_EXECUTABLE(sweep sweep.cpp hv_mesh_gen.cpp gen_mb_funcs.cpp ray_fire.cpp)
ADD_EXECUTABLE(write_obbs write_obbs.cpp obbhexwriter.cpp gen_mb_funcs.cpp)
ADD_EXECUTABLE(count_leaves count_leaves.cpp leafcounter.cpp gen_mb_funcs.cpp)
ADD_EXECUTABLE(hv_cube hv_cube.cpp hv_mesh_gen.cpp gen_mb_funcs.cpp ray_fire.cpp)

TARGET_LINK_LIBRARIES(test_model ${DAGMC_LIBRARIES})
TARGET_LINK_LIBRARIES(sweep ${DAGMC_LIBRARIES} ${MPI_LIBRARIES} )
TARGET_LINK_LIBRARIES(write_obbs ${DAGMC_LIBRARIES})
TARGET_LINK_LIBRARIES(count_leaves ${DAGMC_LIBRARIES})
TARGET_LINK_LIBRARIES(hv_cube ${DAGMC_LIBRARIES})


if(MPI_COMPILE_FLAGS)
  set_target_properties(hello PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(hello PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

ADD_SUBDIRECTORY(single_ray_traversal)

INSTALL(TARGETS test_model sweep write_obbs count_leaves hv_cube DESTINATION ${CMAKE_INSTALL_PREFIX})