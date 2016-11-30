
INCLUDE_DIRECTORIES( ${MOAB_INCLUDE_DIRS} )

ADD_EXECUTABLE(traverse_ray traverse_ray.cpp ray_writer.hpp ray_writer.cpp)

TARGET_LINK_LIBRARIES( traverse_ray ${DAGMC_LIBRARIES} )

INSTALL(TARGETS traverse_ray DESTINATION ${CMAKE_INSTALL_PREFIX})

