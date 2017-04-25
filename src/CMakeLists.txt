
INCLUDE_DIRECTORIES( ${MOAB_INCLUDE_DIRS} )

# ADD_EXECUTABLE(test_model test_model.cpp gen_mb_funcs.cpp ray_fire.cpp)
ADD_EXECUTABLE(sweep sweep.cpp hv_mesh_gen.cpp gen_mb_funcs.cpp ray_fire.cpp)
ADD_EXECUTABLE(write_obbs write_obbs.cpp obbhexwriter.cpp gen_mb_funcs.cpp)
ADD_EXECUTABLE(count_leaves count_leaves.cpp leafcounter.cpp gen_mb_funcs.cpp)
ADD_EXECUTABLE(hv_cube hv_cube.cpp hv_mesh_gen.cpp gen_mb_funcs.cpp ray_fire.cpp)
ADD_EXECUTABLE(traverse_ray traverse_ray.cpp ray_writer.hpp ray_writer.cpp)
ADD_EXECUTABLE(tag_tris_w_norm tag_tris_w_norm.cpp)

# TARGET_LINK_LIBRARIES(test_model ${DAGMC_LIBRARIES})
TARGET_LINK_LIBRARIES(sweep ${DAGMC_LIBRARIES} ${MPI_LIBRARIES})
TARGET_LINK_LIBRARIES(write_obbs ${DAGMC_LIBRARIES})
TARGET_LINK_LIBRARIES(count_leaves ${DAGMC_LIBRARIES})
TARGET_LINK_LIBRARIES(hv_cube ${DAGMC_LIBRARIES})
TARGET_LINK_LIBRARIES(traverse_ray ${DAGMC_LIBRARIES})
TARGET_LINK_LIBRARIES(tag_tris_w_norm ${DAGMC_LIBRARIES})

INSTALL(TARGETS     sweep
		    write_obbs
		    count_leaves
		    hv_cube
		    traverse_ray
		    tag_tris_w_norm
	DESTINATION ${INSTALL_BIN_DIR})
