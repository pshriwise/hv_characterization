

#include "ray_writer.hpp"

using namespace moab;

ErrorCode RayTraversalWriter::visit( EntityHandle node,
				     int depth,
				     bool& descend)
{
  //retrieve this node's bounding box
  OrientedBox box;
  ErrorCode rval = tool->box( node, box );
  assert(MB_SUCCESS == rval);
  if (MB_SUCCESS != rval)
    return rval;

  if (depth > max_depth) max_depth = depth;
  
  //check if the ray intersects the box
  descend = box.intersect_ray( ray_origin, ray_direction, tol, nonneg_ray_len, 
			       neg_ray_len );

  //if the ray does intersect the box
  if (descend) {
    //createa a hex from this box 
    EntityHandle hex;
    rval = box.make_hex(hex, MBI);
    MB_CHK_SET_ERR(rval, "Could not create hex from box.");
    //tag hex with it's depth in the tree
    void* ptr = &depth;
    rval = MBI->tag_set_data(depth_tag, &hex, 1, ptr);
    MB_CHK_SET_ERR(rval, "Could not tag hex with depth value.");
    // add it to the MeshSet we want to write out
    rval = MBI->add_entities(writeSet, &hex, 1);
    MB_CHK_SET_ERR(rval, "Could not add box hex to output set.");
  }
  
};

ErrorCode RayTraversalWriter::leaf( EntityHandle node ) { return MB_SUCCESS; };

ErrorCode RayTraversalWriter::write_single_output_file()
{
  //write out everything in the writeSet in vtk format
  ErrorCode rval = MBI->write_mesh("ray_traversal.vtk", &writeSet, 1);
  MB_CHK_SET_ERR(rval, "Could not write the output mesh file.");
  return rval;
}

ErrorCode RayTraversalWriter::write_ray_file()
{

  ErrorCode rval;
  //create ray_dir tag for vector representation of ray
  Tag ray_dir_tag;
  std::string ray_dir_name = "RAY_DIR";
  rval = MBI->tag_get_handle(ray_dir_name.c_str(), 3, MB_TYPE_DOUBLE, ray_dir_tag, MB_TAG_SPARSE|MB_TAG_CREAT);
  MB_CHK_ERR_CONT(rval);

  EntityHandle temp_set;
  rval = MBI->create_meshset(MESHSET_SET, temp_set);
  MB_CHK_SET_ERR(rval, "Could not create temporary meshset for output.");
  
  //create new verex for ray origin and tag with ray direction
  EntityHandle vert;
  rval = MBI->create_vertex(ray_origin.array(), vert);
  MB_CHK_ERR_CONT(rval);
  10*ray_direction;
  const void *ptr = ray_direction.array();
  rval = MBI->tag_set_data(ray_dir_tag, &vert, 1, ptr);
  MB_CHK_ERR_CONT(rval);
  ray_direction/10;
  //add this vert to the output set
  rval = MBI->add_entities(temp_set, &vert, 1);
  MB_CHK_ERR_CONT(rval);

  //write the ray to its own vtk file
  rval = MBI->write_mesh("RAY.vtk", &temp_set, 1);
  MB_CHK_SET_ERR(rval, "Could not write the ray vector file.");

  return rval;
  
};

ErrorCode RayTraversalWriter::write_vtk_database()
{

  std::string basename = "RAY_BOXES";
  ErrorCode rval;
  EntityHandle temp_set;
  rval = MBI->create_meshset(MESHSET_SET, temp_set);
  MB_CHK_SET_ERR(rval, "Could not create temporary meshset for output.");
  
  for(unsigned int i = 0; i <= max_depth; i++) {
    //ensure this meshset is empty
    rval = MBI->clear_meshset( &temp_set, 1);
    MB_CHK_SET_ERR(rval, "Could not clear out temporary meshset.");
      
    void *ptr[1] = {&i};
    Range hexes;
    //get all hexes with this tag value    
    rval = MBI->get_entities_by_type_and_tag(writeSet, MBHEX, &depth_tag, ptr, 1, hexes);
    MB_CHK_SET_ERR(rval, "Could not get all of the hexes.");
    //we should always find something
    assert(0 != hexes.size());
      
    //add this set of hexes to the meshset
    rval = MBI->add_entities(temp_set, hexes);
    MB_CHK_SET_ERR(rval, "Could not add hexes to temporary meshset.");

    //now write out this meshset with the appropriate filename
    std::string output_filename = basename+"_"+std::to_string(i)+".vtk";
    rval = MBI->write_mesh(output_filename.c_str(), &temp_set, 1);
    MB_CHK_SET_ERR(rval, "Could not write vtk datbase file.");
  }

  return rval;  
};
