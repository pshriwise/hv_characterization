

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

  //check if the ray intersects the box
  descend = box.intersect_ray( ray_origin, ray_direction, tol, nonneg_ray_len, 
			       neg_ray_len );

  //if the ray does intersect the box
  if (descend) {
    //createa a hex from this box 
    EntityHandle hex;
    rval = box.make_hex( hex, MBI );
    MB_CHK_SET_ERR(rval, "Could not create hex from box.");
    // add it to the MeshSet we want to write out
    rval = MBI->add_entities(*writeSet, &hex, 1);
    MB_CHK_SET_ERR(rval, "Could not add box hex to output set.");
  }
  
};


ErrorCode RayTraversalWriter::write_output_file()
{
  //write out everything in the writeSet in vtk format
  ErrorCode rval = MBI->write_mesh("ray_traversal.vtk", writeSet, 1);
  MB_CHK_SET_ERR(rval, "Could not write the output mesh file.");
}
//NO DEFINITION FOR LEAF METHOD, WE HAVE NOTHING TO DO THERE
