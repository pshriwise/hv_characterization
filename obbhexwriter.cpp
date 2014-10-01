
#include "obbhexwriter.hpp"
#include "moab/Types.hpp"
#include "meshkit/MKCore.hpp"

using namespace moab;

OBBHexWriter::OBBHexWriter( OrientedBoxTreeTool *tool_ptr, Interface* interface_ptr )
  : tool(tool_ptr), mbi2(interface_ptr) {}

OBBHexWriter::~OBBHexWriter() {};

// nothing special about the leaves for this op, do nothing
ErrorCode OBBHexWriter::leaf( EntityHandle node ) { return MB_SUCCESS; };

// a visit to a node, will create an OrientedBox object for that node, get a hex from that box and tag it with 
// an integer tag representing it's depth in the tree, and add it to the list of *hexes* for this tree
ErrorCode OBBHexWriter::visit( EntityHandle node,
			   int depth, 
			   bool& descend)
{
  
 ErrorCode rval;
  
  //create a new tag for the tree depth of the obb
  std::string depth_tag_name = "TREE_DEPTH";
  rval = mbi2->tag_get_handle( depth_tag_name.c_str(),  1, MB_TYPE_INTEGER, depth_tag, 
			       MB_TAG_DENSE|MB_TAG_CREAT);
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 


  //get the triangles contained by this node
  std::vector<EntityHandle> tris; 
  
  rval = tool->get_moab_instance()->get_entities_by_type( node, MBTRI, tris);
  assert(MB_SUCCESS == rval );
  if ( MB_SUCCESS != rval ) return rval; 
  
  rval = transfer_tri_inst( tool->get_moab_instance(), mbi2, tris);
  assert( MB_SUCCESS == rval );

  OrientedBox box;
  rval = tool->box( node, box );
  assert(MB_SUCCESS == rval );
  if ( MB_SUCCESS != rval ) return rval; 
  
  rval = box.make_hex( new_hex, mbi2 );
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 
  
  void *ptr = &depth;
  rval = mbi2-> tag_set_data( depth_tag, &new_hex, 1, ptr);
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 
  
  //add the data to the class attributes
  tri_map[new_hex] = tris;
  hexes.push_back(new_hex);
  
  descend = true;    
  
  return MB_SUCCESS;
}


ErrorCode OBBHexWriter::transfer_tri_inst( Interface* orig_inst, Interface* new_inst, std::vector<EntityHandle> &tris)
{

  ErrorCode result;



  std::vector<EntityHandle> new_tris;


  for(std::vector<EntityHandle>::iterator i = tris.begin();
      i != tris.end(); i++)
    {

      std::vector<EntityHandle> old_verts;
      std::vector<EntityHandle> new_verts;
      //get the vertices of the original triangles
      result = orig_inst->get_adjacencies( &(*i), 1, 0, false, old_verts );
      assert( MB_SUCCESS == result );

      assert( old_verts.size() == 3 );

      //make new verts using the old coords (in the new instance)
      
      for(std::vector<EntityHandle>::iterator j = old_verts.begin();
	  j != old_verts.end(); j++)
	{

	  double *coords;
	  result = orig_inst->get_coords( &(*j), 1, coords);
	  assert( MB_SUCCESS == result);

	  EntityHandle new_vert;
	  result = new_inst->create_vertex( coords, new_vert );
	  assert( MB_SUCCESS == result);
	  
	  new_verts.push_back(new_vert);	  

	}

      EntityHandle new_tri; 
      assert( 3 == new_verts.size() );
      result = new_inst->create_element( MBTRI, &(new_verts[0]), 3, new_tri );
      assert( MB_SUCCESS == result);

      new_tris.push_back(new_tri);


      old_verts.clear();
      new_verts.clear();
      
    }

  tris.clear();
  tris = new_tris;

  return MB_SUCCESS;

}
//this function will go through the hexes found for the tree tool of this class
//and write them to file based on their depth in the tree
// filename formats: base_filename + "_" + depth_value + ".h5m"
ErrorCode OBBHexWriter::write_to_files( std::string base_filename )
{
  
  ErrorCode rval; 
  
  EntityHandle temp_set; 
  rval = mbi2->create_meshset( MESHSET_SET, temp_set);
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 
  
  int curr_depth = 0;
  while( !hexes.empty() )
    {
      
      std::vector<EntityHandle> to_write;
      for(std::vector<EntityHandle>::iterator i = hexes.begin(); 
	  i != hexes.end(); i++)
	{
	  
	  EntityHandle this_hex = *i;
	  int hex_depth;
	  void *data = &hex_depth;
	  
	  rval = mbi2->tag_get_data( depth_tag, &this_hex, 1, data);
	  assert( MB_SUCCESS == rval );
	  if( MB_SUCCESS != rval ) return rval; 
	  
	  if( hex_depth == curr_depth )
	    {
	      //add this entity to out list of ents to write
	      to_write.push_back(this_hex);
	      //remove it from the class list of hexes
	      hexes.erase(i);
	      
	      //reset iterator
	      i = hexes.begin()-1;
	    }
	  
	}
      
      //we should now have a list of vectors to write
      if( !to_write.empty() )
	{
	  //add them to the temp_meshset
	  rval = mbi2->add_entities( temp_set, &(to_write[0]), to_write.size() );
	  assert( MB_SUCCESS == rval );
	  if( MB_SUCCESS != rval ) return rval; 
	  
	  //write this meshset to file, appending the current depth then the file suffix
	  
	  std::ostringstream filename;
	  filename << base_filename << "_";
	  filename << curr_depth << ".vtk";


	  rval = mbi2->write_mesh( &(filename.str()[0]) , &temp_set, 1 );
	  assert( MB_SUCCESS == rval );
	  if( MB_SUCCESS != rval ) return rval; 
	  
	  //clear the meshset out and fill it with triangles
	  rval = mbi2->clear_meshset( &temp_set, 1);
	  assert( MB_SUCCESS == rval );

	  //fill w/ triangles
	  for(std::vector<EntityHandle>::iterator i = to_write.begin();
	      i != to_write.end(); i++)
	    {
	      if( 0 != tri_map[*i].size()) 
		{
	      rval = mbi2->add_entities( temp_set, &(tri_map[*i][0]), tri_map[*i].size());
	      assert( MB_SUCCESS == rval );
	      if( MB_SUCCESS != rval ) return rval; 
		}
	      else continue;

	    }

	  std::ostringstream filename1;
	  filename1 << base_filename << "_tris_";
	  filename1 << curr_depth << ".vtk";

	  //check to make sure there's something in the meshset?
	  std::vector<EntityHandle> test_for_ents;
	  rval = mbi2->get_entities_by_type( temp_set,  MBTRI, test_for_ents);
	  assert( MB_SUCCESS == rval );
	  if( 0 != test_for_ents.size() )
	    {
	      rval = mbi2->write_mesh( &(filename1.str()[0]) , &temp_set, 1 );
	      assert( MB_SUCCESS == rval );
	      if( MB_SUCCESS != rval ) return rval; 
	    }


	}
      //now clear out the to_write vector and the meshset
      to_write.clear();
      
      rval = mbi2->clear_meshset( &temp_set, 1);
      assert( MB_SUCCESS == rval );
      if( MB_SUCCESS != rval ) return rval; 
      
      // after writing the boxes for this depth of the tree, move down in depth
      curr_depth++;
      
    }
  return MB_SUCCESS;
  
}//end write_files


