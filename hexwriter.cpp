
#include "hexwriter.hpp"
#include "moab/Types.hpp"
#include "meshkit/MKCore.hpp"

using namespace moab;

HexWriter::HexWriter( OrientedBoxTreeTool *tool_ptr, Interface* interface_ptr )
  : tool(tool_ptr), mbi2(interface_ptr) {}

HexWriter::~HexWriter() {};

// nothing special about the leaves for this op, do nothing
ErrorCode HexWriter::leaf( EntityHandle node ) { return MB_SUCCESS; };

// a visit to a node, will create an OrientedBox object for that node, get a hex from that box and tag it with 
// an integer tag representing it's depth in the tree, and add it to the list of *hexes* for this tree
ErrorCode HexWriter::visit( EntityHandle node,
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
  
  
  OrientedBox box;
  rval = tool->box( node, box );
  assert(MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 
  
  rval = box.make_hex( new_hex, mbi2 );
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 
  
  void *ptr = &depth;
  rval = mbi2-> tag_set_data( depth_tag, &new_hex, 1, ptr);
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 
  
  hexes.push_back(new_hex);
  
  descend = true;    
  
  return MB_SUCCESS;
}

//this function will go through the hexes found for the tree tool of this class
//and write them to file based on their depth in the tree
// filename formats: base_filename + "_" + depth_value + ".h5m"
ErrorCode HexWriter::write_to_files( std::string base_filename )
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


