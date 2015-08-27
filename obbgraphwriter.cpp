
#include "obbgraphwriter.hpp"
#include "moab/Types.hpp"
#include "meshkit/MKCore.hpp"

using namespace moab;




OBBGraphWriter::OBBGraphWriter( OrientedBoxTreeTool *tool_ptr )
  : tool(tool_ptr)
{
    output.open("graph.dot");

}


OBBGraphWriter::~OBBGraphWriter() { output.close(); }

// no need for anything here
ErrorCode OBBGraphWriter::leaf( EntityHandle node ) { return MB_SUCCESS; }


//gathers the number of triangles underneath the node, its child relationships and adds them to the graph
ErrorCode OBBGraphWriter::visit( EntityHandle node,
			   int depth, 
			   bool& descend)
{

  //get the number of triangles under this node
  int num_tris = 0;
  
  Range all_children; 
  ErrorCode rval = tool->get_moab_instance()->get_child_meshsets( node, all_children, 0);
  MB_CHK_SET_ERR(rval, "Failed to get all children of the node.");

  //in case this is a leaf node
  all_children.insert(node);

  for( Range::iterator i = all_children.begin(); i != all_children.end(); i++)
    {
      int dum = 0;
      rval = tool->get_moab_instance()->get_number_entities_by_type( *i, MBTRI, dum);
      MB_CHK_SET_ERR(rval, "Failed to get the triangles under this node.");
      num_tris += dum;
    }

  output << node <<  " " << "[label = " << num_tris << "]" << std::endl;

  Range children; 
  rval = tool->get_moab_instance()->get_child_meshsets( node, children ); 
  MB_CHK_SET_ERR(rval, "Failed to get the children of this node." );

  if (0 != children.size() )
    {
      assert( 2 == children.size() );
      output << node << " -> " << children[0] << std::endl;
      output << node << " -> " << children[1] << std::endl;
    }
  descend = true;


  return MB_SUCCESS;
}


void OBBGraphWriter::open_graph() { output << "digraph {" << std::endl; }


void OBBGraphWriter::close_graph() { output<< std::endl << "}" << std::endl; }
