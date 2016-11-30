

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/CartVect.hpp"
#include "moab/ProgOptions.hpp"
#include <iostream>
#include <assert.h>
#include <string>

int main(int argc, char** argv)
{

  ProgOptions po("A program for tagging MOAB triangles with their normals.");

  std::string filename;
  po.addRequiredArg<std::string>("filename", "MOAB mesh file with triangles to tag.", &filename);

  po.parseCommandLine(argc,argv);

  moab::Interface* mb = new moab::Core();
  moab::ErrorCode result;
  
  result = mb->load_file(filename.c_str());
  MB_CHK_SET_ERR(result,"MOAB FAILURE");  
  moab::Range tris;
  result = mb->get_entities_by_type(0, moab::MBTRI, tris);
  MB_CHK_SET_ERR(result,"MOAB FAILURE");  std::cout << "Found " << tris.size() << " triangles in this model." << std::endl;
  moab::Range::iterator i;

  moab::Tag norm_tag;
  result = mb->tag_get_handle("TRI_NORM",3,moab::MB_TYPE_DOUBLE, norm_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR(result,"MOAB FAILURE");

  moab::EntityHandle vect_meshset;
  result = mb->create_meshset(0,vect_meshset);
  MB_CHK_SET_ERR(result,"Could not create vector meshset.");
  
  for(i = tris.begin(); i != tris.end(); i++)
    {
      
      moab::EntityHandle tri = *i;
      std::vector<moab::EntityHandle> verts;
      result = mb->get_connectivity(&tri, 1, verts);
      MB_CHK_SET_ERR(result,"MOAB FAILURE");
      moab::CartVect coords[3];
      assert(3 == verts.size());

      result = mb->get_coords( &(verts[0]), 1, coords[0].array());
      MB_CHK_SET_ERR(result,"MOAB FAILURE");
      result = mb->get_coords( &(verts[1]), 1, coords[1].array());
      MB_CHK_SET_ERR(result,"MOAB FAILURE");
      result = mb->get_coords( &(verts[2]), 1, coords[2].array());
      MB_CHK_SET_ERR(result,"MOAB FAILURE");
      //create a vertex at the triangle's center
      moab::CartVect tri_center_coords = (coords[0] + coords[1] + coords[2] )/3;

      moab::EntityHandle tri_center;
      result = mb->create_vertex( tri_center_coords.array(), tri_center);
      MB_CHK_SET_ERR(result,"MOAB FAILURE");
      result = mb->add_entities(vect_meshset, &tri_center,1);
      MB_CHK_SET_ERR(result,"Could not add tri center vert to vect meshset.");
      moab::CartVect v1 = coords[1] - coords[0];
      moab::CartVect v2 = coords[2] - coords[0];
      moab::CartVect norm = v1*v2;
      norm.normalize();

      void *ptr = norm.array();
      result = mb->tag_set_data(norm_tag, &tri_center, 1, ptr);
      MB_CHK_SET_ERR(result,"MOAB FAILURE");
    }

  std::cout << "Writing file..." << std::endl;
  std::string basename = filename.substr(0, filename.find_last_of("."));
  std::string outfile = basename + "_norms.h5m";
  result = mb->write_mesh(outfile.c_str(),&vect_meshset,1);
  MB_CHK_SET_ERR(result,"MOAB FAILURE");

  return 0;
}
