
#include <iostream>
#include <sstream>
#include <iomanip> 
#include <limits> 
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "gen.hpp"

using namespace MeshKit;

MKCore *mk;

int main(int argc, char **argv)
{

  mk = new MKCore();
  
  //Load the mesh file
  mk->load_mesh("cube.h5m");

  //Get all of the surface ModelEnts
  MEntVector surfs;
  mk->get_entities_by_dimension(2, surfs);
  
  std::cout << "There are " << surfs.size() << " surfaces in this model." << std::endl;

moab::EntityHandle hv_surf;

  //Now to find one of the surfaces that is constant in z (for convenience)
for( unsigned int i = 0 ; i < surfs.size() ; i++)
    {
      //get the meshset handle for this surface
      moab::EntityHandle sh = surfs[i]->mesh_handle();

      //get the triangles in this meshset
      std::vector<moab::EntityHandle> tris;
      mk->moab_instance()->get_entities_by_type( sh, MBTRI, tris);
      std::cout << "There are " << tris.size() << " triangles in this surface." << std::endl;



//setup a check vector
MBCartVect check;
check [0] = 0; check [1] = 0; check[2] = 1;
//get the normal for the first triangle on each surface

//get the verts for each triangle and their coordinates
std::vector<moab::EntityHandle> verts;
mk->moab_instance()->get_connectivity( &(tris[0]), 1, verts);

MBCartVect coords[3];

mk->moab_instance()->get_coords( &(verts[0]), 3, coords[0].array() );

MBCartVect tri_norm;
gen::triangle_normal( coords[0], coords[1], coords[2], tri_norm);

std::cout << "The normal vector for this surface of the cube is " << tri_norm << std::endl;
      

if( tri_norm == check){
hv_surf = sh;
break;
}
     
    }
  
  
  return 0;

}



