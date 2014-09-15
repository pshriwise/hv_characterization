
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

  //Now to find one of the surfaces that is constant in z (for convenience)
  for( MEntVector::iterator i = surfs.begin(); i != surfs.end() ; i++)
    {
      //get the meshset handle for this surface
      moab::EntityHandle sh = surfs[0]->mesh_handle();
      //get the triangles in this meshset
      std::vector<moab::EntityHandle> tris;
      mk->moab_instance()->get_entities_by_type( sh, MBTRI, tris);
      std::cout << "There are " << tris.size() << " triangles in this surface." << std::endl;


    }
  
  
  return 0;

}
