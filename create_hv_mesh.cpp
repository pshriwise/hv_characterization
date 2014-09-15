
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
#include "../src/tools/measure.hpp"

using namespace MeshKit;

MKCore *mk;

void get_hv_surf( MEntVector surfs, moab::EntityHandle &hv_surf);
void tear_down_surface( moab::EntityHandle surf );
double polygon_area( std::vector<moab::EntityHandle> verts);

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
  moab::EntityHandle hv_surf;

  get_hv_surf( surfs, hv_surf );

  tear_down_surface( hv_surf );
 
  //now its time to create new boxes using the remaining surface verts

  std::vector<moab::EntityHandle> orig_verts;
  mk->moab_instance()->get_entities_by_type( hv_surf, MBVERTEX, orig_verts);

  double surface_area = polygon_area( orig_verts );

  

  mk->save_mesh("cube_mod.h5m");
  return 0;

}

void get_hv_surf( MEntVector surfs, moab::EntityHandle &hv_surf)
{

  for( unsigned int i = 0 ; i < surfs.size() ; i++)
    {
      //get the meshset handle for this surface
      moab::EntityHandle sh = surfs[i]->mesh_handle();
      
      //get the triangles in this meshset
      std::vector<moab::EntityHandle> tris;
      mk->moab_instance()->get_entities_by_type( sh, MBTRI, tris);
      
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
      
      if( tri_norm == check){
	hv_surf = sh;
	break;
      }

    } //end loop  
}


void tear_down_surface( moab::EntityHandle surf)
{

 //get the triangles for this surface and delete them
  std::vector<moab::EntityHandle> tris;
  mk->moab_instance()->get_entities_by_type( surf, MBTRI, tris);

  //remove these triangle from the meshset and destroy them
  std::cout << "Tearing down the surface..." << std::endl; 

  mk->moab_instance()->remove_entities( surf, &(tris[0]), tris.size());
  mk->moab_instance()->delete_entities( &(tris[0]), tris.size());

  std::vector<moab::EntityHandle> test_tris;
  mk->moab_instance()->get_entities_by_type( 0, MBTRI, test_tris, true);
  std::cout << "There are now " << test_tris.size() << " triangles in the model" << std::endl;

}

double polygon_area( std::vector<moab::EntityHandle> verts)
{
  unsigned int num_vertices = verts.size();

  //std::cout << "There are " << num_vertices << " vertices in this surface." << std::endl;

  //get the coordinates of the vertices
  double vert_coords[verts.size()*3];
  mk->moab_instance()->get_coords( &(verts[0]), verts.size(), &(vert_coords[0]) );

  const MBCartVect* coords = reinterpret_cast<const MBCartVect*>(vert_coords);  

  //Calculate the surface's area
  MBCartVect mid(0,0,0);
  for (int i = 0; i < num_vertices; ++i) mid += coords[i];
  mid /= num_vertices;

  double sum = 0.0;
  for (int i = 0; i < num_vertices; i++)
    {
      int j = (i+1)%(num_vertices);
      MBCartVect pnt1, pnt2;
      pnt1 = coords[i];
      pnt2 = coords[j];
      sum += ((mid - pnt1) * (mid - pnt2)).length();
    }
  double  poly_area = sum;

  //std::cout << "The total area of this polygon is: " << poly_area << std::endl;
  
  return poly_area;
}
