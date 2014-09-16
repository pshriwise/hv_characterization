
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
//assumes verts are in the proper order
void add_box_to_surf( moab::EntityHandle surf, std::vector<moab::EntityHandle> verts);
  
int main(int argc, char **argv)
{


  //temp value of A_f, the fraction of the total surface area claimed
  //by high valencies
  double A_f = 0.2;

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
  //going to start taking advantage of knowing the geometry here...
  double surface_side = sqrt(surface_area);
  double cube_area = 6*surface_area;

  //now create the vertices for the new regions
  
  //based on surface area, get the length of one of the center-square sides
  double hv_area = A_f*surface_area; 
  assert( hv_area < surface_area );
  double hv_side = sqrt(hv_area);

  //get the move distance for the given area. 
  double bump_dist = 0.5*(surface_side - hv_side);

  std::vector<moab::EntityHandle> L,R,M,T,B;
  //start creating new verts and adding them to the appropriate lists

  //take the first vertex and determine where it is
  for(std::vector<moab::EntityHandle>::iterator i = orig_verts.begin(); i!=orig_verts.end(); i++)
    {
      MBCartVect coords;
      mk->moab_instance()->get_coords( &(*i), 1, coords.array() );

      //create new vertex vectors
      double x_dist, y_dist;
      if( coords[0] < 0 ) x_dist = bump_dist; else x_dist = -bump_dist;
      if( coords[1] < 0 ) y_dist = bump_dist; else y_dist = -bump_dist;

      MBCartVect x_only, y_only, x_n_y;
      x_only[0]= x_dist; x_only[1] = 0; x_only[2] = 0;
      y_only[0]= 0; y_only[1] = y_dist; y_only[2] = 0;
      x_n_y[0] = x_dist; x_n_y[1] = y_dist; x_n_y[2] = 0;

      //now create new verts at each of these points
      moab::EntityHandle x,y,xy;
      mk->moab_instance()->create_vertex( (coords+x_only).array(), x);
      mk->moab_instance()->create_vertex( (coords+y_only).array(), y);
      mk->moab_instance()->create_vertex( (coords+x_n_y).array(), xy);

      //add the xy vert into the MM set...
      M.push_back(xy);

      std::vector<moab::EntityHandle> *y_box = &L, *x_box=&T;
      // select the proper box to add to based on the corner point we're expanding from
      if( coords[0] < 0 ) y_box = &L; else y_box = &R;
      if( coords[1] < 0 ) x_box = &B; else x_box = &T;

      if( y_box->size() != 0 ) { y_box->push_back(xy); y_box->push_back(y); }
      else { y_box->push_back(y); y_box->push_back(xy); }

      if( x_box->size() != 0 ) { x_box->push_back(xy); x_box->push_back(x); }
      else { x_box->push_back(x); x_box->push_back(xy); }
      
      
      //create a new quad here (for now)
      moab::EntityHandle quad_verts[] = { x, xy, y, *i};
      moab::EntityHandle new_quad;
      mk->moab_instance()->create_element( MBQUAD, &(quad_verts[0]), 4, new_quad);

      //add this quad and verts to the surface set
      mk->moab_instance()->add_entities( hv_surf, &(quad_verts[0]), 4);
      mk->moab_instance()->add_entities( hv_surf, &new_quad, 1);
    }

  //create the L box
  add_box_to_surf( hv_surf, L);
  //create the T box
  add_box_to_surf( hv_surf, T);
  //create the R box
  add_box_to_surf( hv_surf, R);
  //create the B box
  add_box_to_surf( hv_surf, B);

  //create the MM quad 
  moab::EntityHandle M_quad;
  mk->moab_instance()->create_element( MBQUAD, &(M[0]) , 4, M_quad);
  mk->moab_instance()->add_entities( hv_surf, &(M[0]), 4);
  mk->moab_instance()->add_entities( hv_surf, &M_quad, 1);

  std::vector<moab::EntityHandle> *p = &M;

  int foo = 6; 
  int *bar;
  bar = &foo;
  
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

//assumes verts are in the proper order
void add_box_to_surf( moab::EntityHandle surf, std::vector<moab::EntityHandle> verts)
{
  
  //create some storage for the triangles
  std::vector<moab::EntityHandle> tris(2);
  //create an psuedo-wrap in the vector to make triangle creation easy
  verts.push_back(verts[3]); verts.push_back(verts[2]);

  //now create the new triangles and add everything to the surface set
  mk->moab_instance()->create_element( MBTRI, &(verts[0]), 3, tris[0]);
  mk->moab_instance()->create_element( MBTRI, &(verts[3]), 3, tris[1]);

  mk->moab_instance()->add_entities( surf, &(verts[0]), 4);
  mk->moab_instance()->add_entities( surf, &(tris[0]), 2);

}
