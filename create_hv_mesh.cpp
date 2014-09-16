
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


void refacet_surface( moab::EntityHandle surf, double A_f );
void get_hv_surf( MEntVector surfs, moab::EntityHandle &hv_surf);
void tear_down_surface( moab::EntityHandle surf );
double polygon_area( std::vector<moab::EntityHandle> verts);
//assumes verts are in the proper order
void add_box_to_surf( moab::EntityHandle surf, std::vector<moab::EntityHandle> verts);
void replace_curves( std::vector<moab::EntityHandle> orig_curves, std::vector<std::vector<moab::EntityHandle> > new_curve_verts);
  
int main(int argc, char **argv)
{

  //temp value of A_f, the fraction of the total surface area claimed
  //by high valencies
  double A_f = 0.4;

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

  //refacet the surface using the desired area fraction for the hv region
  refacet_surface( hv_surf, A_f );

  mk->save_mesh("cube_mod.h5m");
  return 0;

}

void refacet_surface( moab::EntityHandle surf, double A_f )
{

  //remove all 2-D entities on the surfce
  tear_down_surface( surf );

 
  //now its time to create new boxes using the remaining surface verts
  std::vector<moab::EntityHandle> orig_verts;
  mk->moab_instance()->get_entities_by_type( surf, MBVERTEX, orig_verts);

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
  /*
                    //  NEW SURFACE DESIGN  \\

                          North Curve
   //////////////////////////////////////////////////////////////
   //              //                        //                //
   //              //                        //                //
   //              //            U           //                //
   //              //                        //                //
   //              //                        //                //
 W ////////////////////////////////////////////////////////////// E
 e //              //                        //                // a
 s //              //                        //                // s
 t //              //                        //                // t
   //              //                        //                //
 C //      L       //           M            //       R        // C
 u //              //                        //                // u
 r //              //                        //                // r
 v //              //                        //                // v
 e //              //                        //                // e
   //////////////////////////////////////////////////////////////
   //              //                        //                //
   //              //                        //                //
   //              //           B            //                //
   //              //                        //                //
   //              //                        //                //
   //////////////////////////////////////////////////////////////
                           South Curve
  */

  // Box lists (the corner boxes are created with the new verts)  
  std::vector<moab::EntityHandle> L,R,M,T,B;
  // Curve Lists (this is where we'll place the new curve verts)
  std::vector<moab::EntityHandle> N,S,E,W;
  //start creating new verts and adding them to the appropriate lists

  //create the new areas on this surface
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

      // select the proper box to add to based on the corner point we're expanding from      
      std::vector<moab::EntityHandle> *y_box = &L, *x_box= &T; //dummy setting for the pointers
      if( coords[0] < 0 ) y_box = &L; else y_box = &R;
      if( coords[1] < 0 ) x_box = &B; else x_box = &T;

      // if there are no verts in the list, add xy then y, if verts exist, add y first instead
      if( y_box->size() != 0 ) { y_box->push_back(xy); y_box->push_back(y); }
      else { y_box->push_back(y); y_box->push_back(xy); }
      // same as above but for the x-adj box
      if( x_box->size() != 0 ) { x_box->push_back(xy); x_box->push_back(x); }
      else { x_box->push_back(x); x_box->push_back(xy); }
      

      // this is a similar process as above, but for creating the new curve data
      std::vector<moab::EntityHandle> *y_curve = &W, *x_curve = &S;
      if( coords[0] < 0 ) y_curve = &W; else y_curve = &E;
      if( coords[1] < 0 ) x_curve = &S; else x_curve = &N;

      if( y_curve->size() != 0 ) { y_curve->push_back(y); y_curve->push_back(*i); }
      else { y_curve->push_back(*i); y_curve->push_back(y); }

      if( x_curve->size() != 0 ) { x_curve->push_back(x); x_curve->push_back(*i); }
      else { x_curve->push_back(*i); x_curve->push_back(x); }
      
      //create a new quad here (for now)
      std::vector<moab::EntityHandle> quad_verts(4);
      quad_verts[0] =  x;
      quad_verts[1] =  xy;
      quad_verts[2] =  y;
      quad_verts[3] =  *i;

      moab::EntityHandle new_quad;
      add_box_to_surf( surf, quad_verts);

    }

  //create the L box
  add_box_to_surf( surf, L);
  //create the T box
  add_box_to_surf( surf, T);
  //create the R box
  add_box_to_surf( surf, R);
  //create the B box
  add_box_to_surf( surf, B);


  //re-order the M verts
  std::vector<moab::EntityHandle> new_M(4);
  for(unsigned int i = 0; i < M.size() ; i++)
    {

      MBCartVect v_coords;
      mk->moab_instance()->get_coords( &(M[i]), 1, v_coords.array() );
      double x = v_coords[0]; double y = v_coords[1];

      if ( x < 0 && y < 0) new_M[0] = M[i];
      else if ( x < 0 && y > 0) new_M[1] = M[i];
      else if ( x > 0 && y > 0) new_M[2] = M[i];
      else if ( x > 0 && y < 0) new_M[3] = M[i];
    }

  M = new_M;
  
  //create the middle box
  add_box_to_surf( surf, M);

  //replace the old curves with the new ones
  std::vector<moab::EntityHandle> curves;
  mk->moab_instance()->get_child_meshsets( surf, curves);
  std::vector<std::vector<moab::EntityHandle> > new_curve_verts;
  new_curve_verts.push_back(N);   new_curve_verts.push_back(S);
  new_curve_verts.push_back(E);   new_curve_verts.push_back(W);

  replace_curves( curves, new_curve_verts);

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
  verts.push_back(verts[0]); //verts.push_back(verts[2]);

  //now create the new triangles and add everything to the surface set
  mk->moab_instance()->create_element( MBTRI, &(verts[0]), 3, tris[0]);
  mk->moab_instance()->create_element( MBTRI, &(verts[2]), 3, tris[1]);

  mk->moab_instance()->add_entities( surf, &(verts[0]), 4);
  mk->moab_instance()->add_entities( surf, &(tris[0]), 2);

}

void replace_curves( std::vector<moab::EntityHandle> orig_curves, std::vector<std::vector<moab::EntityHandle> > new_curve_verts)
{
  
  int matches = 0;
  //for each original curve, get the verts and see if there is a match for the beginning and end points
  for(std::vector<moab::EntityHandle>::iterator i = orig_curves.begin(); i != orig_curves.end(); i++)
    {
      
      //get the curve verts
      std::vector<moab::EntityHandle> orig_curve_verts;
      mk->moab_instance()->get_entities_by_type( *i, MBVERTEX, orig_curve_verts);

      for(std::vector<std::vector<moab::EntityHandle> >::iterator j = new_curve_verts.begin();  j != new_curve_verts.end(); j++)
	{
	  //look for a match 
	  bool match = false;
	  std::vector<moab::EntityHandle> these_verts= *j;
	  if( orig_curve_verts.front() == these_verts.front() && orig_curve_verts.back() == these_verts.back() )
	    match = true;
	  //check the reverse sense in case the vert orders are just backwards
	  if( orig_curve_verts.back() == these_verts.front() && orig_curve_verts.front() == these_verts.back() )
	    {
	      std::reverse(these_verts.begin(), these_verts.end());
	      match = true;
	    }

	  //if we have a match we will clear-out the meshset
	  //and add the new verts to replace the previous curve
	  if(match)
	    {
	      matches++;
	      mk->moab_instance()->clear_meshset( &(*i), 1);
	      mk->moab_instance()->add_entities( *i, &(these_verts[0]), these_verts.size());
	      
	    }
	  else continue;
	} // new curve verts loop
    }// orig curves loop

  if(matches != 4) std::cout << "Warning: Not all curves were replaced successfully!" << std::endl;
}
