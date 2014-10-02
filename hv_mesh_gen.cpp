
#include <iostream>
#include <ctime>
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
#include "DagMC.hpp"
#include "moab/OrientedBoxTreeTool.hpp"
#include "hv_mesh_gen.hpp"

//timing includes
#include <fstream>
#include <fcntl.h>
#include <unistd.h>
#include <fcntl.h>
#include <cstdlib>

using namespace MeshKit;

MKCore *mk = new MKCore();

moab::ErrorCode test_hv_cube_mesh( double A_f, double valence, double &ray_fire_time )
{

	  prep_mesh( A_f, valence);

	  ////////////// START OF MOAB STUFF \\\\\\\\\\\\\\\\\\\\\\\	\
	  
	  //now we'll try to load this mesh-file into a dagmc instance
	  moab::DagMC *dag = moab::DagMC::instance();

	  moab::ErrorCode result;
	  //try loading the file 
	  result = dag->load_file( "cube_mod.h5m" );
	  if( MB_SUCCESS != result) return MB_FAILURE;

	  //generate the OBB tree
	  result = dag->init_OBBTree();
	  if( MB_SUCCESS != result) return MB_FAILURE;

	  //get all of the volumes in the dagmc instance
	  moab::Range vols;
	  result = get_volumes( dag->moab_instance(), vols);
	  if( MB_SUCCESS != result) return MB_FAILURE; 
	  
	  //analyze mesh here
	  double avg_fire_time;
	  CartVect source;
	  source[0] = 0; source [1] = 0; source [2] = 0;
	  //call into the new functions for firing random rays and get the avg time
	  fire_rand_rays( dag, vols[0], 100000, avg_fire_time, source);
	  //write time to data file
	  
	  std::cout << "The average fire time for this mesh was: " << avg_fire_time << "s" << std::endl;
	  ray_fire_time = avg_fire_time;
	  
	  dag->moab_instance()->delete_mesh();

	  return MB_SUCCESS;

}

void prep_mesh( double A_f, int valence )
{

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
  refacet_surface( hv_surf, A_f, valence );

  mk->save_mesh("cube_mod.h5m");
  mk->delete_all();

}

void refacet_surface( moab::EntityHandle surf, double A_f, int valence )
{

  //remove all 2-D entities on the surfce
  tear_down_surface( surf );

 
  //now its time to create an empty middle box using the remaining surface verts
  std::vector<moab::EntityHandle> box;
  generate_box_space( surf, A_f, box );

  make_hv_region( surf, box, valence );

}

// returns the verts of an empty box, centered on the origin which is surrounded by triangles
void generate_box_space( moab::EntityHandle surf, double A_f, std::vector<moab::EntityHandle> &box_verts )
{

  std::vector<moab::EntityHandle> corners;
  mk->moab_instance()->get_entities_by_type( surf, MBVERTEX, corners );

  double surface_area = polygon_area( corners );
  //going to start taking advantage of knowing the geometry here...
  double surface_side = sqrt(surface_area);
  double cube_area = 6*surface_area;

  //now create the vertices for the new regions
  
  //based on surface area, get the length of one of the center-square sides
  double hv_area = A_f*cube_area; 
  if ( hv_area >= surface_area ) std::cout << "ERROR: Area fraction must be less than 1/6 for now." << std::endl;
  assert( hv_area < surface_area );
  double hv_side = sqrt(hv_area);

  //get the move distance for the given area. 
  double box_bump_dist = 0.5 * (surface_side - hv_side);
  double dam_bump_short = 0.25 * (surface_side - hv_side);
  double dam_bump_long = 0.5 * surface_side;

  //vectors for the dam nodes and the inner vert points
  
  std::vector<moab::EntityHandle> 
    nd, /*north dam tri verts*/ 
    sd, /*south dam tri verts*/ 
    ed, /*east dam tri verts*/ 
    wd; /*west dam tri verts*/ 
  moab::EntityHandle 
    ndv, /*north dam vertex*/ 
    sdv, /*south dam vertex*/ 
    edv, /*east dam vertex*/ 
    wdv; /*west dam vertex*/

  std::vector<moab::EntityHandle> box; /* inner box vert list */


  //loop over the original verts (corners) and create the needed vertices
  for(std::vector<moab::EntityHandle>::iterator i=corners.begin(); i!=corners.end(); i++)
    {
      //first get the coordinates of this corner
      MBCartVect coords;
      mk->moab_instance()->get_coords( &(*i), 1, coords.array() );
      
      //need three cartesian vectors telling us where to put the verts, two for the dam points, one for the inner point.
      MBCartVect to_ewdam, to_nsdam, to_box;
      to_ewdam[2] = 0; to_nsdam[2] = 0; to_box[2] = 0; //only moving on the x-y plane here
      
      //always want to create this vert
      //using the fact that we're centered on the origin here...
      //x-points
      if( coords[0] < 0 ){ 
	to_box[0] = box_bump_dist; }
      else  { 
	to_box[0] = -box_bump_dist; }
      //y-points
      if( coords[1] < 0 ){ 
	to_box[1]=box_bump_dist; }
      else  { 
	to_box[1] = -box_bump_dist; }

      //create the inner point
      moab::EntityHandle box_vert;
      mk->moab_instance()->create_vertex( (coords+to_box).array(), box_vert);
      box.push_back(box_vert);
      

      //figure out which dams we're adjacent to
      std::vector<moab::EntityHandle> *nsdam = &nd, *ewdam = &ed;
      moab::EntityHandle *nsdam_vert = &ndv, *ewdam_vert = &edv;

      //set the north-south dam based on our y location
      if( coords[1] < 0) { nsdam_vert = &sdv; nsdam = &sd; }
      else { nsdam_vert = &ndv; nsdam = &nd; }

      //set the east-west dam based on our x location
      if( coords[0] < 0) { ewdam_vert = &wdv; ewdam = &wd; }
      else { ewdam_vert = &edv; ewdam = &ed; }

      //////NORTH-SOUTH DAM\\\\\\\

      //if the dams already have info from another corner, no need to create it's point,
      //just add this corner to the dam triangle
      if( 0 != nsdam->size() && NULL != nsdam_vert ) { nsdam->push_back(*i); }
      // if the dam has no info, then we'll create it
      else {
	
	//set the dam coordinates
	if ( coords[0] < 0 ) to_nsdam[0] = dam_bump_long;
	else to_nsdam[0] = -dam_bump_long;
	if ( coords[1] < 0 ) to_nsdam[1] = dam_bump_short;
	else to_nsdam[1] = -dam_bump_short;

	//create the dam vertex 
	mk->moab_instance()->create_vertex( (coords+to_nsdam).array(), *nsdam_vert);
	
	//now add this corner and the dam vert to the dam verts list
	nsdam->push_back(*nsdam_vert); nsdam->push_back(*i);

      } 
	
      //////EAST-WEST DAM\\\\\\\

      //if the dams already have info from another corner, no need to create it's point,
      //just add this corner to the dam triangle
      if( 0 != ewdam->size() && NULL != ewdam_vert ) { ewdam->push_back(*i); }
      // if the dam has no info, then we'll create it
      else {
	
	//set the dam coordinates
	if ( coords[0] < 0 ) to_ewdam[0] = dam_bump_short;
	else to_ewdam[0] = -dam_bump_short;
	if ( coords[1] < 0 ) to_ewdam[1] = dam_bump_long;
	else to_ewdam[1] = -dam_bump_long;

	//create the dam vertex 
	mk->moab_instance()->create_vertex( (coords+to_ewdam).array(), *ewdam_vert);
	
	//now add this corner and the dam vert to the dam verts list
	ewdam->push_back(*ewdam_vert); ewdam->push_back(*i);

      } 
	
      //dam info should now all be set

      //there are always two triangles to create that connect the dam
      //to the corner, we'll do that now
      moab::EntityHandle tri_verts[4] = { *ewdam_vert, *i, box_vert, *nsdam_vert  };
      std::vector<moab::EntityHandle> tris(2);
      mk->moab_instance()->create_element( MBTRI, &tri_verts[0], 3, tris[0] );
      mk->moab_instance()->create_element( MBTRI, &tri_verts[1], 3, tris[1] );
      
      //now we'll add these to the set
      mk->moab_instance()->add_entities( surf, &tri_verts[0], 4 );
      mk->moab_instance()->add_entities( surf, &(tris[0]), tris.size() );
    }

  //now we should have all of the info needed to create our dam triangles
  std::vector<moab::EntityHandle> dam_tris(4);
  assert( nd.size()==3 );
  mk->moab_instance()->create_element( MBTRI, &(nd[0]), nd.size(), dam_tris[0]);
  mk->moab_instance()->create_element( MBTRI, &(sd[0]), sd.size(), dam_tris[1]);
  mk->moab_instance()->create_element( MBTRI, &(ed[0]), ed.size(), dam_tris[2]);
  mk->moab_instance()->create_element( MBTRI, &(wd[0]), wd.size(), dam_tris[3]);

  //add these to the surface
  mk->moab_instance()->add_entities( surf, &(dam_tris[0]), dam_tris.size() );

  
  //re-order the M verts
  std::vector<moab::EntityHandle> ordered_box(4);
  for(unsigned int i = 0; i < box.size() ; i++)
    {

      MBCartVect v_coords;
      mk->moab_instance()->get_coords( &(box[i]), 1, v_coords.array() );
      double x = v_coords[0]; double y = v_coords[1];
      // again, our surface is cenetered on the origin
      if ( x < 0 && y < 0) ordered_box[0] = box[i];
      else if ( x < 0 && y > 0) ordered_box[1] = box[i];
      else if ( x > 0 && y > 0) ordered_box[2] = box[i];
      else if ( x > 0 && y < 0) ordered_box[3] = box[i];
    }

  box_verts = ordered_box;
  //now that the middle box is ordered we can add the triangles that 
  //connect the dams to to the box

  box_verts.push_back(box_verts[0]); // psuedo-loop to make creating these tris easier
  moab::EntityHandle all_dam_verts[4] = { wdv, ndv, edv, sdv };
  std::vector<moab::EntityHandle> last_few_tris(4);

  //create triangle that connect the box to the dams
  for(unsigned int i = 0; i < last_few_tris.size(); i++)
    {
      moab::EntityHandle tri_verts[3] = { box_verts[i], all_dam_verts[i], box_verts[i+1] };
      mk->moab_instance()->create_element( MBTRI, &(tri_verts[0]), 3, last_few_tris[i]);
    }
      
  //now add these to the surface
  mk->moab_instance()->add_entities( surf, &last_few_tris[0], 4);

  box_verts = ordered_box;
  
} // end generate_box_space

//creates a watertight high-valence area with n-valencies per
void make_hv_region( moab::EntityHandle surf, std::vector<moab::EntityHandle> box_verts, int n ) 
{

  //start by making a paramaterization of the box diagonal
  
  //because we know the order of these verts, we can get the diagonal coordinates directly
  MBCartVect sw_coords, ne_coords;

  mk->moab_instance()->get_coords( &box_verts[0], 1, sw_coords.array() );
  mk->moab_instance()->get_coords( &box_verts[2], 1, ne_coords.array() );

  //now create the param vector
  MBCartVect sw_to_ne = ne_coords-sw_coords;

  // now, based on the valence asked for, get the new vertex coords long this line
  std::vector<moab::EntityHandle> diag_verts;
 

  for(unsigned int i = 1; i <= n-1; i++)
    {
      double u = double(i)/double(n);
      MBCartVect new_coords = sw_coords + u*sw_to_ne;
      moab::EntityHandle new_vert;
      mk->moab_instance()->create_vertex( new_coords.array(), new_vert);
      diag_verts.push_back(new_vert);
    }

  //add all of these new verts to the surface
  mk->moab_instance()->add_entities( surf, &(diag_verts[0]), diag_verts.size() );


  // put the box verts at the front and end
  diag_verts.insert( diag_verts.begin(), box_verts[0] );
  diag_verts.push_back( box_verts[2]);

 
  //now create the new triangles
  //each segment along the diagonal will will have two tris
  //one to the nw vert and one to the se vert

  moab::EntityHandle nw_vert = box_verts[1];
  moab::EntityHandle se_vert = box_verts[3];
  for(unsigned int i = 0; i < diag_verts.size()-1; i++)
    {
     
      moab::EntityHandle tri1, tri2;
      //creating new tris
      moab::EntityHandle tri1_verts[3] = { nw_vert, diag_verts[i], diag_verts[i+1] };
      moab::EntityHandle tri2_verts[3] = { se_vert, diag_verts[i], diag_verts[i+1] };
      mk->moab_instance()->create_element( MBTRI, &tri1_verts[0], 3, tri1);
      mk->moab_instance()->create_element( MBTRI, &tri2_verts[0], 3, tri2);

      //now add these to the surface set
      mk->moab_instance()->add_entities( surf, &tri1, 1);
      mk->moab_instance()->add_entities( surf, &tri2, 1);
    }
  
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

 moab::ErrorCode get_volumes( moab::Interface* mb, moab::Range &volumes)
 {

   moab::ErrorCode rval; 

   //get the GEOM_DIM tag
   moab::Tag geom_dim;
   rval = mb->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1,
		       MB_TYPE_INTEGER, geom_dim);
   assert( MB_SUCCESS == rval);
   if( MB_SUCCESS != rval) return rval; 

   // geom dimension of 3 should indicate a volume
   int dim = 3;
   void* ptr = &dim;
   
   //get all volume meshsets (should only be one)
   mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_dim, &ptr, 1, volumes);
   assert( MB_SUCCESS == rval);
   if( MB_SUCCESS != rval) return rval; 

   return MB_SUCCESS;
 }


moab::ErrorCode write_obb_mesh( moab::DagMC *dag, moab::EntityHandle vol, std::string& base_filename) 
{

  moab::ErrorCode rval; 

  moab::EntityHandle root;

  rval = dag->get_root( vol, root );
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 

  moab::OrientedBoxTreeTool *obbtool = dag->obb_tree();

  rval = obbtool->stats( root, std::cout );
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 

  //make a new moab core for the box hexes
  moab::Core mbi2;

  OBBHexWriter hw( obbtool, &mbi2 );

  moab::OrientedBoxTreeTool::TrvStats tree_stats;

  rval = obbtool->preorder_traverse( root, hw, &tree_stats );
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 

  rval = hw.write_to_files( base_filename );
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 


}
  


void fire_rand_rays( moab::DagMC *dagi, moab::EntityHandle vol, int num_rand_rays, double &avg_fire_time, moab::CartVect ray_source)
{

  //always start with the same seed
  srand(randseed);

  moab::CartVect xyz, uvw;

  int random_rays_missed = 0; 

  moab::EntityHandle dum;
  double dist;
  moab::OrientedBoxTreeTool::TrvStats trv_stats;

  //start timer
  std::clock_t start_time, end_time;

  start_time = std::clock();
  for (int j = 0; j < num_random_rays; j++) {
    RNDVEC(uvw, location_az);
    
    xyz = uvw * source_rad + ray_source;
    if (source_rad >= 0.0) {
      RNDVEC(uvw, direction_az);
    }
    
#ifdef DEBUG
    std::cout << "x,y,z,u,v,w,u^2 + v^2 + w^2 = " << xyz 
	      << " " << uvw << " " << uvw%uvw << std::endl;
    uavg += uvw[0]; vavg += uvw[1]; wavg += uvw[2];
#endif
    // added ray orientation
    dagi->ray_fire(vol, xyz.array(), uvw.array(), dum, dist, NULL, 0, 1, &trv_stats );
    
    if( dum == 0){ random_rays_missed++; }
    
  }
  //end timer 
  end_time = std::clock();
  
  if( random_rays_missed ){
    std::cout << "Warning: " << random_rays_missed << " random rays did not hit the target volume" << std::endl;
  }
  
  avg_fire_time = (end_time - start_time) / (double)(CLOCKS_PER_SEC * 100000); 
  std::cout << "Average ray fire time: " << avg_fire_time << "s" << std::endl;
}


