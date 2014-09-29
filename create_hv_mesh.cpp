
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
#include "hexmaker.hpp"


//timing includes
#include <fstream>
#include <fcntl.h>
#include <unistd.h>
#include <fcntl.h>
#include <cstdlib>



#if !defined(_MSC_VER) && !defined(__MINGW32__)
#include <sys/resource.h>
#endif
#ifdef SOLARIS
extern "C" int getrusage(int, struct rusage *);
#ifndef RUSAGE_SELF
#include </usr/ucbinclude/sys/rusage.h>
#endif
#endif


using namespace MeshKit;

MKCore *mk;


////// Functions for generating the hv region \\\\\\\\\\

void prep_mesh( double A_f, int valence );
// creates new facets for the square surface (in const Z-plane) with a high-valence region of size A_f*(surface_area) and a valency of n
void refacet_surface( moab::EntityHandle surf, double A_f, int valence );
// returns the verts of an empty square in the center of the surface and surrounds this square w/ triangles in a watertight fashion
void generate_box_space( moab::EntityHandle surf, double A_f, std::vector<moab::EntityHandle> &box_verts );
// creates the triangles for the high-valence area
void make_hv_region( moab::EntityHandle surf, std::vector<moab::EntityHandle> box_verts, int n );
// returns a surface that is constant in Z
void get_hv_surf( MEntVector surfs, moab::EntityHandle &hv_surf );
// removes and deletes all triangles in the given surface
void tear_down_surface( moab::EntityHandle surf );
// returns the area of a polygon given the ordered verts
double polygon_area( std::vector<moab::EntityHandle> verts );

////// Functions for analyzing the hv region \\\\\\\\\\\\\

static const double PI = acos(-1.0);
static const double denom = 1.0 / ((double) RAND_MAX);
static const double denomPI = PI * denom;


static double facet_tol = 1e-4;
static double source_rad = 0;
static int vol_index = 1;
static int num_random_rays = 1000;
static int randseed = 12345;
static bool do_stat_report = false;
static bool do_trv_stats   = false;
static double location_az = 2.0 * PI;
static double direction_az = location_az;
static const char* pyfile = NULL;

inline void RNDVEC(CartVect& uvw, double &az) 
{
  
  double theta = denom * az * rand();
  double phi = denomPI * rand();
  uvw[0] = cos(theta)*sin(phi);
  uvw[1] = sin(theta)*sin(phi);
  uvw[2] = cos(phi);

}


moab::ErrorCode write_obb_mesh( moab::DagMC *dag, moab::EntityHandle vol, std::string& base_filename);

moab::ErrorCode get_volumes( moab::Interface* mb, moab::Range &volumes);


void get_time_mem(double &tot_time, double &user_time,
                  double &sys_time, double &tot_mem); 

void fire_rand_rays( moab::DagMC *dagi, moab::EntityHandle vol, int num_rand_rays, double &avg_fire_time, moab::CartVect ray_source);

int main(int argc, char **argv)
{
  /*
  //hangle arguments
  if( 3 != argc  ) 
    {
      std::cout << "Please include all necessary arguments! Exiting..." << std::endl;
      std::cout << "To generate a high valence region with area ";
      std::cout << "fraction A_f (double) and verts of valency n:" << std::endl;
      std::cout << "$ ./create_hv_mesh <A_f> <n>" << std::endl;
      return 0;
    }
  */
  double A_f = 0;
  int valence = 0;
  
  std::ofstream data_file, param_file;
  
  param_file.open("params1.dat");
  data_file.open("data1.dat");
  //get the mesh ready_using MeshKit (easier for manipulating meshes)
  mk = new MKCore();

  param_file << 0;

  int area_intervals = 10;
  int valence_intervals = 3*10;
  double max_A_f = 1.0/6.0;
  int max_n = 1e5;

  for(unsigned int i=1; i < area_intervals; i++)
    {

      A_f = (double)i * ( max_A_f / area_intervals);
      //write this value to params file everytime

      for(unsigned int j=1; j < valence_intervals; j++)
	{
	  valence = (double)j * ( max_n / valence_intervals );
	  //the first time we go through the inner loop, 
	  //write all of the valence values
	  if ( 1 == i ) param_file << "\t" << valence << "\t";


	  prep_mesh( A_f, valence);

	  ////////////// START OF MOAB STUFF \\\\\\\\\\\\\\\\\\\\\\\	\
	  
	  //now we'll try to load this mesh-file into a dagmc instance
	  moab::DagMC *dag = moab::DagMC::instance();

	  moab::ErrorCode result;
	  //try loading the file 
	  result = dag->load_file( "cube_mod.h5m" );
	  if( MB_SUCCESS != result) return 1;

	  //generate the OBB tree
	  result = dag->init_OBBTree();
	  if( MB_SUCCESS != result) return 1;

	  //get all of the volumes in the dagmc instance
	  moab::Range vols;
	  result = get_volumes( dag->moab_instance(), vols);
	  if( MB_SUCCESS != result) return 1; 
	  
	  //write the obbs to a new set of files based on depth in the tree
	  std::string dum = "test";
	  result = write_obb_mesh( dag, vols[0], dum);
	  if( MB_SUCCESS != result) return 1;  
	  
	  
	  //analyze mesh here
	  double avg_fire_time;
	  CartVect source;
	  source[0] = 0; source [1] = 0; source [2] = 0;
	  //call into the new functions for firing random rays and get the avg time
	  fire_rand_rays( dag, vols[0], 100000, avg_fire_time, source);
	  //write time to data file
	  
	  std::cout << "The average fire time for this mesh was: " << avg_fire_time << "s" << std::endl;

	  data_file << avg_fire_time << "\t";
	  dag->moab_instance()->delete_mesh();
	}
      //add a new line to the params file
      param_file << std::endl;
      param_file << A_f << "\t";
      data_file << std::endl;
    }

  param_file.close();
  data_file.close();
  
  return 0;

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

  HexMaker op1( obbtool, &mbi2 );

  moab::OrientedBoxTreeTool::TrvStats tree_stats;

  rval = obbtool->preorder_traverse( root, op1, &tree_stats );
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 

  rval = op1.write_to_files( base_filename );
  assert( MB_SUCCESS == rval );
  if( MB_SUCCESS != rval ) return rval; 


}
  


void fire_rand_rays( moab::DagMC *dagi, moab::EntityHandle vol, int num_rand_rays, double &avg_fire_time, moab::CartVect ray_source)
{

  srand(12345);

  moab::CartVect xyz, uvw;

  double ttime1, utime1, stime1, tmem1, ttime2, utime2, stime2, tmem2;
  get_time_mem(ttime1, utime1, stime1, tmem1);

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
    

    //    std::cout << "x,y,z,u,v,w,u^2 + v^2 + w^2 = " << xyz 
    //        << " " << uvw << " " << uvw%uvw << std::endl;
    //uavg += uvw[0]; vavg += uvw[1]; wavg += uvw[2];

    // added ray orientation
    dagi->ray_fire(vol, xyz.array(), uvw.array(), dum, dist, NULL, 0, 1, &trv_stats );
    
    if( dum == 0){ random_rays_missed++; }
    
  }
  //end timer 
  end_time = std::clock();

  get_time_mem(ttime2, utime2, stime2, tmem1);
  double timewith = ttime2 - ttime1;
  
  srand(randseed); // reseed to generate the same values as before
  
  // now without ray fire call, to subtract out overhead
  for (int j = 0; j < num_random_rays; j++) {
    RNDVEC(uvw, location_az);
    
    xyz = uvw * source_rad + ray_source;
    if (source_rad >= 0.0) {
      RNDVEC(uvw, direction_az);
    }
  }
  
  get_time_mem(ttime1, utime1, stime1, tmem2);
  double timewithout = ttime1 - ttime2;
  
  std::cout << " done." << std::endl;
  
  if( random_rays_missed ){
    std::cout << "Warning: " << random_rays_missed << " random rays did not hit the target volume" << std::endl;
  }
  
  if( num_random_rays > 0 ){
    std::cout << "Total time per ray fire: " << timewith/num_random_rays 
	      << " sec" << std::endl;
    std::cout << "Estimated time per call (excluding ray generation): " 
	      << (timewith - timewithout) / num_random_rays << " sec" << std::endl;
  }
 
  avg_fire_time = (end_time - start_time) / (double)(CLOCKS_PER_SEC/1000); 
  std::cout << "My timer says: " << avg_fire_time << "ms" << std::endl;
}


void get_time_mem(double &tot_time, double &user_time,
                  double &sys_time, double &tot_mem) 
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  user_time = (double)r_usage.ru_utime.tv_sec +
    ((double)r_usage.ru_utime.tv_usec/1.e6);
  sys_time = (double)r_usage.ru_stime.tv_sec +
    ((double)r_usage.ru_stime.tv_usec/1.e6);
  tot_time = user_time + sys_time;
  tot_mem = 0;

  // try going to /proc to estimate total memory
    char file_str[4096], dum_str[4096];
    int file_ptr = -1, file_len;
    file_ptr = open("/proc/self/stat", O_RDONLY);
    file_len = read(file_ptr, file_str, sizeof(file_str)-1);
    if (file_len == 0) return;
    
    close(file_ptr);
    file_str[file_len] = '\0';
      // read the preceeding fields and the ones we really want...
    int dum_int;
    unsigned int dum_uint, vm_size, rss;
    int num_fields = sscanf(file_str, 
                            "%d " // pid
                            "%s " // comm
                            "%c " // state
                            "%d %d %d %d %d " // ppid, pgrp, session, tty, tpgid
                            "%u %u %u %u %u " // flags, minflt, cminflt, majflt, cmajflt
                            "%d %d %d %d %d %d " // utime, stime, cutime, cstime, counter, priority
                            "%u %u " // timeout, itrealvalue
                            "%d " // starttime
                            "%u %u", // vsize, rss
                            &dum_int, 
                            dum_str, 
                            dum_str, 
                            &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
                            &dum_uint, &dum_uint, &dum_uint, &dum_uint, &dum_uint,
                            &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
                            &dum_uint, &dum_uint, 
                            &dum_int,
                            &vm_size, &rss);
    if (num_fields == 24)
      tot_mem = ((double)vm_size);

}
