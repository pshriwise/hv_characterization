
#include <iostream>
#include <assert.h>

//Sigma Includes
#include "DagMC.hpp"
#include "moab/OrientedBoxTreeTool.hpp"

//local includes
#include "hv_mesh_gen.hpp"
#include "gen_mb_funcs.hpp"
#include "ray_fire.hpp"

moab::Interface *mbi = new moab::Core();

moab::ErrorCode build_obbs(moab::DagMC* dag, moab::Range surfs, moab::Range vols, moab::OrientedBoxTreeTool::Settings settings) {
  ErrorCode rval = MB_SUCCESS;

  for (Range::iterator i = surfs.begin(); i != surfs.end(); ++i) {
    EntityHandle root;
    Range tris;
    rval = mbi->get_entities_by_dimension( *i, 2, tris );
    if (MB_SUCCESS != rval)
      return rval;
    if (tris.empty())
      std::cerr << "WARNING: Surface " << *i << " has no facets." << std::endl;
    rval = dag->obb_tree()->build( tris, root, &settings );
    if (MB_SUCCESS != rval)
      return rval;
    rval = mbi->add_entities( root, &*i, 1 );
    if (MB_SUCCESS != rval)
      return rval;
    rval = mbi->tag_set_data( dag->obb_tag(), &*i, 1, &root );
    if (MB_SUCCESS != rval)
      return rval;
  }

  for (Range::iterator i = vols.begin(); i != vols.end(); ++i) {
    // get all surfaces in volume
    moab::Range tmp_surfs;
    rval = mbi->get_child_meshsets( *i, tmp_surfs );
    if (MB_SUCCESS != rval)
      return rval;

    // get OBB trees for each surface
    moab::EntityHandle root;
    moab::Range trees;
    for (Range::iterator j = tmp_surfs.begin();  j != tmp_surfs.end(); ++j) {
      // skip any surfaces that are non-manifold in the volume
      // because point containment code will get confused by them
      int sense = 0;
      rval = dag->surface_sense( *i, *j, sense );
      if (MB_SUCCESS != rval) {
        std::cerr << "Surface/Volume sense data missing." << std::endl;
        return rval;
      }
      if (!sense)
        continue;
      rval = mbi->tag_get_data( dag->obb_tag(), &*j, 1, &root );
      if (MB_SUCCESS != rval || !root) return MB_FAILURE;
      trees.insert( root );
    }
    
    // build OBB tree for volume
    rval = dag->obb_tree()->join_trees( trees, root );
    if (MB_SUCCESS != rval) return rval;

    rval = mbi->tag_set_data( dag->obb_tag(), &*i, 1, &root );
    if (MB_SUCCESS != rval) return rval;
  }

  return MB_SUCCESS;
}

moab::ErrorCode test_hv_cube_mesh( double A_f, double valence, double &ray_fire_time, double worst_split_ratio )
{

  std::cout << "TESTING MESH - Valence: " << valence << " Area Fraction: " << A_f << std::endl;
	  prep_mesh( A_f, valence);
	  
	  ////////////// START OF MOAB STUFF \\\\\\\\\\\\\\\\\\\\\\\	\

	  moab::OrientedBoxTreeTool::Settings settings; 
	  settings.max_leaf_entities = 8; 
	  settings.max_depth = 0;
	  settings.worst_cost = worst_split_ratio; 
	  settings.best_cost = 0.4;
	  settings.set_options = MESHSET_SET;

	  //now we'll try to load this mesh-file into a dagmc instance
	  moab::DagMC *dag = new moab::DagMC( mbi );
	  
	  moab::ErrorCode result;
	  
	  //load the mesh data from the moab instance into dag
	  result = dag->load_existing_contents(); 
	  if( MB_SUCCESS != result) return MB_FAILURE;

	  result = dag->setup_impl_compl();
	  MB_CHK_SET_ERR(result, "Could not setup implicit compliment.");

	  moab::Range surfs,vols;
	  result = dag->setup_geometry(surfs,vols);
	  MB_CHK_SET_ERR(result,"Failed to setup the geometry.");
	  // build obbs
	  result = build_obbs(dag,surfs,vols,settings);
	  MB_CHK_SET_ERR(result, "Failed to setup the OBBs");
	  
	  // setup indices
	  result = dag->setup_indices();
	  MB_CHK_SET_ERR(result, "Failed to setup problem indices");
	  
	  //analyze mesh here
	  double avg_fire_time;
	  CartVect source;
	  source[0] = 0; source [1] = 0; source [2] = 0;
	  //call into the new functions for firing random rays and get the avg time
	  fire_rand_rays( dag, vols[0], 100000, avg_fire_time, source);
	  //write time to data file

	  mbi->write_mesh("hvcube.h5m");
	  std::cout << "The average fire time for this mesh was: " << avg_fire_time << "s" << std::endl;
	  ray_fire_time = avg_fire_time;
	  
	  mbi->delete_mesh();
	  return MB_SUCCESS;

}


void create_cube()
{
  // Define a 2x2x2 cube centered at orgin
  // with concavity in +Z face.
  const double coords[] = {
    50, -50, -50, 
    50,  50, -50,
   -50,  50, -50,
   -50, -50, -50,
    50, -50,  50, 
    50,  50,  50,
   -50,  50,  50,
   -50, -50,  50};
  const int connectivity[] = {
    0, 3, 1,  3, 2, 1, // -Z
    0, 1, 4,  5, 4, 1, // +X
    1, 2, 6,  6, 5, 1, // +Y
    6, 2, 3,  7, 6, 3, // -X
    0, 4, 3,  7, 3, 4, // -Y
    4, 5, 7,  5, 6, 7, // +Z
  };
  unsigned tris_per_surf[] = { 2, 2, 2, 2, 2, 2 };

   // Create the geometry
  const unsigned num_verts = sizeof(coords) / (3*sizeof(double));
  const unsigned num_tris = sizeof(connectivity) / (3*sizeof(int));
  const unsigned num_surfs = sizeof(tris_per_surf) / sizeof(unsigned);
  EntityHandle verts[num_verts], tris[num_tris], surfs[num_surfs];
  moab::ErrorCode rval;
  for (unsigned i = 0; i < num_surfs; ++i) {
    rval = mbi->create_meshset( MESHSET_SET, surfs[i] );
    MB_CHK_ERR_CONT(rval);
  }
  for (unsigned i = 0; i < num_verts; ++i) {
    rval = mbi->create_vertex( coords + 3*i, verts[i] ); 
    MB_CHK_ERR_CONT(rval);
  }
  EntityHandle *surf_iter = surfs;
  unsigned *tps = tris_per_surf;
  for (unsigned i = 0, j = 0; i < num_tris; ++i, ++j) {
    const EntityHandle conn[] = { verts[connectivity[3*i  ]], 
                                    verts[connectivity[3*i+1]], 
                                    verts[connectivity[3*i+2]] };
    rval = mbi->create_element( MBTRI, conn, 3, tris[i] );
    MB_CHK_ERR_CONT(rval);
    rval = mbi->add_entities(*surf_iter, &tris[i], 1);
    MB_CHK_ERR_CONT(rval);
    rval = mbi->add_entities(*surf_iter, &conn[0], 3);
    MB_CHK_ERR_CONT(rval);
    if (j == *tps) {
      surf_iter++;
      tps++;
      j = 0;
    }
  }
    
  Tag dim_tag, id_tag, sense_tag, category_tag;
  rval = mbi->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 
                              1, MB_TYPE_INTEGER, 
                              dim_tag,
                              MB_TAG_SPARSE|MB_TAG_CREAT );
  MB_CHK_ERR_CONT(rval);
  rval = mbi->tag_get_handle( GLOBAL_ID_TAG_NAME, 
                              1, MB_TYPE_INTEGER, 
                              id_tag,
                              MB_TAG_DENSE|MB_TAG_CREAT );
  MB_CHK_ERR_CONT(rval);
  rval = mbi->tag_get_handle( CATEGORY_TAG_NAME, 
                              CATEGORY_TAG_SIZE, MB_TYPE_OPAQUE, 
                              category_tag,
                              MB_TAG_SPARSE|MB_TAG_CREAT );
  MB_CHK_ERR_CONT(rval);
  rval = mbi->tag_get_handle( "GEOM_SENSE_2", 
                              2, MB_TYPE_HANDLE, 
                              sense_tag,
                              MB_TAG_SPARSE|MB_TAG_CREAT );
  MB_CHK_ERR_CONT(rval);

  std::vector<int> dims( num_surfs, 2 );
  rval = mbi->tag_set_data( dim_tag, surfs, num_surfs, &dims[0] );
  MB_CHK_ERR_CONT(rval);
  char surf[CATEGORY_TAG_SIZE] = "Surface";
  for( unsigned int i = 0; i < num_surfs; i++) {
    rval = mbi->tag_set_data( category_tag, &surfs[i], num_surfs, &surf[0]);
    MB_CHK_ERR_CONT(rval);
  }
  std::vector<int> ids( num_surfs );
  for (size_t i = 0; i < ids.size(); ++i) ids[i] = i+1;
  rval = mbi->tag_set_data( id_tag, surfs, num_surfs, &ids[0] );
  MB_CHK_ERR_CONT(rval);

  EntityHandle volume;
  rval = mbi->create_meshset( MESHSET_SET, volume );
  MB_CHK_ERR_CONT(rval);
  for (unsigned i = 0; i < num_surfs; ++i) {
    rval = mbi->add_parent_child( volume, surfs[i] );
    MB_CHK_ERR_CONT(rval);
  }
  
  std::vector<EntityHandle> senses( 2*num_surfs, 0 );
  for (size_t i = 0; i < senses.size(); i += 2)
    senses[i] = volume;
  rval = mbi->tag_set_data( sense_tag, surfs, num_surfs, &senses[0] );
  MB_CHK_ERR_CONT(rval);
  
  const int three = 3;
  const int one = 1;
  char volume_cat[CATEGORY_TAG_SIZE] = "Volume";
  rval = mbi->tag_set_data( dim_tag, &volume, 1, &three );
  MB_CHK_ERR_CONT(rval);
  rval = mbi->tag_set_data( id_tag, &volume, 1, &one );
  MB_CHK_ERR_CONT(rval);
  rval = mbi->tag_set_data( category_tag, &volume, 1, &volume_cat[0]);
  MB_CHK_ERR_CONT(rval);
}


void prep_mesh( double A_f, int valence )
{

  //create the initial geometry and mesh it
  create_cube();

  //Get all of the surface ModelEnts
  
  //Now to find one of the surfaces that is constant in z (for convenience)
  moab::EntityHandle hv_surf;
  get_hv_surf( hv_surf );

  //refacet the surface using the desired area fraction for the hv region
  refacet_surface( hv_surf, A_f, valence );

}

void refacet_surface( moab::EntityHandle surf, double A_f, int valence )
{

  //remove all 2-D entities on the surfce
  tear_down_surface( surf );

 
  //now its time to create an empty middle box using the remaining surface verts
  std::vector<moab::EntityHandle> box;
  generate_box_space( surf, A_f, box, 1 );

  make_hv_region( surf, box, valence );

}

// returns the verts of an empty box, centered on the origin which is surrounded by triangles
void generate_box_space( moab::EntityHandle surf, double A_f, std::vector<moab::EntityHandle> &box_verts, int axis )
{
  int x_idx, y_idx;
  //set the indices of the cartesian vectors based on the constant axis for this surface 
  switch(axis){
  case 0:
    x_idx = 1; 
    y_idx = 2; 
    break;
  case 1:
    x_idx = 0;
    y_idx = 2; 
    break; 
  case 2:
    x_idx = 0;
    y_idx = 1;
    break;
  default:
    std::cout << "Value of axis muxt be 0, 1, or 2" << std::endl; 
    assert(false); 
  }

  std::vector<moab::EntityHandle> corners;
  moab::ErrorCode rval = mbi->get_entities_by_type( surf, MBVERTEX, corners );
  MB_CHK_ERR_CONT(rval);

  double surface_area = polygon_area( corners );
  //going to start taking advantage of knowing the geometry here...
  double surface_side = sqrt(surface_area);
  double cube_area = 6*surface_area;

  //now create the vertices for the new regions
  
  //based on surface area, get the length of one of the center-square sides
  if( A_f == 1 )
    {
      std::vector<moab::EntityHandle> ordered_corners;
      order_corner_verts(corners, ordered_corners);
      box_verts = ordered_corners;
      return;
    }

  double hv_area = A_f*cube_area/6; //divide by 6 because this surface represents all surfaces of the cube
  if ( hv_area >= surface_area ) std::cout << "ERROR: Area fraction must be less than 1/6 for now. " << surf << std::endl;
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
      moab::CartVect coords;
      mbi->get_coords( &(*i), 1, coords.array() );
      
      //need three cartesian vectors telling us where to put the verts, two for the dam points, one for the inner point.
      moab::CartVect to_ewdam, to_nsdam, to_box;
      to_ewdam[axis] = 0; to_nsdam[axis] = 0; to_box[axis] = 0; //only moving on the x-y plane here
      
      //always want to create this vert
      //using the fact that we're centered on the origin here...
      //x-points
      if( coords[x_idx] < 0 ){ 
	to_box[x_idx] = box_bump_dist; }
      else  { 
	to_box[x_idx] = -box_bump_dist; }
      //y-points
      if( coords[y_idx] < 0 ){ 
	to_box[y_idx]= box_bump_dist; }
      else  { 
	to_box[y_idx] = -box_bump_dist; }

      //create the inner point
      moab::EntityHandle box_vert;
      mbi->create_vertex( (coords+to_box).array(), box_vert);
      box.push_back(box_vert);
      

      //figure out which dams we're adjacent to
      std::vector<moab::EntityHandle> *nsdam = &nd, *ewdam = &ed;
      moab::EntityHandle *nsdam_vert = &ndv, *ewdam_vert = &edv;

      //set the north-south dam based on our y location
      if( coords[y_idx] < 0) { nsdam_vert = &sdv; nsdam = &sd; }
      else { nsdam_vert = &ndv; nsdam = &nd; }

      //set the east-west dam based on our x location
      if( coords[x_idx] < 0) { ewdam_vert = &wdv; ewdam = &wd; }
      else { ewdam_vert = &edv; ewdam = &ed; }

      //////NORTH-SOUTH DAM\\\\\\\

      //if the dams already have info from another corner, no need to create it's point,
      //just add this corner to the dam triangle
      if( 0 != nsdam->size() && NULL != nsdam_vert ) { nsdam->push_back(*i); }
      // if the dam has no info, then we'll create it
      else {
	
	//set the dam coordinates
	if ( coords[x_idx] < 0 ) to_nsdam[x_idx] = dam_bump_long;
	else to_nsdam[x_idx] = -dam_bump_long;
	if ( coords[y_idx] < 0 ) to_nsdam[y_idx] = dam_bump_short;
	else to_nsdam[y_idx] = -dam_bump_short;

	//create the dam vertex 
	mbi->create_vertex( (coords+to_nsdam).array(), *nsdam_vert);
	
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
	if ( coords[x_idx] < 0 ) to_ewdam[x_idx] = dam_bump_short;
	else to_ewdam[x_idx] = -dam_bump_short;
	if ( coords[y_idx] < 0 ) to_ewdam[y_idx] = dam_bump_long;
	else to_ewdam[y_idx] = -dam_bump_long;

	//create the dam vertex 
	mbi->create_vertex( (coords+to_ewdam).array(), *ewdam_vert);
	
	//now add this corner and the dam vert to the dam verts list
	ewdam->push_back(*ewdam_vert); ewdam->push_back(*i);

      } 
	
      //dam info should now all be set

      //there are always two triangles to create that connect the dam
      //to the corner, we'll do that now
      moab::CartVect corner_coords;
      mbi->get_coords(&(*i), 1, corner_coords.array());
      moab::EntityHandle tri1_verts[3] = {*i, *nsdam_vert, box_vert};
      moab::EntityHandle tri2_verts[3] = {*i, *ewdam_vert, box_vert};
      if( coords[x_idx] < 0 ) {
	if ( coords[y_idx] > 0 ) {
	  tri1_verts[1] = *nsdam_vert;
	  tri1_verts[2] = box_vert;
	  tri2_verts[1] = box_vert;
	  tri2_verts[2] = *ewdam_vert;
	}
	else {
	  tri1_verts[1] = box_vert;
	  tri1_verts[2] = *nsdam_vert;
	}
      }
      else {
	if ( coords[y_idx] > 0 ) {
	  tri1_verts[1] = box_vert;
	  tri1_verts[2] = *nsdam_vert;
	}
	else {
	  tri2_verts[1] = box_vert;
	  tri2_verts[2] = *ewdam_vert;
	}
      }
      std::vector<moab::EntityHandle> tris(2);
      mbi->create_element( MBTRI, &tri1_verts[0], 3, tris[0] );
      mbi->create_element( MBTRI, &tri2_verts[0], 3, tris[1] );
      
      //now we'll add these to the set
      mbi->add_entities( surf, &tri1_verts[0], 3 );
      mbi->add_entities( surf, &tri2_verts[0], 3 );      
      mbi->add_entities( surf, &(tris[0]), tris.size() );
    }

  //now we should have all of the info needed to create our dam triangles
  std::vector<moab::EntityHandle> dam_tris(4);
  assert( nd.size()==3 );
  moab::EntityHandle temp = nd[1];
  nd[1] = nd[2];
  nd[2] = temp;
  temp = wd[1];
  wd[1] = wd[2];
  wd[2] = temp;
  mbi->create_element( MBTRI, &(nd[0]), nd.size(), dam_tris[0]);
  mbi->create_element( MBTRI, &(sd[0]), sd.size(), dam_tris[1]);
  mbi->create_element( MBTRI, &(ed[0]), ed.size(), dam_tris[2]);
  mbi->create_element( MBTRI, &(wd[0]), wd.size(), dam_tris[3]);

  //add these to the surface
  mbi->add_entities( surf, &(dam_tris[0]), dam_tris.size() );

  
  //re-order the M verts
  std::vector<moab::EntityHandle> ordered_box(4);
  order_corner_verts(box, ordered_box);

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
      mbi->create_element( MBTRI, &(tri_verts[0]), 3, last_few_tris[i]);
    }
      
  //now add these to the surface
  mbi->add_entities( surf, &last_few_tris[0], 4);

  box_verts = ordered_box;
  
} // end generate_box_space

//creates a watertight high-valence area with n-valencies per
void make_hv_region( moab::EntityHandle surf, std::vector<moab::EntityHandle> box_verts, int n ) 
{


  //we're asked for 0 triangles in the hv area, return at least one
  if ( 0 == n )
    {
      make_hv_region( surf, box_verts, 1); 
      return; 
    }

  //start by making a paramaterization of the box diagonal
  //because we know the order of these verts, we can get the diagonal coordinates directly
  moab::CartVect sw_coords, ne_coords;

  mbi->get_coords( &box_verts[0], 1, sw_coords.array() );
  mbi->get_coords( &box_verts[2], 1, ne_coords.array() );

  //now create the param vector
  moab::CartVect sw_to_ne = ne_coords-sw_coords;

  // now, based on the valence asked for, get the new vertex coords long this line
  std::vector<moab::EntityHandle> diag_verts;
 

  for(unsigned int i = 1; i <= n-1; i++)
    {
      double u = double(i)/double(n);
      moab::CartVect new_coords = sw_coords + u*sw_to_ne;
      moab::EntityHandle new_vert;
      mbi->create_vertex( new_coords.array(), new_vert);
      diag_verts.push_back(new_vert);
    }

  //add all of these new verts to the surface
  mbi->add_entities( surf, &(diag_verts[0]), diag_verts.size() );


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
      moab::EntityHandle tri1_verts[3] = { nw_vert, diag_verts[i+1], diag_verts[i] };
      moab::EntityHandle tri2_verts[3] = { se_vert, diag_verts[i], diag_verts[i+1] };
      mbi->create_element( MBTRI, &tri1_verts[0], 3, tri1);
      mbi->create_element( MBTRI, &tri2_verts[0], 3, tri2);

      //now add these to the surface set
      mbi->add_entities( surf, &tri1, 1);
      mbi->add_entities( surf, &tri2, 1);
    }
  
}

void get_hv_surf(moab::EntityHandle &hv_surf)
{
  char category[CATEGORY_TAG_SIZE] = "Surface";
   //get the name tag from the moab instance
  moab::Tag category_tag; 
  moab::ErrorCode result = mbi->tag_get_handle(CATEGORY_TAG_NAME, category_tag); 
  MB_CHK_ERR_CONT(result); 

  //create void pointer for tag data match
  const void *dum = &(category[0]);
  moab::Range surfs;
  result = mbi->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &category_tag, &dum, 1, surfs);
  MB_CHK_ERR_CONT(result);

  std::cout << "There are " << surfs.size() << " surfaces in this model." << std::endl;
  assert( 6 == surfs.size() ); 

  for( unsigned int i = 0 ; i < surfs.size() ; i++)
    {
      //get the meshset handle for this surface
      moab::EntityHandle sh = surfs[i];
      
      //get the triangles in this meshset
      std::vector<moab::EntityHandle> tris;
      mbi->get_entities_by_type( sh, MBTRI, tris);
      
      //setup a check vector
      moab::CartVect check;
      check [0] = 0; check [1] = 1; check[2] = 0;
      //get the normal for the first triangle on each surface
      
      //get the verts for each triangle and their coordinates
      std::vector<moab::EntityHandle> verts;
      mbi->get_connectivity( &(tris[0]), 1, verts);
      
      moab::CartVect coords[3];
      
      mbi->get_coords( &(verts[0]), 3, coords[0].array() );
      
      moab::CartVect tri_norm;
      tri_norm = (coords[1]-coords[0])*(coords[2]-coords[0]);
      tri_norm.normalize();
      
      if( tri_norm == check){
	hv_surf = sh;
	break;
      }

    } //end loop  
} //end get_hv_surf

void tear_down_surface( moab::EntityHandle surf)
{

 //get the triangles for this surface and delete them
  std::vector<moab::EntityHandle> tris;
  moab::ErrorCode rval = mbi->get_entities_by_type( surf, MBTRI, tris);
  MB_CHK_ERR_CONT(rval);
  
  //remove these triangle from the meshset and destroy them
  std::cout << "Tearing down the surface..." << std::endl; 

  rval = mbi->remove_entities( surf, &(tris[0]), tris.size());
  MB_CHK_ERR_CONT(rval);

  rval = mbi->delete_entities( &(tris[0]), tris.size());
  MB_CHK_ERR_CONT(rval);

  std::vector<moab::EntityHandle> test_tris;
  rval = mbi->get_entities_by_type( 0, MBTRI, test_tris, true);
  MB_CHK_ERR_CONT(rval);
  std::cout << "There are now " << test_tris.size() << " triangles in the model" << std::endl;

}

double polygon_area( std::vector<moab::EntityHandle> verts)
{
  unsigned int num_vertices = verts.size();

  //std::cout << "There are " << num_vertices << " vertices in this surface." << std::endl;

  //get the coordinates of the vertices
  double vert_coords[verts.size()*3];
  mbi->get_coords( &(verts[0]), verts.size(), &(vert_coords[0]) );

  const moab::CartVect* coords = reinterpret_cast<const moab::CartVect*>(vert_coords);  

  //Calculate the surface's area
  moab::CartVect mid(0,0,0);
  for (int i = 0; i < num_vertices; ++i) mid += coords[i];
  mid /= num_vertices;

  double sum = 0.0;
  for (int i = 0; i < num_vertices; i++)
    {
      int j = (i+1)%(num_vertices);
      moab::CartVect pnt1, pnt2;
      pnt1 = coords[i];
      pnt2 = coords[j];
      sum += ((mid - pnt1) * (mid - pnt2)).length();
    }
  double  poly_area = sum;

  //std::cout << "The total area of this polygon is: " << poly_area << std::endl;
  
  return poly_area;
}

void save_mesh(std::string filename)
{

  mbi->write_mesh( filename.c_str() );

}

void order_corner_verts(std::vector<moab::EntityHandle> corners, std::vector<moab::EntityHandle>& ordered_corners)
{
  ordered_corners.clear(); ordered_corners.resize(4);
  for(unsigned int i = 0; i < corners.size() ; i++)
    {
      moab::CartVect v_coords;
      mbi->get_coords( &(corners[i]), 1, v_coords.array() );
      double x = v_coords[0]; double y = v_coords[2];
      // again, our surface is cenetered on the origin
      if ( x < 0 && y < 0) ordered_corners[0] = corners[i];
      else if ( x < 0 && y > 0) ordered_corners[1] = corners[i];
      else if ( x > 0 && y > 0) ordered_corners[2] = corners[i];
      else if ( x > 0 && y < 0) ordered_corners[3] = corners[i];
    }
}
