

#include "moab/Types.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "MBTagConventions.hpp"
#include "moab/ProgOptions.hpp"
#include "ray_writer.hpp"
#include "DagMC.hpp"

int main (int argc, char** argv)
{

  ProgOptions po("Traverse Ray: A program for writing out all OBBs encoutnered by a single ray traversal");

  int vol_id;
  double ray_origin[3], ray_dir[3];
  std::string fn;
  po.addRequiredArg<double>("x", "Ray origin's x position", &(ray_origin[0]));
  po.addRequiredArg<double>("y", "Ray origin's y position", &(ray_origin[1]));
  po.addRequiredArg<double>("z", "Ray origin's z position", &(ray_origin[2]));
  po.addRequiredArg<double>("u", "Ray direction u value", &(ray_dir[0]));
  po.addRequiredArg<double>("v", "Ray direction v value", &(ray_dir[1]));
  po.addRequiredArg<double>("w", "Ray direction w value", &(ray_dir[2]));
  po.addRequiredArg<int>("vol_id", "GLOBAL ID of the volume to fire the ray in", &vol_id);
  po.addRequiredArg<std::string>("filename", "DAGMC model", &fn);

  po.parseCommandLine( argc, argv );

  //start up MOAB and DAGMC
  moab::DagMC *DAG = moab::DagMC::instance();
  moab::Interface *MBI = DAG->moab_instance();
  //load the file
  moab::ErrorCode result;
  result = DAG->load_file( fn.c_str() );
  MB_CHK_SET_ERR(result, "DAGMC could not load this file");
  //create the OBB tree (should really be able to do this for only one volume...)
  result = DAG->init_OBBTree();
  MB_CHK_SET_ERR(result, "Could not create the OBB tree");
  //get the tags needed to retrieve the volume we want
  moab::Tag global_id_tag;
  result = MBI->tag_get_handle(GLOBAL_ID_TAG_NAME, global_id_tag);
  MB_CHK_SET_ERR(result, "Could not get the GLOBAL_ID tag");
  moab::Tag category_tag;
  result = MBI->tag_get_handle(CATEGORY_TAG_NAME, 32, moab::MB_TYPE_OPAQUE, category_tag);
  //MB_CHK_SET_ERR(result, "Could not get the CATEGORY tag");
  //now retrieve the volume with that id
  char cat_val[CATEGORY_TAG_SIZE] = "Volume";
  moab::Tag tags[2] = {global_id_tag,category_tag};
  void* cat_ptr = &cat_val;
  void* id_ptr = &vol_id;
  void* vals[2] = {id_ptr,cat_ptr};
  moab::Range sets;
  result = MBI->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &tags[0], &vals[0], 2, sets, moab::Interface::INTERSECT, true);
  MB_CHK_SET_ERR(result, "Could not retrieve the user-specified volume");
  //make sure we get only one volume
  if ( 1 != sets.size() ) {
    MB_CHK_SET_ERR(moab::MB_FAILURE, "Incorrect number of volumes found. Size of sets: " << sets.size());
  }
  moab::EntityHandle vol = sets[0];
  //retrieve that volume's OBB root node
  moab::EntityHandle root;
  result = DAG->get_root( vol, root );
  MB_CHK_SET_ERR( result, "Could not retrieve the root of this volume's OBB tree");
  //create the ray_writer op
  moab::OrientedBoxTreeTool *treeTool = DAG->obb_tree();
  double ray_len = 1.e8, neg_ray_len = -1.e-6;

  RayTraversalWriter rtw( treeTool, &(ray_origin[0]), &(ray_dir[0]), &ray_len, &neg_ray_len, DAG->numerical_precision(), 0, &root);

  //do the traversal, adding hexes to an EntitySet as we go
  result = treeTool->preorder_traverse( root, rtw );
  MB_CHK_SET_ERR(result, "Something went wrong while performing traversal");

  //write the desired output file
  result = rtw.write_output_file();
  MB_CHK_SET_ERR(result, "Could not write the desired output file");
  
  return 0;

}
