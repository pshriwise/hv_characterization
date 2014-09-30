
#include <iostream>
#include "hv_mesh_gen.hpp"

int main(int argc, char** argv)

{

  if( 3 != argc)
    {
      std::cout << "Please include all necessary arguments! Exiting..." << std::endl;
      std::cout << "To generate a high valence region with area ";
      std::cout << "fraction A_f (double) and verts of valency n:" << std::endl;
      std::cout << "$ ./create_hv_mesh <A_f> <n>" << std::endl;
    }


  double A_f = atof(argv[1]);

  int valence = atoi(argv[2]);

  //create the hv cube and write it to file
  prep_mesh( A_f, valence );

  return 0;

}
