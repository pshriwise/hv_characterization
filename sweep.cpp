// The purpose of this file is to creat a bash script which will sweep
// the available parameter spaces of A_f and n

#include <iostream>
#include <fstream>
#include "hv_mesh_gen.hpp"

void gnuplot_script( int A_f, int valence);

int main(int argc, char** argv)
{

  double A_f = 0;
  int valence = 0;
  
  std::ofstream data_file, param_file;
  
  param_file.open("params1.dat");
  data_file.open("data1.dat");
  //get the mesh ready_using MeshKit (easier for manipulating meshes)

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

	  double fire_time = 0;
	  std::cout << "Area fraction: " << A_f << std::endl << "Valence: " << valence << std::endl;

	  test_hv_cube_mesh( A_f, valence, fire_time);
	  data_file << fire_time << "\t";

	}

      param_file << std::endl;
      param_file << A_f << "\t";
      data_file << std::endl;
    }

  param_file.close();
  data_file.close();
  
  return 0;


}


void gnuplot_script( int area_intervals, int valence_intervals)
{

  std::ofstream gp_file; 
  gp_file.open("timing_plot.p");
    
  gp_file << "set dgrid3d " << area_intervals << ", " << valence_intervals << std::endl;
  gp_file << "set style data lines" << std::endl;

  //gp_file << "set hidden3d" << std::endl;

  gp_file << "set xlabel 'Area Fraction (A_f)' " << std::endl;
  gp_file << "set ylabel 'Valence (n)'" << std::endl;
  gp_file << "set zlabel 'Avg Ray Fire Time (micro-s)' " << std::endl;

  gp_file << "splot 'data.dat' u 1:2:($3 * 1000000) notitle" << std::endl;

  gp_file << "pause -1" << std::endl;

}





