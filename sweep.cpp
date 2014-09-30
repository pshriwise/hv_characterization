// The purpose of this file is to creat a bash script which will sweep
// the available parameter spaces of A_f and n

#include <iostream>
#include <fstream>
#include "hv_mesh_gen.hpp"
#include "ProgOptions.hpp"

void gnuplot_script();

int main(int argc, char** argv)
{

  int area_intervals = 10;
  int valence_intervals = 4*area_intervals;
  double max_A_f;
  int max_n;


  ProgOptions po("sweep: a program for sweeping two parameter spaces for a manually modified high-valence mesh");

  po.addRequiredArg<double>( "af", "Fraction of the cube surface area that should be made high-valence.", &max_A_f);

  po.addRequiredArg<int>( "v", "Valence of the high-valence region", &max_n);

  po.addOpt<int>( "af_ints", "Number of intervals for the Area fraction param", &area_intervals );

  po.addOpt<int>( "v_ints", "Number of intervals for the valence param", &valence_intervals);

  po.parseCommandLine( argc, argv );

  double A_f = 0;
  int valence = 0;

  gnuplot_script();
  std::ofstream data_file, param_file;
  
  param_file.open("params.dat");
  data_file.open("data.dat");
  //start the data file with an empty line to make the future paste easier
  data_file << std::endl;
  //get the mesh ready_using MeshKit (easier for manipulating meshes)

  param_file << 0;


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


void gnuplot_script()
{

  std::ofstream gp_file; 
  gp_file.open("timing_plot.p");
    
  gp_file << "set xlabel 'Area Fraction (A_f)' " << std::endl;
  gp_file << "set ylabel 'Valence (n)'" << std::endl;
  gp_file << "set zlabel 'Avg Ray Fire Time (ms)' " << std::endl;

  gp_file << "splot 'all.dat' matrix nonuniform w lines" << std::endl;

  gp_file << "pause -1" << std::endl;

}





