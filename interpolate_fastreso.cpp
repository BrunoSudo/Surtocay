#include <iostream>
#include <vector>
#include <fstream>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <cmath>
#include "interpolate_fastreso.h"

#define n 402 //Number of points in the file


//Function for reading input from FastReso and store data in arrays (vectors)
Interpolated_data TakeInput_FastReso(int id) {

  // Create filename based on id
  std::string filename = "PDGid_" + std::to_string(id) + "_total_T0.1450_Fj.out";

  //Read file
  std::ifstream file(filename);  //Variable storing the file (must not have the header)

  //Error if the file can't be read
  if (!file) {
      std::cout << "WARNING: Cannot open file" << std::endl;
  }

  // Ignore the header
  std::string line;
  for(int i = 0; i < 3; i++) {
      std::getline(file, line);
  }

  //Set vectors and variable for storing the data
  double m;
  double *pbar = new double[n];
  double Feq1[n], Feq2[n];
  double Fshear1[n], Fshear2[n], Fshear3[n];
  double Fbulk1[n], Fbulk2[n];
  double Ftemp1[n], Ftemp2[n];
  double Fvel1[n], Fvel2[n], Fvel3[n];

  //Assume the value of f(0) is 0 for all functions
  pbar[0] = 0;
  Feq1[0] = 0;
  Feq2[0] = 0;
  Fshear1[0] = 0;
  Fshear2[0] = 0;
  Fshear3[0] = 0;
  Fbulk1[0] = 0;  
  Fbulk2[0] = 0;
  Ftemp1[0] = 0;
  Ftemp2[0] = 0;
  Fvel1[0] = 0;
  Fvel2[0] = 0;
  Fvel3[0] = 0;
  

  //Read each line and store each column in its respective array inside the structure until the end of file
  int N = 1;
  while (true) {
    if (file.eof()) break;

    file >> pbar[N] >> m >> Feq1[N] >> Feq2[N] >> Fshear1[N] >> Fshear2[N] >> Fshear3[N]
          >> Fbulk1[N] >> Fbulk2[N] >> Ftemp1[N] >> Ftemp2[N] >> Fvel1[N] >> Fvel2[N] >> Fvel3[N];

    N++;
  }
  file.close();

  //structure that will be retuned when the function is called
  Interpolated_data i_data;
  
  i_data.m = m;

  //interpolating all arrays for values of pbar using gsl
  i_data.Feq1 = gsl_spline_alloc(gsl_interp_steffen, N);
  i_data.Feq2 = gsl_spline_alloc(gsl_interp_steffen, N);
  i_data.Fshear1 = gsl_spline_alloc(gsl_interp_steffen, N);
  i_data.Fshear2 = gsl_spline_alloc(gsl_interp_steffen, N);
  i_data.Fshear3 = gsl_spline_alloc(gsl_interp_steffen, N);
  i_data.Fbulk1 = gsl_spline_alloc(gsl_interp_steffen, N);
  i_data.Fbulk2 = gsl_spline_alloc(gsl_interp_steffen, N);
  i_data.Ftemp1 = gsl_spline_alloc(gsl_interp_steffen, N);
  i_data.Ftemp2 = gsl_spline_alloc(gsl_interp_steffen, N);
  i_data.Fvel1 = gsl_spline_alloc(gsl_interp_steffen, N);
  i_data.Fvel2 = gsl_spline_alloc(gsl_interp_steffen, N);
  i_data.Fvel3 = gsl_spline_alloc(gsl_interp_steffen, N);
  
  gsl_spline_init(i_data.Feq1, pbar, Feq1, N);
  gsl_spline_init(i_data.Feq2, pbar, Feq2, N);
  gsl_spline_init(i_data.Fshear1, pbar, Fshear1, N);
  gsl_spline_init(i_data.Fshear2, pbar, Fshear2, N);
  gsl_spline_init(i_data.Fshear3, pbar, Fshear3, N);
  gsl_spline_init(i_data.Fbulk1, pbar, Fbulk1, N);
  gsl_spline_init(i_data.Fbulk2, pbar, Fbulk2, N);
  gsl_spline_init(i_data.Ftemp1, pbar, Ftemp1, N);
  gsl_spline_init(i_data.Ftemp2, pbar, Ftemp2, N);
  gsl_spline_init(i_data.Fvel1, pbar, Fvel1, N);
  gsl_spline_init(i_data.Fvel2, pbar, Fvel2, N);
  gsl_spline_init(i_data.Fvel3, pbar, Fvel3, N);

  return i_data;
}


//Call to evaluate the interpolation of a given function in a given pbar
Evaluated_data Evaluate_interpolation(double pbar, Interpolated_data i_data) {
  gsl_interp_accel* acc = gsl_interp_accel_alloc();
  Evaluated_data eval;

  eval.m = i_data.m;
  eval.Feq1 = gsl_spline_eval(i_data.Feq1, pbar, acc);
  eval.Feq2 = gsl_spline_eval(i_data.Feq2, pbar, acc);
  eval.Fshear1 = gsl_spline_eval(i_data.Fshear1, pbar, acc);
  eval.Fshear2 = gsl_spline_eval(i_data.Fshear2, pbar, acc);
  eval.Fshear3 = gsl_spline_eval(i_data.Fshear3, pbar, acc);
  eval.Fbulk1 = gsl_spline_eval(i_data.Fbulk1, pbar, acc);
  eval.Fbulk2 = gsl_spline_eval(i_data.Fbulk2, pbar, acc);
  eval.Ftemp1 = gsl_spline_eval(i_data.Ftemp1, pbar, acc);
  eval.Ftemp2 = gsl_spline_eval(i_data.Ftemp2, pbar, acc);
  
  gsl_interp_accel_free(acc);
  return eval;
};