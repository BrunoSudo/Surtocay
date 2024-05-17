//interpolate_fastreso.h

#ifndef INTERPOLATE_FASTRESO_H 
#define INTERPOLATE_FASTRESO_H

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

//For storing the interpolated functions
struct Interpolated_data {
  double m;
  gsl_spline* Feq1;
  gsl_spline* Feq2;
  gsl_spline* Fshear1;
  gsl_spline* Fshear2;
  gsl_spline* Fshear3;
  gsl_spline* Fbulk1;
  gsl_spline* Fbulk2;
  gsl_spline* Ftemp1;
  gsl_spline* Ftemp2;
  gsl_spline* Fvel1;
  gsl_spline* Fvel2;
  gsl_spline* Fvel3;
};

//For storing the evaluated functions in a point of Epbar
struct Evaluated_data {
  double m;
  double Feq1 , Feq2;
  double Fshear1 , Fshear2 , Fshear3;
  double Fbulk1 , Fbulk2;
  double Ftemp1 , Ftemp2;
  double Fvel1 , Fvel2 , Fvel3;
};


Interpolated_data TakeInput_FastReso(int id);

Evaluated_data Evaluate_interpolation(double Epbar,Interpolated_data i_data);

#endif