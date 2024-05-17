#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include "calculations.h"
#include "read_MUSIC.h"
#include "interpolate_fastreso.h"
#include "output.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <vector>
using namespace std;

#define PI M_PI


/////////////////////////////////////////////////////////////////////////////////
//Operations
double contraction(vector<double> A_a, vector<double> Bb) {
  double dot = 0;
  for (int i = 0; i < 4; i++) {
    dot += A_a[i]*Bb[i];
  }
  return dot;
}

vector<vector<double>> lower_index_pi(vector<vector<double>> pi, double tau, int N_index) {

  //metric tensor g_{\mu\nu}
  vector<vector<double>> g = {{1.0, 0.0, 0.0, 0.0}, 
                              {0.0, -1.0, 0.0, 0.0}, 
                              {0.0, 0.0, -1.0, 0.0}, 
                              {0.0, 0.0, 0.0, -(tau*tau)}};
  
  //result tensor
  vector<vector<double>> result(4, vector<double>(4, 0.0));

  //Lower first index
  for (int mu = 0; mu < 4; mu++) {
    for (int nu = 0; nu < 4; nu++) {
      for (int alpha = 0; alpha < 4; alpha++) {
        result[mu][nu] += g[mu][alpha] * pi[alpha][nu];
      }
    }
  }

  //Lower only first index 
  if (N_index == 1) {
    return result;
  }


  //lower both indices
  else {
    // Copy the result to pi for the second lowering operation
    pi = result;

    // Reset result to zero
    for (int mu = 0; mu < 4; mu++) {
      for (int nu = 0; nu < 4; nu++) {
        result[mu][nu] = 0.0;
      }
    }

    // Lower the second index
    for (int mu = 0; mu < 4; mu++) {
      for (int nu = 0; nu < 4; nu++) {
        for (int alpha = 0; alpha < 4; alpha++) {
          result[mu][nu] += g[nu][alpha] * pi[mu][alpha];
        }
      }
    }
    
    return result;
  }
}

double contraction_matrix(vector<double> A_a, vector<double> B_b, vector<vector<double>> Cab) {
  double dot = 0;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      dot += A_a[i]*Cab[i][j]*B_b[j];
    }
  }
  return dot;
}

///////////////////////////////////////////////////////////////////////////////////////////
////Calculations of the functions for the integral

//Calculate each step of the integral while reading each line of the surface.dat file
void Energy_Density( vector<double> p_transverses, vector<double> p_azimuthals, vector<double> p_rapidities) {
  //Make 3D grid for the momentums that the expression will be evaluated for for each particle
  vector<vector<vector<double>>> grid_id_211(p_transverses.size(), vector<vector<double>>(p_azimuthals.size(), vector<double>(p_rapidities.size(), 0)));
  vector<vector<vector<double>>> grid_id_321(p_transverses.size(), vector<vector<double>>(p_azimuthals.size(), vector<double>(p_rapidities.size(), 0)));
  vector<vector<vector<double>>> grid_id_2212(p_transverses.size(), vector<vector<double>>(p_azimuthals.size(), vector<double>(p_rapidities.size(), 0)));
  vector<vector<vector<double>>> grid_id_3122(p_transverses.size(), vector<vector<double>>(p_azimuthals.size(), vector<double>(p_rapidities.size(), 0)));
  vector<vector<vector<double>>> grid_id_3312(p_transverses.size(), vector<vector<double>>(p_azimuthals.size(), vector<double>(p_rapidities.size(), 0)));
  vector<vector<vector<double>>> grid_id_3334(p_transverses.size(), vector<vector<double>>(p_azimuthals.size(), vector<double>(p_rapidities.size(), 0)));

  //Interpolated data for each particle
  Interpolated_data id_211 = TakeInput_FastReso(211);
  Interpolated_data id_321 = TakeInput_FastReso(321);
  Interpolated_data id_2212 = TakeInput_FastReso(2212);
  Interpolated_data id_3122 = TakeInput_FastReso(3122);
  Interpolated_data id_3312 = TakeInput_FastReso(3312);
  Interpolated_data id_3334 = TakeInput_FastReso(3334);
  

  string surface_filename = "surface.dat";
  string path_ = "../.."; 


  ostringstream surfdat_stream;
  surfdat_stream << path_ << "/" << surface_filename;
  ifstream surfdat;
  surfdat.open(surfdat_stream.str().c_str(), std::ios::binary);
  
  if (!surfdat.is_open()) {
    cout << "Error: file " << surfdat_stream.str() << "could not be opened." << endl;
    exit(1);
  }

  
  //Go through the surface.dat file
  int line = 0;

  while (!surfdat.eof()) {
    FO_surf surf = read_MUSIC_line(surfdat);

    //break the loop if the values in the file are zero
    if (surf.Edec == 0) {
      break;
    }

    double c_s = 1./3.;

    //Calculate the energy density for each momentum vector inputed
    for (int i = 0; i < p_transverses.size(); i++) { // Magnitude of the momentum
      for (int j = 0; j < p_azimuthals.size(); j++) {  //Angle of the momentum
        for (int k = 0; k < p_rapidities.size(); k++) {

          //Calculate the magnitude of the momentum squared
          double pbar_squared = pow(p_transverses[i],2) + pow(p_rapidities[j],2);
          
          //Calculate the momentum components
          double p_x = p_transverses[i]*cos(p_azimuthals[j]);
          double p_y = p_transverses[i]*sin(p_azimuthals[j]); 

          //Calculate the longitudinal momentum
          double p_z_211 = sqrt(pow(1.40000e-01, 2) * pbar_squared) * (exp(2*p_rapidities[k]) -1) / (exp(2*p_rapidities[k]) + 1);
          double p_z_321 = sqrt(pow(4.94000e-01, 2) * pbar_squared) * (exp(2*p_rapidities[k]) -1) / (exp(2*p_rapidities[k]) + 1);
          double p_z_2212 = sqrt(pow(9.38300e-01, 2) * pbar_squared) * (exp(2*p_rapidities[k]) -1) / (exp(2*p_rapidities[k]) + 1);
          double p_z_3122 = sqrt(pow(1.11568e+00, 2) * pbar_squared) * (exp(2*p_rapidities[k]) -1) / (exp(2*p_rapidities[k]) + 1);
          double p_z_3312 = sqrt(pow(1.32131e+00, 2) * pbar_squared) * (exp(2*p_rapidities[k]) -1) / (exp(2*p_rapidities[k]) + 1);
          double p_z_3334 = sqrt(pow(1.67245e+00, 2) * pbar_squared) * (exp(2*p_rapidities[k]) -1) / (exp(2*p_rapidities[k]) + 1);



          //Calculate the magnitude of the momentum
          double pbar_211 = sqrt(pbar_squared);
          double pbar_321 = sqrt(pbar_squared);
          double pbar_2212 = sqrt(pbar_squared);
          double pbar_3122 = sqrt(pbar_squared);
          double pbar_3312 = sqrt(pbar_squared);
          double pbar_3334 = sqrt(pbar_squared);

          cout << "pbar_211: " << pbar_211 << endl;
          cout << "pbar_321: " << pbar_321 << endl;
          cout << "pbar_2212: " << pbar_2212 << endl;
          cout << "pbar_3122: " << pbar_3122 << endl;
          cout << "pbar_3312: " << pbar_3312 << endl;
          cout << "pbar_3334: " << pbar_3334 << endl;

          //Evaluate the interpolated data for each particle
          Evaluated_data eval_id_211 = Evaluate_interpolation(pbar_211, id_211);
          Evaluated_data eval_id_321 = Evaluate_interpolation(pbar_321, id_321);
          Evaluated_data eval_id_2212 = Evaluate_interpolation(pbar_2212, id_2212);
          Evaluated_data eval_id_3122 = Evaluate_interpolation(pbar_3122, id_3122);
          Evaluated_data eval_id_3312 = Evaluate_interpolation(pbar_3312, id_3312);
          Evaluated_data eval_id_3334 = Evaluate_interpolation(pbar_3334, id_3334);

          //Calculate the integral for each particle
          double step_id_211 = Expression(surf, eval_id_211, p_x, p_y, p_z_211, pbar_squared, pbar_211, c_s);
          double step_id_321 = Expression(surf, eval_id_321, p_x, p_y, p_z_321, pbar_squared, pbar_321, c_s);
          double step_id_2212 = Expression(surf, eval_id_2212, p_x, p_y, p_z_2212, pbar_squared, pbar_2212, c_s);
          double step_id_3122 = Expression(surf, eval_id_3122, p_x, p_y, p_z_3122, pbar_squared, pbar_3122, c_s);
          double step_id_3312 = Expression(surf, eval_id_3312, p_x, p_y, p_z_3312, pbar_squared, pbar_3312, c_s);
          double step_id_3334 = Expression(surf, eval_id_3334, p_x, p_y, p_z_3334, pbar_squared, pbar_3334, c_s);
          
          //Add the step to the grid
          grid_id_211[i][j][k] += step_id_211; 
          grid_id_321[i][j][k] += step_id_321;
          grid_id_2212[i][j][k] += step_id_2212;
          grid_id_3122[i][j][k] += step_id_3122;
          grid_id_3312[i][j][k] += step_id_3312;
          grid_id_3334[i][j][k] += step_id_3334;
        }
      }
    }
  line++;
  }

  surfdat.close();
  //multiply by term in front of integral considering the spin factor of each particle
    for (int i = 0; i < p_transverses.size(); i++) { // Transverse momentum
      for (int j = 0; j < p_azimuthals.size(); j++) {  // Azimuthal angle of the momentum
        for (int k = 0; k < p_rapidities.size(); k++) {  //Pseudo-rapidity of the momentum
          grid_id_211[i][j][k] *= 1/pow(2*PI,3);  
          grid_id_321[i][j][k] *= 2/pow(2*PI,3);
          grid_id_2212[i][j][k] *= 2/pow(2*PI,3);
          grid_id_3122[i][j][k] *= 1/pow(2*PI,3);
          grid_id_3312[i][j][k] *= 2/pow(2*PI,3);
          grid_id_3334[i][j][k] *= 4/pow(2*PI,3);
        }
      }
    }

  generate_output(grid_id_211, p_transverses, p_azimuthals, p_rapidities, 211);
  generate_output(grid_id_321, p_transverses, p_azimuthals, p_rapidities, 321);
  generate_output(grid_id_2212, p_transverses, p_azimuthals, p_rapidities, 2212);
  generate_output(grid_id_3122, p_transverses, p_azimuthals, p_rapidities, 3122);
  generate_output(grid_id_3312, p_transverses, p_azimuthals, p_rapidities, 3312);
  generate_output(grid_id_3334, p_transverses, p_azimuthals, p_rapidities, 3334);
}




double F(Evaluated_data eval, vector<double>p, vector<double>da, double momentum_pi_term, double bulkPi_term) {
  double term1 = eval.Feq1;
  
  double term2 = eval.Fshear1*momentum_pi_term;

  double term3 = eval.Fbulk1*bulkPi_term;

  double term4 = contraction(p, da);

  double f = (term1 + term2 + term3)*term4;

  return f;
}

double G(Evaluated_data eval, vector<double>da, vector<double> u, double momentum_pi_term, double bulkPi_term, double Epbar) {
  double term1 = eval.Feq2 - eval.Feq1;
  
  double term2 = (eval.Fbulk2 - eval.Fbulk1)*bulkPi_term;

  double term3 = (eval.Fshear3 - eval.Fshear1)*momentum_pi_term;

  double term4 = contraction(da, u);

  double g = (term1 + term2 + term3)* Epbar * term4;

  return g;
}

double H(Evaluated_data eval, double momentum_freezeout_pi_term) {
  double term1 = (eval.Fshear2 - eval.Fshear1)*(2./5.);

  double term2 = momentum_freezeout_pi_term;

  double h = term1*term2;
  return h;
}

//Calculates the entire expression inside the integral
double Expression(FO_surf surf, Evaluated_data eval, double p_x, double p_y, double p_z, double pbar_squared, double Epbar, double c_s) {
  //decided value of b_{\Pi}
  double b_Pi = 1./15.;

  // Epbar from pbar
  Epbar = sqrt(pbar_squared + eval.m*eval.m);

  //make vectors for the freeze-out surface, momentum and \pi
  vector<double> p = {eval.m*eval.m , p_x, p_y, p_z};
  vector<double> u = {surf.u0, surf.u1, surf.u2, surf.u3};
  vector<double> da = {surf.da0, surf.da1, surf.da2, surf.da3};
  vector<vector<double>> pi = {{surf.pi00, surf.pi01, surf.pi02, surf.pi03},
                               {surf.pi01, surf.pi11, surf.pi12, surf.pi13},
                               {surf.pi02, surf.pi12, surf.pi22, surf.pi23},
                               {surf.pi03, surf.pi13, surf.pi23, surf.pi33}};
  
  //{p^{\nu}p^{\mu}(\pi_{\nu})_{\mu}
  vector<vector<double>> pi_lower_mu_lower_nu = lower_index_pi(pi, surf.tau, 2);

  double dot_pxpxpi = contraction_matrix(p, p, pi_lower_mu_lower_nu);

  double denominator = (2*(surf.Edec + surf.Pdec)*surf.Tdec*surf.Tdec);

  double momentum_pi_term = dot_pxpxpi/denominator;

  //{p^{\nu}da_{\mu}(\pi_{\mu})^{\nu}
  vector<vector<double>> pi_lower_mu_upper_nu = lower_index_pi(pi, surf.tau, 1);

  double dot_daxpxpi = contraction_matrix(da, p, pi_lower_mu_upper_nu);

  double momentum_freezeout_pi_term = dot_daxpxpi/denominator;
  

  
  //- \frac{b_{\Pi}\Pi}{\left(\left(\frac{1}{3} - c_{s}^2\right)^2\right) \cdot (E_{\text{dec}} + P_{\text{dec}})}
  double bulkPi_term = (-1)*b_Pi/((pow(1./3.-pow(1,2),2) * (surf.Edec + surf.Pdec))); //missing the \Pi

  double f = F(eval, p, da, momentum_pi_term, bulkPi_term);
  double g = G(eval, da, u, momentum_pi_term, bulkPi_term, Epbar);
  double h = H(eval, momentum_freezeout_pi_term);

  double integral = (f + g + h);

  
  
  return integral;
}