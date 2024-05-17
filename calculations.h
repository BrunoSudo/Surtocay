#ifndef CALCULATIONS_H 
#define CALCULATIONS_H

#include "read_MUSIC.h"
#include "interpolate_fastreso.h"

using namespace std;





double contraction(vector<double> A_a, vector<double> Bb);

vector<vector<double>> lower_index_pi(vector<vector<double>> pi);

double contraction_matrix(vector<double> A_a, vector<double> B_b, vector<vector<double> > Cab);

void Energy_Density( vector<double> p_transverses, vector<double> p_rapidities, vector<double> p_azimuthals);

double F(Evaluated_data eval,  vector<double>p, vector<double>da, double momentum_pi_term, double bulkPi_term);

double G(Evaluated_data eval, vector<double>da, vector<double> u, double momentum_pi_term, double bulkPi_term, double Epbar);

double H(Evaluated_data eval, double momentum_freezeout_pi_term);


double Expression(FO_surf surf, Evaluated_data eval, double p_x, double p_y, double p_z, double pbar_squared, double Epbar, double c_s);

#endif