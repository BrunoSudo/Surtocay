#include "interpolate_fastreso.h"
#include "read_MUSIC.h"
#include "calculations.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#define PI M_PI

using namespace std;

int main() {
  
  //Momenta considered for calculations

  vector<double> p_transverses = {0,3,7, 10}; //transverse momentum
  vector<double> p_rapidities = {-4,2,0,2,4};  //rapidity
  vector<double> p_azimuthals =   {0,3.14, 2*3.14}; //azimuthal angle

  cout << p_transverses.size() << endl;
  
  //Calculate the energy density for each momentum vector inputed
  for (int i = 0; i < p_transverses.size(); i++) { // Magnitude of the momentum
    for (int j = 0; j < p_azimuthals.size(); j++) {  //Angle of the momentum
      for (int k = 0; k < p_rapidities.size(); k++) {
        
        //Calculate the magnitude of the momentum squared
        double pbar_squared = pow(p_transverses[i],2) + pow(p_rapidities[j],2);

        double pbar_211 = sqrt(pbar_squared);
        double pbar_321 = sqrt(pbar_squared);
        double pbar_2212 = sqrt(pbar_squared);
        double pbar_3122 = sqrt(pbar_squared);
        double pbar_3312 = sqrt(pbar_squared);
        double pbar_3334 = sqrt(pbar_squared);

        cout << "pbar_211: " << pbar_211 << endl;


      }
    }
  }
  //Energy_Density(p_transverses, p_rapidities, p_azimuthals);


  

  return 0;
}
