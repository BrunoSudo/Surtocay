#include "output.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

void generate_output(vector<vector<vector<double>>> grid, vector<double> p_mags, vector<double> p_thetas, vector<double> p_phis, int id) {

  ofstream output;
  string filename = "Observables_for_PDGid_" + to_string(id) + ".dat";
  output.open(filename, ios::out);

  output << fixed << setprecision(10);

  //write header with angles
  output << "     Angles theta: ";
  for (int i = 0; i < grid[0].size(); ++i) {
    output << p_thetas[i] << " ";
  }
  output << "\n";

  output << "Mags.\\Angles phi: ";
  for (int i = 0; i < grid[0].size(); ++i) {
    output << p_phis[i] << " ";
  }
  output << "\n";

  //writes magnitude and each line of the grid
  for (int i = 0; i < grid.size(); ++i) {
    output << p_mags[i] << " ";

    for (int j = 0; j < grid[0].size(); ++j) {
      for (int k = 0; k < grid[0][0].size(); ++k) {
        output << grid[i][j][k] << " ";
      }
    }
    output << "\n";
  } 
  output.close();

  cout << "Observables_for_PDGid_" << id << ".dat has been written." << endl;
}