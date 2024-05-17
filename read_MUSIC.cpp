#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "read_MUSIC.h"
#include "../eos_hotQCD.h"

#define  hbarC 0.197327053

using namespace std;


void print_FO_surf(const std::vector<FO_surf>& surf_ptr, const std::string& filename) {
  ofstream file(filename);

  // Print the headers
  file << "tau\txpt\typt\teta_s\tdat\tdax\tday\tdan\tux\tuy\tun\tut\tE\tP\tT\tmuB\tmuS\tmuC\tpi00\tpi01\tpi02\tpi03\tpi11\tpi12\tpi13\tpi22\tpi23\tpi33\tbulkPi\tBn\tqmu0\tqmu1\tqmu2\tqmu3\n";

  // Print the data
  for (const auto& surf : surf_ptr) {
    file << surf.tau << "\t" << surf.xpt << "\t" << surf.ypt << "\t" << surf.eta << "\t" << surf.da0 << "\t" << surf.da1 << "\t" << surf.da2 << "\t" << surf.da3 << "\t" << surf.u0 << "\t" << surf.u1 << "\t" << surf.u2 << "\t" << surf.u3 << "\t" << surf.Edec << "\t" << surf.Pdec << "\t" << surf.Tdec << "\t" << surf.muB << "\t" << surf.muS << "\t" << surf.muQ << "\t" << surf.pi00 << "\t" << surf.pi01 << "\t" << surf.pi02 << "\t" << surf.pi03 << "\t" << surf.pi11 << "\t" << surf.pi12 << "\t" << surf.pi13 << "\t" << surf.pi22 << "\t" << surf.pi23 << "\t" << surf.pi33 << "\t" << surf.bulkPi << "\t" << surf.Bn << "\t" << surf.qmu0 << "\t" << surf.qmu1 << "\t" << surf.qmu2 << "\t" << surf.qmu3 << "\n";
  }

  file.close();
}

FO_surf read_MUSIC_line(std::ifstream& surfdat) {

  FO_surf surf_elem;
  float array[34];
  for (int i = 0; i < 34; i++) {
    float temp = 0.;
    surfdat.read(reinterpret_cast<char*>(&temp), sizeof(float));
    array[i] = temp;
  }

  surf_elem.tau = array[0];
  surf_elem.xpt = array[1];
  surf_elem.ypt = array[2];
  surf_elem.eta = array[3];
  surf_elem.da0 = array[4];
  surf_elem.da1 = array[5];
  surf_elem.da2 = array[6];
  surf_elem.da3 = array[7];
  surf_elem.u0  = array[8];
  surf_elem.u1  = array[9];
  surf_elem.u2  = array[10];
  surf_elem.u3  = array[11];

  surf_elem.Edec = array[12]*hbarC;
  surf_elem.Tdec = array[13]*hbarC;
  surf_elem.muB  = array[14]*hbarC;
  surf_elem.muS  = array[15]*hbarC;
  surf_elem.muQ  = array[16]*hbarC;
  surf_elem.Pdec = array[17]*surf_elem.Tdec - surf_elem.Edec;

  surf_elem.pi00 = array[18]*hbarC;  // GeV/fm^3
  surf_elem.pi01 = array[19]*hbarC;  // GeV/fm^3
  surf_elem.pi02 = array[20]*hbarC;  // GeV/fm^3
  surf_elem.pi03 = array[21]*hbarC;  // GeV/fm^3
  surf_elem.pi11 = array[22]*hbarC;  // GeV/fm^3
  surf_elem.pi12 = array[23]*hbarC;  // GeV/fm^3
  surf_elem.pi13 = array[24]*hbarC;  // GeV/fm^3
  surf_elem.pi22 = array[25]*hbarC;  // GeV/fm^3
  surf_elem.pi23 = array[26]*hbarC;  // GeV/fm^3
  surf_elem.pi33 = array[27]*hbarC;  // GeV/fm^3

  surf_elem.bulkPi = array[28]*hbarC;   // GeV/fm^3

  surf_elem.Bn   = array[29];             // 1/fm^3
  surf_elem.qmu0 = array[30];
  surf_elem.qmu1 = array[31];
  surf_elem.qmu2 = array[32];
  surf_elem.qmu3 = array[33];
    
  // Close the file after reading

  return surf_elem;
}


double get_velocity_of_sound_sq(double Edec) {
  EOS_base EOS1;
  EOS1.initialize_eos(1);
  EOS_hotQCD EOS(1);
  EOS.initialize_eos();

  double c_s = EOS.get_cs2(Edec, 0);

  return c_s;
}
