// read_music.h

#ifndef READ_MUSIC_H
#define READ_MUSIC_H


#include <vector>
#include <string>

typedef struct {
    float tau, xpt, ypt, eta;
    float da0, da1, da2, da3;
    float u0, u1, u2, u3;
    float Edec, Tdec, Pdec;
    float Bn, muB, muS, muQ;
    float pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33;
    float bulkPi;
    float qmu0, qmu1, qmu2, qmu3;
    std::vector<float> particle_mu_PCE;
} FO_surf;

FO_surf read_MUSIC_line(std::ifstream& surfdat);

void print_FO_surf(const std::vector<FO_surf>& surf_ptr, const std::string& filename);

double get_velocity_of_sound_sq(double Edec);

#endif