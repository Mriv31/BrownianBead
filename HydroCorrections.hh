#ifndef HYDROCOR_H
#define HYDROCOR_H

#include <math.h>
#include <ostream>
#include <iostream>


// Series for computation of transverse translationnal and rotationnal coefficients. Physica A 189 (1992) 447-477  G.S. Perkins and R.B. Jones
extern float Pn_tt[21];
extern float Pn_rr[21];
extern float Pn_tr[21];
extern float eta_water; // viscosity of water g.nm/s far from the surface

struct gamma_coeffs
{
  float z_tt;  //perpendicular to the wall, translation
  float z_rr;  //rotation along the axis perpendicular to the wall
  float xy_tt; //translation parallel to the wall
  float xy_tr; //coupling between translation and rotation
  float xy_rr; //rotation along an axis parallel to the wall
} ;


class Hydrocorr
{
private:
  const float R; // diameter of the bead, nanometers
  const float eta; // viscosity of fluid g.nm/s
  float d; // distance between the center of the sphere and the wall
  // Coefficients of hydrodynamic drag
  struct gamma_coeffs coeffs;
  int compute_coefficients(void);
  float compute_z_tt(int nmax = 30);
  float compute_z_rr(int nmax = 30);
  float compute_xy_tt(void);
  float compute_xy_tr(void);
  float compute_xy_rr(void);

public:
  Hydrocorr():R(500),eta(eta_water),d(540) {compute_coefficients();}
  Hydrocorr(float Ra, float etaa, float da):R(Ra),eta(etaa),d(da) {compute_coefficients();}
  Hydrocorr(float Ra, float da):R(Ra),eta(eta_water),d(da) {compute_coefficients();}
  int update_distance_and_coeffs(float di) {d = di; return compute_coefficients();}
  struct gamma_coeffs get_coeffs(void) const {return coeffs;} // return a copy of the structyre hydro_coeffs
};

std::ostream& operator << (std::ostream& out, const Hydrocorr& H);

#endif
