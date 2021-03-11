#ifndef BRSIMULATOR_H
#define BRSIMULATOR_H

#include <math.h>
#include <ostream>
#include <iostream>
#include <Eigen/Eigenvalues>
#include "HydroCorrections.hh"
#include <chrono>
#include <random>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

const double kbT = 4.11; //pn.nm
const double alpha_WLC[6]={-0.5164228,-2.737418,16.07497,-38.87607,39.49944,-14.17718};
const double L0 = 0.34; //parameter of WLC


class BrSimulator
{
private:
  // parameters of the simulation
  const double R; // radius of the bead
  const double force; //magnetic force applied on the bead (pN)
  const double angle; //angle between anchoring point and main magnetic point
  const double Tt; //whole time simulation (second)
  const double dt; //time step
  const double alpha; // magnetisation of the bead
  const double Nb; //number of bases in DNA
  const int Three_D; // if false no movement in x
  const int no_coupling; //if true no hydrodynamic coupling





  // current parameters of the simulation
  double *x; // X-position of the bead
  double *y; // Y- position of the bead
  double *z; // Z-position of the bead

  double sqrtdt;
  double sigmat0;
  double sigmat1;
  double sigmath;
  double sigmaz;

  int exported;

  //gamma Coefficients
  struct gamma_coeffs gam;

  Eigen::Vector3d mol; //current position of the molecule compared to the center of the bead
  Eigen::Vector3d pole; //current position of the magnetic pole compared to the center of the bead
  Eigen::Vector3d mag_field;
  double fm; //current force applied by the molecule
  double Ct; //torque
  int n; //number of steps in simulation
  double l; //current length of the molecule
  int ii; //current point
  double Lpeff; //effective persitence length for worm like chain

  Eigen::Matrix2d P;
  Eigen::Matrix2d P_1;

  Eigen::Vector2d D;


  int init_simulation();
  int compute_coupling_matrices();
  int compute_WLC_effective_persistence_length();



  double WLC(); //return the current force applied by the molecule on the bead
  double WLC_inversed(double force);

  std::mt19937 * generator;
  std::normal_distribution<double>* distribution;


public:
  BrSimulator():R(500),force(10),angle(1.57),Tt(1),dt(1e-7),alpha(0.01),Nb(120),Three_D(1), no_coupling(0){init_simulation();}
  BrSimulator(double Ra, double fra, double aa, double Ta, double dta, double aaa, double Nba, int Threea, int noa):R(Ra),force(fra),angle(aa),Tt(Ta),dt(dta),alpha(aaa),Nb(Nba),Three_D(Threea),no_coupling(noa){init_simulation();}
  int iterate();
  int perform_whole_simulation();
  int export_simulation(double **xe, double **ye, double **ze) {*xe=x; *ye=y; *ze=z; exported = 1; return n;}
  int export_simulation(float **xe, float **ye, float **ze) {  *xe = new float[n];*ye = new float[n];*ze = new float[n];for  (int i = 0;i<n;i++) {(*xe)[i]=(float)(x[i]); (*ye)[i]=(float)(y[i]); (*ze)[i]=(float)(z[i]);} return n;}

  ~BrSimulator(){delete generator; delete distribution; if (exported == 0) {delete x; delete y; delete z;}};

};

struct WLC_params
{
  double Fapp;
  double Lpeff;

} ;
double WLC_GSL(double ln, void *params); //takes length returns force






#endif
