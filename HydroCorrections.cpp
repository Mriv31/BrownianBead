#ifndef HYDROCOR_C
#define HYDROCOR_C

#include "HydroCorrections.hh"

float Pn_tt[21]={1,0.5625,0.3164036,0.05297852,0.135292,0.1979151,0.0004110932,0.06747958,0.1289003,-0.01470387,0.08618274,0.07600641,-0.01460087,0.07790617,0.03443769,0.01174160,0.05623434,0.02651131,0.01847407,0.03826369,0.02841849};
float Pn_tr[21]={0,0,0,0,-0.09375,0.03515625,0.01977539,-0.08702087,0.04004002,-0.001967847,-0.06136864,0.04651307,-0.02224856,-0.03267993,0.03247387,-0.02923364,-0.008874159,0.01258875,-0.01817025,-0.005859199,0.005789533};
float Pn_rr[21]={1,0,0,0.3125,0,0,0.15625,0,0.01171875,0.1066895,-0.004211426,0.01411057,0.08651996,-0.01874109,0.03269177,0.05358306,-0.01315024,0.03888607,0.02790070,0.002833608,0.02777366};

float eta_water = 1e-9; // viscosity of water g.nm/s far from the surface


int Hydrocorr::compute_coefficients(void)
{

  coeffs.z_tt = 6*M_PI*eta*R*compute_z_tt();
  coeffs.z_rr = 8*M_PI*eta*R*R*R*compute_z_rr();
  coeffs.xy_tr = 8*M_PI*eta*R*R*compute_xy_tr();
  coeffs.xy_tt = 6*M_PI*eta*R*compute_xy_tt();
  coeffs.xy_rr = 8*M_PI*eta*R*R*R*compute_xy_rr();
  return 0;
}

float Hydrocorr::compute_z_tt(int nmax)
{
   double res = 0, fac = 0, num = 0, den = 0, sum = 0;
   double alpha = acosh(d/R);
   for (int n = 1;n<nmax;n++)
   {
     fac = n*(n+1)/(float)((2*n-1)*(2*n+3));
     num = 2*sinh(alpha*(2*n+1))+(2*n+1)*sinh(2*alpha);
     den = 4*sinh(alpha*(n+0.5))*sinh(alpha*(n+0.5))-(2*n+1)*(2*n+1)*sinh(alpha)*sinh(alpha);
     sum += fac*(num/den-1);
     // std::cout << n << " " << sum << " " << num << std::endl;
   }
   return static_cast<float>((4/3.0)*sinh(alpha)*sum);
}

float Hydrocorr::compute_z_rr(int nmax)
{
   float alpha = acosh(d/R), sum = 0;
   for (int n = 1;n<nmax;n++)
   {
    sum += 1/pow(sinh(n*alpha),3);
   }
   return pow(sinh(alpha),3)*sum;
}

float Hydrocorr::compute_xy_tt(void)
{
  float t = R/d,sum = 1;
  int nmax = 21;
  for (int n =1;n<nmax;n++)
  {
    sum+=(Pn_tt[n]-8.0/15.0/n)*pow(t,n);
  }
  return -8/15.0 * log(1-t) + sum;
}

float Hydrocorr::compute_xy_tr(void)
{
  float t = R/d,sum = 0;
  int nmax = 21;
  for (int n =1;n<nmax;n++)
  {
    sum+=(Pn_tr[n]+1.0/10.0/n)*pow(t,n);
  }
  return 1/10.0 * log(1-t) + sum;
}

float Hydrocorr::compute_xy_rr(void)
{
  float t = R/d,sum = 1;
  int nmax = 21;
  for (int n =1;n<nmax;n++)
  {
    sum+=(Pn_rr[n]-2/5.0/n)*pow(t,n);
  }
  return -2/5.0 * log(1-t) + sum;
}

std::ostream& operator << (std::ostream& out, const Hydrocorr& H)  // << overloading in order to print a graph
{
  struct gamma_coeffs g = H.get_coeffs();
  out << "gamma_z_tt " << g.z_tt << std::endl;
  out << "gamma_z_rr " << g.z_rr << std::endl;
  out << "gamma_xy_tt " << g.xy_tt << std::endl;
  out << "gamma_xy_tr " << g.xy_tr << std::endl;
  out << "gamma_xy_rr " << g.xy_rr << std::endl;
  return out;
}





#endif
