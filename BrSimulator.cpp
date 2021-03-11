#ifndef BRSIMULATOR_C
#define BRSIMULATOR_C

#include "BrSimulator.hh"
#ifndef M_PI
#define M_PI 3.1415926535
#endif
using namespace std::chrono;



int BrSimulator::init_simulation()
{

  n = static_cast<int>(Tt/dt);
  x = new double[n];
  y = new double[n];
  z = new double[n];

  ii = 0;
  compute_WLC_effective_persistence_length();
  l = WLC_inversed(force)*Nb;

  mol << 0,cos(angle),-sin(angle);
  pole << 0,1,0;
  mag_field << 0,1,0;
  y[0] = -R*cos(angle); //
  z[0] = l + R*sin(angle); //
  x[0] = 0;

  Ct = 300000*force*alpha; // pour un gap de 300 microns

  exported = 0;

  sqrtdt = sqrt(dt);

  Hydrocorr H(R,z[0]); //compute hydrodynamic coefficients at the current situation
  // std::cout << H << std::endl;
  gam = H.get_coeffs(); // coppy coefficients in the structure
  if (no_coupling) gam.xy_tr = 0;
  compute_coupling_matrices();//compute coupling matrix between rotation and translation


  sigmaz = sqrt(2*kbT*gam.z_tt);
  sigmat0 = sqrt(2*kbT*D(0));
  sigmat1 = sqrt(2*kbT*D(1));
  sigmath = sqrt(2*kbT*gam.z_rr);

  distribution=new std::normal_distribution<double>;
  std::normal_distribution<double> dtemp(0,1);
  distribution->param(dtemp.param());
  generator=new std::mt19937;

  return 0;

}

int BrSimulator::iterate()
{
  if ( ii<0 || ii>=n) return 0;



  // compute force applied by the molecule
  fm = WLC(); // norm of the force
  Eigen::Vector3d fd(-mol(0)-x[ii]/R,-mol(1)-y[ii]/R,-mol(2)-z[ii]/R); //direction of the force
  fd = fd * fm / fd.norm();  // renormalize vector so that its norm is fm;
  // std::cout << fd << std::endl;

  Eigen::Vector3d Mom(mol);
  Mom = R*Mom.cross(fd); // mol needs to be multiplied by R
  Eigen::Vector3d Mommag(pole);
  Mommag = Mommag.cross(mag_field);
  Mom = Mom + Mommag*Ct;

  Eigen::Vector3d omega;

  double Fyt0 = P_1(0,0)*fd(1) + P_1(0,1)*Mom(0); //coupling between vy and omegax
  double Fyt1 = P_1(1,0)*fd(1) + P_1(1,1)*Mom(0);

  double Fxt0 = P(0,0)*fd(0) + P(0,1)*Mom(1);//coupling between vx and omegay
  double Fxt1 = P(1,0)*fd(0) + P(1,1)*Mom(1);

  double vxt0 =  Fxt0 / D(0) +
       sigmat0 / sqrtdt * (*distribution)(*generator) / D(0);

  double vyt0 =  Fyt0 / D(0) +
          sigmat0 / sqrtdt * (*distribution)(*generator) / D(0);

  double vxt1 =  Fxt1 / D(1) +
               sigmat1 / sqrtdt * (*distribution)(*generator) / D(1);

  double vyt1 =  Fyt1 / D(1) +
                  sigmat1 / sqrtdt * (*distribution)(*generator) / D(1);

  double vz =  (force+fd(2)) / gam.z_tt +
                  sigmaz / sqrtdt * (*distribution)(*generator) / gam.z_tt;

// std::cout << vxt0 << " " << vxt1 <<" " << vyt0 << " " <<vyt1 <<  " " << vz << std::endl;



  double vx = P_1(0,0)*vxt0+P_1(0,1)*vxt1;
  double vy = P(0,0)*vyt0+P(0,1)*vyt1;
  if (Three_D)
  {
    omega(1) = P_1(1,0)*vxt0+P_1(1,1)*vxt1;
    omega(2) =  (Mom(2)) / gam.z_rr +
                                  sigmath / sqrtdt * (*distribution)(*generator) / gam.z_rr;
  }
  else
  {
    omega(1) = 0;
    omega(2) = 0;
    }
  omega(0) = P(1,0)*vyt0+P(1,1)*vyt1;
  if (Three_D)   x[ii+1] = x[ii] + dt*vx;
  y[ii+1] = y[ii] + dt*vy;
  z[ii+1] = z[ii] + dt*vz;

  pole = pole + dt*omega.cross(pole);
  pole = pole/pole.norm();

  mol = mol + dt*omega.cross(mol);
  mol = mol/mol.norm();

  //correct pole to preserve angle

  double cp = mol.dot(pole);
  double alpha0 = acos(cp);

  pole = pole*cos(alpha0-angle)+sin(alpha0-angle)/sin(alpha0) * mol - sin(alpha0-angle)*cp/sin(alpha0)*pole; //correct pole to avoid long term drift
  omega(0) = R*mol(0) + x[ii+1];
  omega(1) = R*mol(1) + y[ii+1];
  omega(2) = R*mol(2) + z[ii+1];

  l = omega.norm();
  //std::cout << pole << std::endl;
  // std::cout << x[ii] << " " << y[ii] <<" " << z[ii] << " " <<fm <<  " " << fd(2) << " " << force << "" << l << std::endl;





  ii++;
  if (ii == n) return 1;

  return 0;

}

double BrSimulator::WLC() //takes length returns force
{
    double ln = l/Nb;
    double res =    (4.114/Lpeff) * (0.25/(1-ln/L0)/(1-ln/L0) - 0.25 + ln/L0);
    double poww = ln*ln/L0/L0;
    for (int i= 2;i<8;i++)
    {
        res += (4.114/Lpeff) *(alpha_WLC[i-2] *poww);
        poww = poww * ln/L0;
      }
    return res;
}

double WLC_GSL(double ln, void *params) //takes length returns force
{
    struct WLC_params *p = (struct WLC_params *) params;
    double Fapp = p->Fapp;
    double Lpeff  = p->Lpeff;
    double res =    (4.114/Lpeff) * (0.25/(1-ln/L0)/(1-ln/L0) - 0.25 + ln/L0);
    double poww = ln*ln/L0/L0;
    for (int i= 2;i<8;i++)
    {
        res += (4.114/Lpeff) *(alpha_WLC[i-2] *poww);
        poww = poww * ln/L0;
      }
    return res-Fapp;
}


double BrSimulator::WLC_inversed(double ff) // takes force return length
{

  double x_lo = 0.0, x_hi = L0-0.001;
  double r;
  gsl_function F;
  const gsl_root_fsolver_type *T;
  int iter = 0;
  int max_iter = 1e4;
  gsl_root_fsolver *s;
  int status;
  struct WLC_params params = {ff,Lpeff};

  F.function = &WLC_GSL;
  F.params = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);


  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,
                                       0, 0.001);

    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);

  return r;
}

int BrSimulator::compute_WLC_effective_persistence_length()
{

    double p_inf = 51.7;
    double L_c = L0*Nb;
    double a =2.78;
    Lpeff = p_inf/(1+a*p_inf/L_c);
    return 0;

}

int BrSimulator::compute_coupling_matrices()
{
    Eigen::Matrix2d CM;
    CM << gam.xy_tt, gam.xy_tr,
          gam.xy_tr, gam. xy_rr;
    Eigen::EigenSolver<Eigen::Matrix2d> eigensolver(CM);
      if (eigensolver.info() != Eigen::Success) abort();

    P = eigensolver.eigenvectors().real();
    D =  eigensolver.eigenvalues().real();
    P_1 = P.inverse();
    // std::cout<< P(0,0) <<" "<<P(0,1)<< " " <<P(1,0) << " " <<P(1,1) << std::endl;
    // std::cout<< P_1(0,0) <<" "<<P_1(0,1)<< " " <<P_1(1,0) << " " <<P_1(1,1) << std::endl;

    return 0;
}

int main () {
  int n = 0; // contains size of simulation

 float *x = NULL, *y = NULL, *z = NULL;
 BrSimulator B(500,10,1.57,0.1,1e-8,0,25200,1,0);
 while (B.iterate() == 0)
 {
 ;
 }
 n = B.export_simulation(&x,&y,&z);
// while (B.iterate() == 0)
// {
// ;
// }


}






#endif
