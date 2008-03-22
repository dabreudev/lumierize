#include "modulos.h"
void NumericalDerivCovars_STY_p_M(int n,double *magn,double *z,double *par,double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf);

int main()
{
  int n;
  double *magn;
  double *z;
  double *par;
  double *sigpar;

  double mlim;
  struct cosmo_param cosmo;
  struct Schlf_M lf;
  
  n = 100;
  mlim = 10;
  cosmo.H0 = 75;
  cosmo.q0 = 1;
  magn = vector_d(100);
  z = vector_d(100);
  
  par = vector_d(3);
  sigpar = vector_d(3);
  
  par[0] = 0;
  par[1] = 1;
  par[2] = 2;
  
 NumericalDerivCovars_STY_p_M(n,magn,z,par,sigpar,mlim,cosmo,&lf);
 
 return 0; 
}
