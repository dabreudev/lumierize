#include "modulos.h"
#include "gsl_hessian.h"
#include <gsl/gsl_machine.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>

#define ZMIN 0.00001
#define FTOL  1e-12
#define FTOL2 1e-6
#define FTOL3 1e-7
#define MAXITER  1000
#define MAXITER2 120
#define NBSNORMA  10
#define MAXTRIES   5
#define VERBOSE 0
#define DEBUG  0
#define DEBUG2 0
#define DEBUG3 0
#define DEBUGPLOT 0
#define CONTOURPLOT 0
#define NCONFL 20
#define TOLERR 0.07 

/* #define TOLERR 0.0001 */

/* Ahora mismo esta con el pcut definido, de modo que se aleja de valores 
   cercanos a 0. Puede afectar al calculo de las covarianzas!! */

/* Estructura para contener los parámetros */
struct Amoe_Funk_gsl_param_STY_p_M
{
  int nData;
  double *magn;
  double *z;
};

void prepareGlobalVars_STY_p_M(double *z, double *magn);
double Amoe_Funk_STY_p_M_main(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_p_M_main_gsl_multimin(const gsl_vector *x, void *params);
void   NumericalHessianCovars_STY_p_M(int n,double *magn,double *z,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf);

struct cosmo_param *_cosmo_STY_p_M;
double _mlim_STY_p_M;

int _ndata_STY_p_M;
int _iter_m_STY_p_M;
int _iter_c_STY_p_M;
int _nconfl_STY_p_M;
double *_pp_STY_p_M;
double _MLmax_STY_p_M;
double _xf_STY_p_M,_Tf_STY_p_M;
double _sigi_STY_p_M;
double _xtmp_STY_p_M;
double _zlow_STY_p_M,_zup_STY_p_M,_strrad_STY_p_M;
double *_Mabsn_STY_p_M;

int  MLA_STY_p_M(int n,double *magn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf, struct MLProcessInfo *mlinfo)
{

  double par[3];
  double sigpar[3];
  int i;
  int iter_amo;

  /* Variables to use for vvmax fitting */
  struct Steplf_M  lfvvmax;
  double minMabs,maxMabs;
  struct Schlf_M  lffit;
  double chisq;

  /* Copias globales de estas variables para que se vean en toda la subrutina */
  _ndata_STY_p_M=n;
  _cosmo_STY_p_M=&cosmo;
  _mlim_STY_p_M=mlim;
  _strrad_STY_p_M=strrad;
  _zlow_STY_p_M=zlow;
  _zup_STY_p_M=zup;

  _iter_m_STY_p_M=0;

  iter_amo=MAXITER+1;

  /* Los límites en z se reajustan */
  _zlow_STY_p_M = (_zlow_STY_p_M < ZMIN ? ZMIN : _zlow_STY_p_M);

  if(DEBUG3) printf(" iter_amo %d\n",iter_amo);

  /* VVmax computing as a first solution */
  lfvvmax.nbin=(int)(n/3);
  if(lfvvmax.nbin>30) lfvvmax.nbin=30;

  lfvvmax.magni     =vector_d(lfvvmax.nbin+1);
  lfvvmax.errmagni  =vector_d(lfvvmax.nbin+1);
  lfvvmax.lnlf      =vector_d(lfvvmax.nbin);
  lfvvmax.errlnlf   =vector_d(lfvvmax.nbin);
  lfvvmax.lf        =vector_d(lfvvmax.nbin);
  lfvvmax.errlf     =vector_d(lfvvmax.nbin);
  lfvvmax.covarlnlf =matrix_d(lfvvmax.nbin,lfvvmax.nbin);

  prepareGlobalVars_STY_p_M(z,magn); /* inicializa _Mabsn_STY_p_M */

  MinMax_d(n,_Mabsn_STY_p_M,&minMabs,&maxMabs);
  for(i=0;i<=lfvvmax.nbin;i++) lfvvmax.magni[i]=minMabs+i*(maxMabs-minMabs)/lfvvmax.nbin;
  VVmax_M(n,magn,magn,z,mlim,strrad,zlow,zup,cosmo,&lfvvmax);

  /* Fit of the VVmax solution to a Schechter function */
  FitSch2StepLF_M(lfvvmax,&lffit, &chisq);

  /* Free vvmax related things */
  free(lfvvmax.magni);
  free(lfvvmax.errmagni);
  free(lfvvmax.lnlf);
  free(lfvvmax.errlnlf);
  free(lfvvmax.lf);
  free(lfvvmax.errlf);
  free_matrix_d(lfvvmax.covarlnlf,lfvvmax.nbin,lfvvmax.nbin);

  /* Feed the initial solution*/
  par[0]=lffit.alfa;  /* Alpha */
  par[1]=lffit.Mstar;
  par[2]=log(lffit.phistar);
  sigpar[0]=10.*lffit.erralfa;
  sigpar[1]=10.*lffit.errMstar;
  sigpar[2]=10.*lffit.errphistar/lffit.phistar;

  printf(" Computing LF...\n");

  iter_amo=Amoeba_d(n,magn,z,3,par,sigpar,FTOL,MAXITER,Amoe_Funk_STY_p_M_main);
/*   printf(" FINAL par0 %.15g par1 %.15g \n",par[0],par[1]); */
  
  /* Info that will be output in mlinfo */
  _MLmax_STY_p_M=Amoe_Funk_STY_p_M_main(n,magn,z,par);
  mlinfo->nIter = iter_amo;
  mlinfo->MLmax = _MLmax_STY_p_M;

  /* Meto la solucion en la salida */
  lf->alfa=par[0];
  lf->Mstar=par[1];
  lf->phistar=exp(par[2]);

  /* Estimacion de los errores en los parametros */

  printf(" Computing errors in LF parameters...\n");

  NumericalHessianCovars_STY_p_M(n,magn,z,par,sigpar,mlim,cosmo,lf); 
  if(DEBUG) printf(" Solucion final: Mstar %.15g +/- %.15g alpha %.4g +/- %.4g\n",lf->Mstar,lf->errMstar,lf->alfa,lf->erralfa);  
  
  if(DEBUG) printf(" MLF %g\n",_MLmax_STY_p_M); 

  free(_Mabsn_STY_p_M);

  if(iter_amo>=MAXITER-1) return(2);
  return(0);
}

void prepareGlobalVars_STY_p_M(double *z, double *magn)
{
  size_t i;

   _Mabsn_STY_p_M = vector_d(_ndata_STY_p_M);
  for(i=0;i<_ndata_STY_p_M;i++) _Mabsn_STY_p_M[i]=Mag(z[i],magn[i],*_cosmo_STY_p_M);
}

double Amoe_Funk_STY_p_M_main(int n, double *x, double *y, double *p)
{
  int i;
  double logL=0.;
  struct Schlf_M lfamo;
  double Mabs;
  double log_gamma_int;
  double Mlow;
  double Lstar;
  double Llow;
  double Ntot;

  lfamo.alfa=p[0];
  lfamo.Mstar=p[1];
  lfamo.phistar=exp(p[2]);

  logL=0.;

  Lstar=pow(10.,-0.4*lfamo.Mstar);

  for(i=0;i<_ndata_STY_p_M;i++)
  {
    Mabs=_Mabsn_STY_p_M[i];
    Mlow=Mag(y[i],_mlim_STY_p_M,*_cosmo_STY_p_M);
    Llow=pow(10.,-0.4*Mlow);

    /* debido a un underflow, tuvimos que poner este if
       El 0.25 es debido a que gsl_sfi_gamma_inc llama a gsl_sf_gamma_inc_CF
       para x > 0.25 y nos devolvía un underflow cuando se cumplía la segunda
       condición -> los fuentes de gsl están en:
       /net/gladiolo/scratch/dabreu/local/SOURCES/gsl-1.8/specfunc/exp.c
       /net/gladiolo/scratch/dabreu/local/SOURCES/gsl-1.8/specfunc/gamma_inc.c
    */
    if(Llow/Lstar > 0.25 && (lfamo.alfa*log(Llow/Lstar) - Llow/Lstar) <= GSL_LOG_DBL_MIN)
    {
      log_gamma_int=GSL_LOG_DBL_MIN;
    }
    else
    {
      log_gamma_int=log(gsl_sf_gamma_inc(1+lfamo.alfa,Llow/Lstar));
    }
    /* log(lfamo.phistar) + log_gamma_int -> integral de la función de Schecter
    entre Llow e inf */
    logL-= log(Schechter_M(Mabs,lfamo)) - log(lfamo.phistar) - log_gamma_int;
    
  }

  if(DEBUG) printf(" logL %g\n",logL);

  /* Aqui viene la parte de la poissoniana de npob */
  Ntot=Int_sch_M(lfamo,_zlow_STY_p_M,_zup_STY_p_M,_mlim_STY_p_M,*_cosmo_STY_p_M)*_strrad_STY_p_M/4./M_PI;
  logL-= (_ndata_STY_p_M*log(Ntot) - Ntot - gammln((double)_ndata_STY_p_M+1.));

  _iter_m_STY_p_M++;

  if(DEBUG) printf(" Ntot %g ndata %d rad %g zlw %f zup %f mag %f lfamophi %f\n",Ntot, _ndata_STY_p_M, _strrad_STY_p_M,_zlow_STY_p_M,_zup_STY_p_M,_mlim_STY_p_M,lfamo.phistar);

  if(DEBUG) printf(" Iter %d  logL %f par0 %g par1 %g par2 %g (%.10g)\n",_iter_m_STY_p_M,logL,p[0],p[1],exp(p[2]),p[2]);
  return(logL);
}

/* Errores utilizando derivadas numéricas con GSL */

void NumericalHessianCovars_STY_p_M(int n,double *magn,double *z,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf)
{
  gsl_matrix *gsl_hessian;
  gsl_vector *param;
  double **covar;
  double **bb;
  size_t i,j;

  double hessianStep;
  hessianStep = GSL_ROOT4_DBL_EPSILON;
  /* derivStep = 0.01; */

  gsl_hessian=gsl_matrix_alloc(3,3); 
  param=gsl_vector_alloc(3);
  gsl_vector_set(param, 0, par[0]);
  gsl_vector_set(param, 1, par[1]);
  gsl_vector_set(param, 2, par[2]);

  covar=matrix_d(3 ,3 );
  bb=matrix_d(3,1);

  struct Amoe_Funk_gsl_param_STY_p_M deriv_param;
  
  deriv_param.nData = _ndata_STY_p_M;
  deriv_param.magn = magn;
  deriv_param.z = z;

  gsl_multimin_function LogLFunction;
  LogLFunction.f = &Amoe_Funk_STY_p_M_main_gsl_multimin;
  LogLFunction.params = &deriv_param;
  LogLFunction.n = 3;
  
  printf("Vamos a por el gsl_hessian...\n");
  gsl_hessian_central(&LogLFunction, param, hessianStep, gsl_hessian);

  for(i=0;i<3;++i) 
  {
    for(j=0;j<3;++j) 
    {
      covar[i][j]=gsl_matrix_get(gsl_hessian,i,j);
    }
  }
  gaussj_d(covar,3,bb,1);

  lf->erralfa=sqrt(covar[0][0]);
  lf->errMstar=sqrt(covar[1][1]);
  lf->errphistar=sqrt(covar[2][2])*lf->phistar; /* Ya que p[2] esta en logaritmos */
  lf->covaralfaMstar=covar[0][1];
  lf->covaralfaphistar=covar[0][2]*lf->phistar;  /* Por la misma razon */
  lf->covarphistarMstar=covar[1][2]*lf->phistar;

  gsl_matrix_free(gsl_hessian);
  free_matrix_d(bb,3,1);
}

double Amoe_Funk_STY_p_M_main_gsl_multimin(const gsl_vector *x, void *params)
{
  int nData           = ((struct Amoe_Funk_gsl_param_STY_p_M *)params)->nData;
  double* magn        = ((struct Amoe_Funk_gsl_param_STY_p_M *)params)->magn;
  double* z           = ((struct Amoe_Funk_gsl_param_STY_p_M *)params)->z;
  double value;
  double* lf_param = vector_d(3);

  lf_param[0] = gsl_vector_get(x,0);
  lf_param[1] = gsl_vector_get(x,1);
  lf_param[2] = gsl_vector_get(x,2);

  value=Amoe_Funk_STY_p_M_main(nData, magn, z, lf_param);
  free(lf_param);
  return value;
}

