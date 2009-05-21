#include "modulos.h"
#include "gsl_hessian.h"
#include <gsl/gsl_machine.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_interp.h>
//#define FTOL  1e-10
#define FTOL  1e-9
#define FTOL2 1e-6
#define FTOL3 1e-7
#define MAXITER  1000
#define MAXITER2 120
#define MAXITERVEGAS 10
#define VERBOSE 0
#define DEBUG  1
#define DEBUG2 1
#define DEBUG3 1
#define DEBUG4 1
#define DEBUGPLOT 0
#define TOLERR 0.07 

#define ZMAX 10
#define NRHO 200

/* too much output for DEBUG4 (> 1Gb) */

/* #define TOLERR 0.0001 */
/* Estructura para contener los parámetros de Amo..main_gsl_multimin */
struct Amoe_Funk_gsl_param_STY_gmz_p_f_M_wC
{
  int nData;
  double *magn;
  double *z;
};

void   PrepareGlobalVars_STY_gmz_p_f_M_wC(double *z, double *magn);
double Amoe_Funk_STY_gmz_p_f_M_wC_main(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_gmz_p_f_M_wC_main_gsl_multimin(const gsl_vector *x, void *params);
double vegas_funk_numerator_STY_gmz_p_f_M_wC(double *x, size_t dim, void *params);
double vegas_funk_denominator_STY_gmz_p_f_M_wC(double *x, size_t dim, void *params);
double compute_rho_STY_gmz_p_f_M_wC();
double vegas_integrate_STY_gmz_p_f_M_wC(gsl_monte_function * f,
                                        double xl[], double xu[],
                                        size_t dim, size_t calls,
                                        double *error);
void   NumericalHessianCovars_STY_gmz_p_f_M_wC(int n,double *magn,double *errmagn, double *z, double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf);


struct Schlf_M *_lf_STY_gmz_p_f_M_wC;
double *_errmagn_STY_gmz_p_f_M_wC;
struct cosmo_param *_cosmo_STY_gmz_p_f_M_wC;
//double _mlim_STY_gmz_p_f_M_wC;
struct fermifsel_M _fsel_STY_gmz_p_f_M_wC;
double _z_i_STY_gmz_p_f_M_wC;
double _errz_i_STY_gmz_p_f_M_wC;
double *_errz_STY_gmz_p_f_M_wC;
double _magDistn_i_STY_gmz_p_f_M_wC;
double _errmagDistn_i_STY_gmz_p_f_M_wC;
double _magSeln_i_STY_gmz_p_f_M_wC;
double *_magSeln_STY_gmz_p_f_M_wC;
double *_magDistn_STY_gmz_p_f_M_wC;
double *_errmagDistn_STY_gmz_p_f_M_wC;
double _color_mean_STY_gmz_p_f_M_wC;
double _color_stddev_STY_gmz_p_f_M_wC;
double _color_i_STY_gmz_p_f_M_wC;

int _ndata_STY_gmz_p_f_M_wC;
int _iter_m_STY_gmz_p_f_M_wC;
int _iter_c_STY_gmz_p_f_M_wC;
double _MLmax_STY_gmz_p_f_M_wC;
double _zlow_STY_gmz_p_f_M_wC;
double _zup_STY_gmz_p_f_M_wC;
double _strrad_STY_gmz_p_f_M_wC;

double _zstep;
double _nRho=200.;
double _rho_STY_gmz_p_f_M_wC[NRHO];
double _zRho_STY_gmz_p_f_M_wC[NRHO];
double _dVdz_STY_gmz_p_f_M_wC[NRHO];
gsl_interp * _rho_interp_STY_gmz_p_f_M_wC;
gsl_interp_accel * _rho_interp_accel_STY_gmz_p_f_M_wC;
gsl_interp * _dVdz_interp_STY_gmz_p_f_M_wC;
gsl_interp_accel * _dVdz_interp_accel_STY_gmz_p_f_M_wC;
gsl_rng * _random_gen_STY_gmz_p_f_M_wC;

int  MLA_STY_gmz_p_f_M_wC(int n,double *magSeln, double *magDistn, double *errmagDistn, double color_mean, double color_stddev, double *z, double *errz, struct fermifsel_M fsel, double strrad, double zlow, double zup, struct cosmo_param cosmo, struct Schlf_M *lf, struct MLProcessInfo *mlinfo)
  {

  double par[3];
  double sigpar[3];
  int i;
  int iter_amo;
  const gsl_rng_type *T_rng;

  /* Variables to use for vvmax fitting */
  struct Steplf_M  lfvvmax;
  double minMabs, maxMabs;
  double *Mabs;
  double chisq;
  struct Schlf_M  lffit;

  Mabs=vector_d(n);

  /* Copias globales de estas variables para que se vean en toda la subrutina */
  
  _magSeln_STY_gmz_p_f_M_wC=magSeln;
  _magDistn_STY_gmz_p_f_M_wC=magDistn;
  _errmagDistn_STY_gmz_p_f_M_wC=errmagDistn;
  _color_mean_STY_gmz_p_f_M_wC=color_mean;
  _color_stddev_STY_gmz_p_f_M_wC=color_stddev;
  _errz_STY_gmz_p_f_M_wC=errz;
  _ndata_STY_gmz_p_f_M_wC=n;
  _cosmo_STY_gmz_p_f_M_wC=&cosmo;
  //_mlim_STY_gmz_p_f_M_wC=mlim;
  _fsel_STY_gmz_p_f_M_wC=fsel;
  _strrad_STY_gmz_p_f_M_wC=strrad;
  _zlow_STY_gmz_p_f_M_wC=zlow;
  _zup_STY_gmz_p_f_M_wC=zup;

  /* Init random number generator */
  gsl_rng_env_setup ();
  T_rng = gsl_rng_default;
  _random_gen_STY_gmz_p_f_M_wC = gsl_rng_alloc (T_rng);
  

  if(_zup_STY_gmz_p_f_M_wC==0) _zup_STY_gmz_p_f_M_wC=ZMAX;

  _iter_m_STY_gmz_p_f_M_wC=0;

  if(DEBUG3) {
    for(i=0;i<n;i++) printf(" Entrada x %g\n",magDistn[i]);
  }

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

  for(i=0;i<n;i++)   Mabs[i]=Mag(z[i],magDistn[i],cosmo);
  MinMax_d(n,Mabs,&minMabs,&maxMabs);
  for(i=0;i<=lfvvmax.nbin;i++) lfvvmax.magni[i]=minMabs+i*(maxMabs-minMabs)/lfvvmax.nbin;
  VVmax_M(n,magSeln,magDistn,z,fsel.magcut,strrad,zlow,zup,cosmo,&lfvvmax);
/*   cpgopen("?"); */
/*   PlotStepLF_M(lfvvmax); */
  if(DEBUG) printf(" Salida VVmax\n");
  if(DEBUG) for(i=0;i<lfvvmax.nbin;i++) printf(" Mabs %g - %g LF %g\n",lfvvmax.magni[i],lfvvmax.magni[i+1],lfvvmax.lnlf[i]/log(10));
  FitSch2StepLF_M(lfvvmax,&lffit, &chisq);
  if(DEBUG) {
    printf(" Después ajuste MRQ\n");
    printf(" Schechter FIT\n");
    printf(" Mstar (:  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",lffit.Mstar,lffit.alfa,lffit.phistar,log10(lffit.phistar));
    printf(" E_Mstar:    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lffit.errMstar,lffit.erralfa,lffit.errphistar,lffit.errphistar/lffit.phistar/log(10.));
    printf(" Covar(Lstar,alpha): %g Covar(Mstar,alpha): %g\n",lffit.covaralfaMstar,lffit.covaralfaMstar);
    printf(" Covar(Mstar,Phistar): %g Covar(Mstar,log(Phistar)): %g\n",lffit.covarphistarMstar,lffit.covarphistarMstar/lffit.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lffit.covaralfaphistar,lffit.covaralfaphistar/lffit.phistar/log(10.));
  }
  free(Mabs);
  free(lfvvmax.magni);
  free(lfvvmax.errmagni);
  free(lfvvmax.lnlf);
  free(lfvvmax.errlnlf);
  free(lfvvmax.lf);
  free(lfvvmax.errlf);
  free_matrix_d(lfvvmax.covarlnlf,lfvvmax.nbin,lfvvmax.nbin);

/*   PlotStepSchLF_M(lfvvmax,lffit); */

/*   cpgclos(); */

  /* dabreu */
  /* struct MLProcessInfo mlinfo2;
  i=MLA_STY_p_M(n,magn,z,mlim,strrad,zlow,zup,cosmo,&lffit, &mlinfo2);
  printf(" STY as a first solution.\n");
  printf(" STY -> Mstar: %g alpha: %g phistar: %g\n",lffit.Mstar,lffit.alfa,lffit.phistar); */

  printf(" Computing LF (method STY_gmz_p_f_m_wC)...\n");

  iter_amo=0;
  while(iter_amo==0)
  {
    par[0]=lffit.alfa;
    par[1]=lffit.Mstar;
    par[2]=log(lffit.phistar);
   /* to feed the code with the initial solution uncomment following lines */
    /* par[0]=-1.3;
    par[1]=-20.4;
    par[2]=log(0.0033); */

    sigpar[0]=10.*lffit.erralfa;
    sigpar[1]=10.*lffit.errMstar;
    sigpar[2]=10.*lffit.errphistar/lffit.phistar;
    iter_amo=Amoeba_d(n,magDistn,z,3,par,sigpar,FTOL,MAXITER,Amoe_Funk_STY_gmz_p_f_M_wC_main);
    if(DEBUG) printf(" iteramo %d\n",iter_amo);
    lf->alfa=par[0];
    lf->Mstar=par[1];
    lf->phistar=exp(par[2]);
    if(DEBUG)
    {
      printf(" Solucion MALA\n");
      printf(" Mstar :  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",lf->Mstar,lf->alfa,lf->phistar,log10(lf->phistar));
      printf(" E_log(Lstar):    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lf->errMstar,lf->erralfa,lf->errphistar,lf->errphistar/lf->phistar/log(10.));
    } 
  }
  /* Info that will be output in mlinfo */
  _MLmax_STY_gmz_p_f_M_wC=Amoe_Funk_STY_gmz_p_f_M_wC_main(n,magDistn,z,par);
  mlinfo->nIter = iter_amo;
  mlinfo->MLmax = _MLmax_STY_gmz_p_f_M_wC;

  /* Meto la solucion en la salida */
 
  lf->alfa=par[0];
  lf->Mstar=par[1];
  lf->phistar=exp(par[2]);

  /* Estimacion de los errores en los parametros */

  printf(" Computing errors in LF parameters...\n");

  NumericalHessianCovars_STY_gmz_p_f_M_wC(n,magDistn,errmagDistn,z,par,sigpar,_fsel_STY_gmz_p_f_M_wC.magcut,cosmo,lf);
  if(DEBUG) printf(" Solucion final: Mstar %.15g +/- %.15g alpha %.4g +/- %.4g\n",lf->Mstar,lf->errMstar,lf->alfa,lf->erralfa);

  if(DEBUG) printf(" MLF %g\n",_MLmax_STY_gmz_p_f_M_wC);

  /* Code exit status */
  gsl_rng_free(_random_gen_STY_gmz_p_f_M_wC);
  if(iter_amo>=MAXITER-1) return(2);
  return(0);
}

double Amoe_Funk_STY_gmz_p_f_M_wC_main(int n, double *x, double *y, double *p) {

  int i; 
  double logL=0.;
  double logLold=0.;
  struct Schlf_M lfamo;
  double Lstar;
  double Mabs;
  double Mlow;
  double xmin;
/*   double Mup=-30; */
/*   double Lup; */
  double Ntot;
  double x1,x2;
  /* int npa=21;
  int npb=21;
  int npc=31; */
  double probarriba, errprobarriba;
  double probabajo, errprobabajo;

  double colori;
  double logColor=0.0;

  double intsup;

  size_t dim_den=4, dim_num=2;
  double xl_den[4], xu_den[4];
  double xl_num[2], xu_num[2];


  gsl_monte_function G_den = { &vegas_funk_denominator_STY_gmz_p_f_M_wC, dim_den, 0 };
  gsl_monte_function G_num = { &vegas_funk_numerator_STY_gmz_p_f_M_wC, dim_num, 0 };

//  too slow (test)
//  size_t calls_den = 5000000;
//  size_t calls_num = 50000;
//  regular
  size_t calls_den = 500000;
  size_t calls_num = 5000;
//  fast (test)
//  size_t calls_den = 50000;
//  size_t calls_num = 500;
//  very fast (test)
//  size_t calls_den = 500;
//  size_t calls_num = 500;


  lfamo.alfa=p[0];
  lfamo.Mstar=p[1];
  lfamo.phistar=exp(p[2]);

  _lf_STY_gmz_p_f_M_wC=&lfamo;

  logL=0.;
  logLold=0.;

  Lstar=pow(10.,-0.4*lfamo.Mstar);
  
  if(DEBUG3) printf(" Entra con  par0 %g par1 %g par2 %g\n",p[0],p[1],p[2]);

  intsup=incom(1+lfamo.alfa,100.);
  for(i=0;i<_ndata_STY_gmz_p_f_M_wC;i++) 
  {
    /* Contribucion del color en el numerador */
    colori=_magDistn_STY_gmz_p_f_M_wC[i] - _magSeln_STY_gmz_p_f_M_wC[i];
    _color_i_STY_gmz_p_f_M_wC = colori;
    if (_color_stddev_STY_gmz_p_f_M_wC==0)
    {
      logColor=0.0;
    }
    else
    {
      logColor=lngaussian(colori, _color_mean_STY_gmz_p_f_M_wC, _color_stddev_STY_gmz_p_f_M_wC);
    }

    _z_i_STY_gmz_p_f_M_wC=y[i];
    _errz_i_STY_gmz_p_f_M_wC=_errz_STY_gmz_p_f_M_wC[i];
    _magDistn_i_STY_gmz_p_f_M_wC=x[i];
    _errmagDistn_i_STY_gmz_p_f_M_wC=_errmagDistn_STY_gmz_p_f_M_wC[i];
    if(_errmagDistn_i_STY_gmz_p_f_M_wC==0)
    {
      /* _errmagn_i_STY_gm_p_M=GSL_SQRT_DBL_MIN; */
      _errmagDistn_i_STY_gmz_p_f_M_wC=0.0001;
      /* printf("Using GSL_SQRT_DBL_MIN instead of 0 for errmag.\n"); */
    }
    if(_errz_i_STY_gmz_p_f_M_wC==0) _errz_i_STY_gmz_p_f_M_wC=0.00001;

    Mabs=Mag(y[i],x[i],*_cosmo_STY_gmz_p_f_M_wC);
    Mlow=Mag(y[i],_fsel_STY_gmz_p_f_M_wC.magcut+colori,*_cosmo_STY_gmz_p_f_M_wC);
    xmin =pow(10.,-0.4*Mlow)/Lstar;
    if(xmin>=100) 
    {
      logL+=GSL_LOG_DBL_MAX;
      if(DEBUG3) printf("Pasando por logL+=GSL_LOG_DBL_MAX\n");
      if(DEBUG3) printf("y[i] %g\n",y[i]);
    }
    else 
    {
      x1=lfamo.Mstar-5; /* de dónde sale este 5??? */
      x2=Mag(y[i],_fsel_STY_gmz_p_f_M_wC.magcut+6*_errmagDistn_i_STY_gmz_p_f_M_wC,*_cosmo_STY_gmz_p_f_M_wC);
      if(DEBUG4) printf("x1 %g x2 %g y[i] %g _errmagDistn_i %g\n",x1,x2,y[i],_errmagDistn_i_STY_gmz_p_f_M_wC);

      /* DENOMINATOR */
      if(DEBUG3) printf(" Computing denominator.\n");
      /* compute limits */
      /* xl_den={_zlow_STY_gmz_p_f_M_wC,-6*_errz_i_STY_gmz_p_f_M_wC,x1,-6*_errmagDistn_i_STY_gmz_p_f_M_wC}; 
      xu_den={_zup_STY_gmz_p_f_M_wC,6*_errz_i_STY_gmz_p_f_M_wC,x2,6*_errmagDistn_i_STY_gmz_p_f_M_wC}; */
      xl_den[0]=_zlow_STY_gmz_p_f_M_wC;
      xl_den[1]=-6*_errz_i_STY_gmz_p_f_M_wC;
      xl_den[2]=x1;
      xl_den[3]=-6*_errmagDistn_i_STY_gmz_p_f_M_wC; 
      xu_den[0]=_zup_STY_gmz_p_f_M_wC;
      xu_den[1]=6*_errz_i_STY_gmz_p_f_M_wC;
      xu_den[2]=x2;
      xu_den[3]=6*_errmagDistn_i_STY_gmz_p_f_M_wC;
      if(DEBUG3) printf(" Limits xl: %g %g %g %g\n",xl_den[0],xl_den[1],xl_den[2],xl_den[3]);
      if(DEBUG3) printf(" Limits xu: %g %g %g %g\n",xu_den[0],xu_den[1],xu_den[2],xu_den[3]);
      /* integration */
      probabajo = vegas_integrate_STY_gmz_p_f_M_wC
            (&G_den, xl_den, xu_den, dim_den, calls_den,
             &errprobabajo);

      if(DEBUG2) printf(" Abajo %g +- %g\n",probabajo,errprobabajo);

      /* We have to compute them again because Funk? may have overriden them */
      _z_i_STY_gmz_p_f_M_wC=y[i];
      _errz_i_STY_gmz_p_f_M_wC=_errz_STY_gmz_p_f_M_wC[i];//+GSL_SQRT_DBL_MIN;
      _magDistn_i_STY_gmz_p_f_M_wC=x[i];
      _errmagDistn_i_STY_gmz_p_f_M_wC=_errmagDistn_STY_gmz_p_f_M_wC[i];//+GSL_SQRT_DBL_MIN;

      if(_errmagDistn_i_STY_gmz_p_f_M_wC==0) _errmagDistn_i_STY_gmz_p_f_M_wC=GSL_DBL_EPSILON*1000;
      if(_errz_i_STY_gmz_p_f_M_wC==0) _errz_i_STY_gmz_p_f_M_wC=GSL_DBL_EPSILON*1000;

      /* NUMERATOR*/
      if(DEBUG4) printf(" errors: %g %g\n",_errz_i_STY_gmz_p_f_M_wC,_errmagDistn_i_STY_gmz_p_f_M_wC);
      if(DEBUG3) printf(" Computing numerator.\n");
      /* x[0] are zreal integral limits */
      xl_num[0]=_z_i_STY_gmz_p_f_M_wC-6*_errz_i_STY_gmz_p_f_M_wC;
      xu_num[0]=_z_i_STY_gmz_p_f_M_wC+6*_errz_i_STY_gmz_p_f_M_wC;
      /* x[1] are mreal-mobs integral limits */
      if(4 * _fsel_STY_gmz_p_f_M_wC.deltamag >  _errmagDistn_i_STY_gmz_p_f_M_wC)
      {
             /* The Fermi factor has smooth gradients along the gaussian (the factor 4 has been tested) */
             xl_num[1] = -6*_errmagDistn_i_STY_gmz_p_f_M_wC; 
             xu_num[1] = +6*_errmagDistn_i_STY_gmz_p_f_M_wC; 
             if(DEBUG3) printf(" Caso 1Limits xl_num: %g %g\n",xl_num[0],xl_num[1]);
             if(DEBUG3) printf(" Caso 1Limits xu_num: %g %g\n",xu_num[0],xu_num[1]);
             /* integration */
             probarriba = vegas_integrate_STY_gmz_p_f_M_wC
                 (&G_num, xl_num, xu_num, dim_num, calls_num,
                  &errprobarriba);
      }
      else
      {
             /* The Fermi factor changes abrouptly within the gaussian */
             double magcut_sigma =
                   (_fsel_STY_gmz_p_f_M_wC.magcut -
                    _magDistn_i_STY_gmz_p_f_M_wC) /
		     _errmagDistn_i_STY_gmz_p_f_M_wC;
             if(DEBUG3) printf(" magcut_sigma %g\n", magcut_sigma);
             if(magcut_sigma < -6)
             {
                 /* The integral will be almost 0. Since the gaussian
                    is ~0 when Fermi = 1 and Fermi~0 when gaussian 
                    is relevant. */
                 /* The -10 * deltamag factor is to allow a little bit
                    of integration into the Fermi~1 regime */
                 xl_num[1] = -6*_errmagDistn_i_STY_gmz_p_f_M_wC -
                              10 * _fsel_STY_gmz_p_f_M_wC.deltamag;; 
                 xu_num[1] = +6*_errmagDistn_i_STY_gmz_p_f_M_wC; 
                 if(DEBUG3) printf(" Caso 2Limits xl_num: %g %g\n",
                                   xl_num[0],xl_num[1]);
                 if(DEBUG3) printf(" Caso 2Limits xu_num: %g %g\n",
                                   xu_num[0],xu_num[1]);
                 /* integration */
                 probarriba = vegas_integrate_STY_gmz_p_f_M_wC
                     (&G_num, xl_num, xu_num, dim_num, calls_num/10,
                      &errprobarriba);

             }
             else if(magcut_sigma > 6)
             {
                 /* The gaussian is all in the Fermi ~=1 regime */
                 xl_num[1] = -6*_errmagDistn_i_STY_gmz_p_f_M_wC; 
                 xu_num[1] = +6*_errmagDistn_i_STY_gmz_p_f_M_wC; 
                 if(DEBUG3) printf(" Caso 3Limits xl_num: %g %g\n",
                                   xl_num[0],xl_num[1]);
                 if(DEBUG3) printf(" Caso 3Limits xu_num: %g %g\n",
                                   xu_num[0],xu_num[1]);
                 /* integration */
                 probarriba = vegas_integrate_STY_gmz_p_f_M_wC
                     (&G_num, xl_num, xu_num, dim_num, calls_num,
                      &errprobarriba);
             }
             else
             {
                 /* The Fermi cut is in the middle of the gaussian.
                    There will be two integration intervals */
                 probarriba = 0;
                 /* First interval */
                 /* The -6 *deltmag factor ensures that xl[1] < xu[1] */
                 xl_num[1] = -6 * _errmagDistn_i_STY_gmz_p_f_M_wC
                             -10 * _fsel_STY_gmz_p_f_M_wC.deltamag; 
                 xu_num[1] = magcut_sigma *
                     _errmagDistn_i_STY_gmz_p_f_M_wC -
                     10 * _fsel_STY_gmz_p_f_M_wC.deltamag; 
                 if(DEBUG3) printf(" Caso 4.1Limits xl_num: %g %g\n",
                                   xl_num[0],xl_num[1]);
                 if(DEBUG3) printf(" Caso 4.1Limits xu_num: %g %g\n",
                                   xu_num[0],xu_num[1]);
                 /* integration of first half */
                 probarriba += vegas_integrate_STY_gmz_p_f_M_wC
                     (&G_num, xl_num, xu_num, dim_num, calls_num,
                      &errprobarriba);
                 /* Second interval */
                 if(_fsel_STY_gmz_p_f_M_wC.deltamag != 0)
                 {
                     xl_num[1] = magcut_sigma *
                          _errmagDistn_i_STY_gmz_p_f_M_wC -
                   	  10 * _fsel_STY_gmz_p_f_M_wC.deltamag; 
                     xu_num[1] = magcut_sigma * 
                          _errmagDistn_i_STY_gmz_p_f_M_wC +
                          10 * _fsel_STY_gmz_p_f_M_wC.deltamag; 
                     if(DEBUG3) printf("Caso 4.2Limits xl_num: %g %g\n",
                                       xl_num[0],xl_num[1]);
                     if(DEBUG3) printf("Caso 4.2Limits xu_num: %g %g\n",
                                       xu_num[0],xu_num[1]);
                     /* integration of second half */
                     probarriba += vegas_integrate_STY_gmz_p_f_M_wC
                         (&G_num, xl_num, xu_num, dim_num, calls_num,
                          &errprobarriba);
           }
        }
      }
      if(DEBUG2) printf(" Calculo arriba %g +- %g magn %g err %g magnl %g\n",probarriba,errprobarriba,_magDistn_i_STY_gmz_p_f_M_wC,_errmagDistn_i_STY_gmz_p_f_M_wC,_fsel_STY_gmz_p_f_M_wC.magcut);

      if(probarriba==0 || probabajo==0) logL+=GSL_LOG_DBL_MAX; 
      else logL-= log(probarriba) - log(probabajo);   /* Perfectamente testado */
      logL-= logColor;

      if(DEBUG3) printf(" iobj %d logL %f loglold %f      sch %g    pa %g pb %g (%g)  xmin %g x1 %g x2 %g magn %g err %g\n",i,logL,logLold,Schechter_M(Mabs,lfamo),log(probarriba),log(probabajo),probabajo,xmin,x1,x2,_magDistn_i_STY_gmz_p_f_M_wC,_errmagDistn_i_STY_gmz_p_f_M_wC);
      if(DEBUG3) printf(" logColor %g colori %g\n",logColor,colori);
      if(DEBUG3) printf(" #LPA %g %g\n",log(probarriba),log(probabajo));
    }
  } /* end of objects loop */

  /* free gsl interp */
  //gsl_interp_free(_rho_interp_STY_gmz_p_f_M_wC);
  //gsl_interp_accel_free(_rho_interp_accel_STY_gmz_p_f_M_wC);
  //gsl_interp_free(_dVdz_interp_STY_gmz_p_f_M_wC);
  //gsl_interp_accel_free(_dVdz_interp_accel_STY_gmz_p_f_M_wC);

  /* Aqui viene la parte de la poissoniana de npob */
  /* TODO: Change to Int_sch_f_M_wC (when done) */
  /* TODO: Change integration method in Int_sch_f_M_wC */
  Ntot=Int_sch_M_wC(lfamo,_zlow_STY_gmz_p_f_M_wC,_zup_STY_gmz_p_f_M_wC,_color_mean_STY_gmz_p_f_M_wC,_color_stddev_STY_gmz_p_f_M_wC,_fsel_STY_gmz_p_f_M_wC.magcut,*_cosmo_STY_gmz_p_f_M_wC)*_strrad_STY_gmz_p_f_M_wC/4./M_PI; 
  logL-=    (_ndata_STY_gmz_p_f_M_wC*log(Ntot) - Ntot - gammln((double)_ndata_STY_gmz_p_f_M_wC+1.)); 
  
  if(DEBUG2) printf(" NTOT %f ndata*log(Ntot) %f gamm %f\n",Ntot,_ndata_STY_gmz_p_f_M_wC*log(Ntot),gammln((double)_ndata_STY_gmz_p_f_M_wC)+1.);
  if(DEBUG3) printf(" NTOT: color_mean %g color_stddev %g\n",_color_mean_STY_gmz_p_f_M_wC,_color_stddev_STY_gmz_p_f_M_wC);

  _iter_m_STY_gmz_p_f_M_wC++;
  if(DEBUG) printf(" Iter %d  logL %g logLold %g par0 %g par1 %g par2 %g\n",_iter_m_STY_gmz_p_f_M_wC,logL,logLold,p[0],p[1],p[2]);
  return(logL); 
}

double vegas_funk_numerator_STY_gmz_p_f_M_wC (double *x, size_t dim, void *params)
{
  double res;
  double zreal;
  double zobs;
  double mDistReal;
  double mDistObs;
  double Mabs;
  double logfacerrz;
  double logfacerrm;
  double logfacLF;
  double facsel;
  double dVdzreal, rhoz;

  if(DEBUG4) printf(" inside vegas_funk_numerator.\n");
  /* x[0] -> zreal
     x[1] -> mreal - mobs */

  zreal = x[0];
  zobs  = _z_i_STY_gmz_p_f_M_wC;
  mDistObs  = _magDistn_i_STY_gmz_p_f_M_wC;
  mDistReal = x[1] + mDistObs;
  Mabs = Mag(zreal,mDistReal,*_cosmo_STY_gmz_p_f_M_wC);

  if ( zreal < _zlow_STY_gmz_p_f_M_wC || zreal > _zup_STY_gmz_p_f_M_wC)
  {
    if(DEBUG4) printf(" out petando zreal: zreal<0\n");
    return(0);
  } 

  /* if (mreal > _fsel_STY_gmz_p_f_M_wC.magcut+3*_fsel_STY_gmz_p_f_M_wC.deltamag)
  {
    if(DEBUG4) printf(" out petando mreal: mreal %g magcut %g\n",mreal,_fsel_STY_gmz_p_f_M_wC.magcut);
    return(0);
  } Pertenece a los límites */

  logfacLF = log(Schechter_M(Mabs,*_lf_STY_gmz_p_f_M_wC));
  logfacerrz = lngaussian(zobs, zreal, 
                          _errz_i_STY_gmz_p_f_M_wC);
  logfacerrm = lngaussian(mDistObs, mDistReal, 
                          _errmagDistn_i_STY_gmz_p_f_M_wC);
//  dVdz=Lagr2_d(_zRho_STY_gmz_p_f_M_wC,_dVdz_STY_gmz_p_f_M_wC,NRHO,zreal);
//  rhoz=Lagr2_d(_zRho_STY_gmz_p_f_M_wC,_rho_STY_gmz_p_f_M_wC,NRHO,zreal);
  //dVdz=gsl_interp_eval(_dVdz_interp_STY_gmz_p_f_M_wC,_zRho_STY_gmz_p_f_M_wC,_dVdz_STY_gmz_p_f_M_wC,zreal,_dVdz_interp_accel_STY_gmz_p_f_M_wC);
  dVdzreal = dVdz(zreal,*_cosmo_STY_gmz_p_f_M_wC);
  //rhoz=gsl_interp_eval(_rho_interp_STY_gmz_p_f_M_wC,_zRho_STY_gmz_p_f_M_wC,_rho_STY_gmz_p_f_M_wC,zreal,_rho_interp_accel_STY_gmz_p_f_M_wC);
  /* facsel = Fermi(mreal,_fsel_STY_gmz_p_f_M_wC.magcut+_color_i_STY_gmz_p_f_M_wC,_fsel_STY_gmz_p_f_M_wC.deltamag); */
  facsel = Fermi(mDistReal-_color_i_STY_gmz_p_f_M_wC,
                 _fsel_STY_gmz_p_f_M_wC.magcut,
                 _fsel_STY_gmz_p_f_M_wC.deltamag);
  res=exp(logfacerrm+logfacerrz+logfacLF);
  /* test without rhoz */
  rhoz=1.0;
  res=res*dVdzreal*rhoz*facsel;

  if(DEBUG4) printf("Num todopati: zreal %10g zobs %10g mreal %10g mobs %10g Mabs %10g\n",zreal,zobs,mDistReal,mDistObs,Mabs);
  if(DEBUG4) printf("Num morralla: logfacLF %g logfacerrz %g logfacerrm %g dVdzreal %g rhoz %g  facsel %g\n",logfacLF,logfacerrz,logfacerrm,dVdzreal,rhoz,facsel);
  if(DEBUG4) printf("Num result: %g\n",res);
  return(res);
}

double vegas_funk_denominator_STY_gmz_p_f_M_wC (double *x, size_t dim, void *params)
{
  double res;
  double zreal;
  double zobs;
  double mreal;
  double mobs;
  double Mabs;
  double logfacerrz;
  double logfacerrm;
  double logfacLF;
  double facsel;
  double dVdzreal, rhoz;

  if(DEBUG4) printf(" inside vegas_funk_denominator.\n");

  zobs = x[0];
  zreal = zobs + x[1];
  mobs = mag(zreal, x[2],*_cosmo_STY_gmz_p_f_M_wC);
  mreal = mobs + x[3];
  Mabs = Mag(zreal,mreal,*_cosmo_STY_gmz_p_f_M_wC);

  if ( zreal < 0)
  {
    if(DEBUG4) printf(" out petando zreal: zreal<0\n");
    return(0);
  } 

  if (mreal > _fsel_STY_gmz_p_f_M_wC.magcut+3*_fsel_STY_gmz_p_f_M_wC.deltamag)
  {
    if(DEBUG4) printf(" out petando mreal: mreal %g magcut %g\n",mreal,_fsel_STY_gmz_p_f_M_wC.magcut);
    return(0);
  } 

  logfacLF = log(Schechter_M(Mabs,*_lf_STY_gmz_p_f_M_wC));
  logfacerrz = lngaussian(zobs, zreal, _errz_i_STY_gmz_p_f_M_wC);
  logfacerrm = lngaussian(mobs, mreal, _errmagDistn_i_STY_gmz_p_f_M_wC);
  //dVdz=Lagr2_d(_zRho_STY_gmz_p_f_M_wC,_dVdz_STY_gmz_p_f_M_wC,NRHO,zreal);
  //rhoz=Lagr2_d(_zRho_STY_gmz_p_f_M_wC,_rho_STY_gmz_p_f_M_wC,NRHO,zreal);
  //dVdz=gsl_interp_eval(_dVdz_interp_STY_gmz_p_f_M_wC,_zRho_STY_gmz_p_f_M_wC,_dVdz_STY_gmz_p_f_M_wC,zreal,_dVdz_interp_accel_STY_gmz_p_f_M_wC);
  dVdzreal = dVdz(zreal,*_cosmo_STY_gmz_p_f_M_wC);
  //rhoz=gsl_interp_eval(_rho_interp_STY_gmz_p_f_M_wC,_zRho_STY_gmz_p_f_M_wC,_rho_STY_gmz_p_f_M_wC,zreal,_rho_interp_accel_STY_gmz_p_f_M_wC);
  facsel = Fermi(mreal,_fsel_STY_gmz_p_f_M_wC.magcut,_fsel_STY_gmz_p_f_M_wC.deltamag); /* TODO: seguro que no va colori? */
  res=exp(logfacerrm+logfacerrz+logfacLF);
  /* test without rhoz */
  rhoz=1.0;
  res=res*dVdzreal*rhoz*facsel;
  if(DEBUG4) printf(" lasx: x[0] %g x[1] %g x[2] %g x[3] %g\n",x[0],x[1],x[2],x[3]);
  if(DEBUG4) printf("Den todopati: zreal %10g zobs %10g mreal %10g mobs %10g Mabs %10g\n",zreal,zobs,mreal,mobs,Mabs);
  if(DEBUG4) printf("Den morralla: logfacLF %g logfacerrz %g logfacerrm %g dVdzreal %g rhoz %g  facsel %g\n",logfacLF,logfacerrz,logfacerrm,dVdzreal,rhoz,facsel);
  if(DEBUG4) printf("Den result: %g\n",res);
  return(res); 
}

/* Errores utilizando derivadas numéricas con GSL */

void NumericalHessianCovars_STY_gmz_p_f_M_wC(int n,double *magn,double *errmagn,double *z,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf)
{
  gsl_matrix *gsl_hessian;
  gsl_vector *param;
  double **covar;
  double **bb;
  size_t i,j;

  double hessianStep = GSL_ROOT4_DBL_EPSILON;
  /* derivStep = 0.01; */

  gsl_hessian=gsl_matrix_alloc(3,3);
  param=gsl_vector_alloc(3);
  gsl_vector_set(param, 0, par[0]);
  gsl_vector_set(param, 1, par[1]);
  gsl_vector_set(param, 2, par[2]);

  covar=matrix_d(3 ,3 );
  bb=matrix_d(3,1);

  struct Amoe_Funk_gsl_param_STY_gmz_p_f_M_wC deriv_param;

  deriv_param.nData = _ndata_STY_gmz_p_f_M_wC;
  deriv_param.magn = magn;
  deriv_param.z = z;

  gsl_multimin_function LogLFunction;
  LogLFunction.f = &Amoe_Funk_STY_gmz_p_f_M_wC_main_gsl_multimin;
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

double Amoe_Funk_STY_gmz_p_f_M_wC_main_gsl_multimin(const gsl_vector *x, void *params)
{
  int nData           = ((struct Amoe_Funk_gsl_param_STY_gmz_p_f_M_wC *)params)->nData;
  double* magn        = ((struct Amoe_Funk_gsl_param_STY_gmz_p_f_M_wC *)params)->magn;
  double* z           = ((struct Amoe_Funk_gsl_param_STY_gmz_p_f_M_wC *)params)->z;
  double value;
  double* lf_param = vector_d(3);

  lf_param[0] = gsl_vector_get(x,0);
  lf_param[1] = gsl_vector_get(x,1);
  lf_param[2] = gsl_vector_get(x,2);


  value=Amoe_Funk_STY_gmz_p_f_M_wC_main(nData, magn, z, lf_param);
  free(lf_param);
  return value;
}

double vegas_integrate_STY_gmz_p_f_M_wC(gsl_monte_function * f,
                                        double xl[], double xu[],
                                        size_t dim, size_t calls,
                                        double *error)

{
  int itervegas = 0;
  double result, abserr;
  gsl_monte_vegas_state *state = gsl_monte_vegas_alloc (dim);
  /* Warm-up */
  gsl_monte_vegas_integrate (f, xl, xu, dim, calls/MAXITERVEGAS/5, 
                             _random_gen_STY_gmz_p_f_M_wC, 
                             state, &result, &abserr);
  if(DEBUG2) printf (" first iter: %g +- %g",result,abserr);
  if(DEBUG2) printf (" converging...\n");

  do
  {
    itervegas+=1;
    gsl_monte_vegas_integrate (f, xl, xu, dim, calls/MAXITERVEGAS, 
                               _random_gen_STY_gmz_p_f_M_wC,
                               state, &result, &abserr);
    if(DEBUG3) printf (" result = % .9f sigma = % .9f "
                       "chisq/dof = %.1f\n", result, abserr, state->chisq);
    if(DEBUG3) printf(" itervegas: %d\n",itervegas);
  }
  while (fabs (state->chisq - 1.0) > 0.5 && itervegas < MAXITERVEGAS);

  gsl_monte_vegas_free (state);
   
  return result;
}

