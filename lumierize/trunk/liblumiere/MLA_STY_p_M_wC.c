#include "modulos.h"
#include "gsl_hessian.h"
#include <gsl/gsl_machine.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

#define ZMIN 0.00001
#define FTOL  1e-12
#define FTOL2 1e-6
#define FTOL3 1e-7
#define MAXITER  1000
#define VERBOSE 0
#define DEBUG  0
#define DEBUG2 0
#define DEBUG3 0
#define DEBUGPLOT 0
#define MAXITERVEGAS 10

/* Estructura para contener los parámetros */
struct Amoe_Funk_gsl_param_STY_p_M_wC
{
  int nData;
  double *magn;
  double *z;
};

void prepareGlobalVars_STY_p_M_wC(double *z, double *magn);
double Amoe_Funk_STY_p_M_wC_main(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_p_M_wC_main_gsl_multimin(const gsl_vector *x, void *params);
void   NumericalHessianCovars_STY_p_M_wC(int n,double *magn,double *z,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf);
double vegas_funk_denominator_STY_p_M_wC (double *x, size_t dim, void *params);
double vegas_integrate_STY_p_M_wC(gsl_monte_function * f,
                                  double xl[], double xu[],
                                  size_t dim, size_t calls,
                                  double *error);

struct cosmo_param *_cosmo_STY_p_M_wC;
double _mlim_STY_p_M_wC;

int _ndata_STY_p_M_wC;
int _iter_m_STY_p_M_wC;
int _iter_c_STY_p_M_wC;
double *_pp_STY_p_M_wC;
double _MLmax_STY_p_M_wC;
double _xf_STY_p_M,_Tf_STY_p_M_wC;
double _sigi_STY_p_M_wC;
double _xtmp_STY_p_M_wC;
double _zlow_STY_p_M_wC,_zup_STY_p_M_wC,_strrad_STY_p_M_wC;
double *_magSeln_STY_p_M_wC;
double *_magDistn_STY_p_M_wC;
double _color_mean_STY_p_M_wC;
double _color_stddev_STY_p_M_wC;
double *_Mabsn_STY_p_M_wC;
double _z_i_STY_p_M_wC;
double _magSel_i_STY_p_M_wC;
struct Schlf_M *_lf_STY_p_M_wC;
gsl_rng * _random_gen_STY_gmz_p_f_M_wC;

int MLA_STY_p_M_wC(int n,double *magSeln, double *magDistn, double color_mean, double color_stddev, double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf, struct MLProcessInfo *mlinfo)
{

  double *y;
  double par[3];
  double sigpar[3];
  int i;
  int iter_amo;
  const gsl_rng_type *T_rng;

  /* Variables to use for vvmax fitting */
  struct Steplf_M  lfvvmax;
  double minMabs,maxMabs;
  struct Schlf_M  lffit;
  double chisq;

  /* Copias globales de estas variables para que se vean en toda la subrutina */
  _ndata_STY_p_M_wC=n;
  _cosmo_STY_p_M_wC=&cosmo;
  _mlim_STY_p_M_wC=mlim;
  _strrad_STY_p_M_wC=strrad;
  _zlow_STY_p_M_wC=zlow;
  _zup_STY_p_M_wC=zup;
  _magSeln_STY_p_M_wC=magSeln;
  _magDistn_STY_p_M_wC=magDistn;
  _color_mean_STY_p_M_wC=color_mean;
  _color_stddev_STY_p_M_wC=color_stddev;
  if(color_stddev==0)
  {
    _color_stddev_STY_p_M_wC=GSL_DBL_EPSILON*1000;
    printf("Using 1000*GSL_DBL_EPSILON instead of 0 for color_stddev.\n");
  }

  y=vector_d(n);
  _iter_m_STY_p_M_wC=0;

  /* Init random number generator */
  gsl_rng_env_setup ();
  T_rng = gsl_rng_default;
  _random_gen_STY_gmz_p_f_M_wC = gsl_rng_alloc (T_rng);

  iter_amo=MAXITER+1;

  /* Los límites en z se reajustan */
  _zlow_STY_p_M_wC = (_zlow_STY_p_M_wC < ZMIN ? ZMIN : _zlow_STY_p_M_wC);

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
  lfvvmax.ngalbin   =vector_i(lfvvmax.nbin);
  lfvvmax.covarlnlf =matrix_d(lfvvmax.nbin,lfvvmax.nbin);

  prepareGlobalVars_STY_p_M_wC(z,magDistn); /* inicializa _Mabsn_STY_p_M_wC */

  MinMax_d(n,_Mabsn_STY_p_M_wC,&minMabs,&maxMabs);
  for(i=0;i<=lfvvmax.nbin;i++) lfvvmax.magni[i]=minMabs+i*(maxMabs-minMabs)/lfvvmax.nbin;
  VVmax_M(n,magSeln,magDistn,z,mlim,strrad,zlow,zup,cosmo,&lfvvmax);

  /* Fit of the VVmax solution to a Schechter function */
  FitSch2StepLF_M(lfvvmax,&lffit, &chisq);

  /* Free vvmax related things */
  free(lfvvmax.magni);
  free(lfvvmax.errmagni);
  free(lfvvmax.lnlf);
  free(lfvvmax.errlnlf);
  free(lfvvmax.lf);
  free(lfvvmax.errlf);
  free(lfvvmax.ngalbin);
  free_matrix_d(lfvvmax.covarlnlf,lfvvmax.nbin,lfvvmax.nbin);

  /* Feed the initial solution*/
  par[0]=lffit.alfa;  /* Alpha */
  par[1]=lffit.Mstar;
  par[2]=log(lffit.phistar);
  sigpar[0]=10.*lffit.erralfa;
  sigpar[1]=10.*lffit.errMstar;
  sigpar[2]=10.*lffit.errphistar/lffit.phistar;
  
  /* Initial solution */
  printf("par0 %g par1 %g par2 %g\n",par[0],par[1],par[2]);

  printf(" Computing LF...\n");

  iter_amo=Amoeba_d(n,magDistn,z,3,par,sigpar,FTOL,MAXITER,Amoe_Funk_STY_p_M_wC_main);
/*   printf(" FINAL par0 %.15g par1 %.15g \n",par[0],par[1]); */

  /* Info that will be output in mlinfo */
  _MLmax_STY_p_M_wC=Amoe_Funk_STY_p_M_wC_main(n,magDistn,z,par);
  mlinfo->nIter = iter_amo;
  mlinfo->MLmax = _MLmax_STY_p_M_wC;

  /* Meto la solucion en la salida */
 
  lf->alfa=par[0];
  lf->Mstar=par[1];
  lf->phistar=exp(par[2]);

  if(DEBUG) printf(" Solucion sin errores: Mstar %.15g alpha %.4g ",lf->Mstar,lf->alfa);  

  /* Estimacion de los errores en los parametros */

  NumericalHessianCovars_STY_p_M_wC(n,magDistn,z,par,sigpar,mlim,cosmo,lf); 
  if(DEBUG) printf(" Solucion final: Mstar %.15g +/- %.15g alpha %.4g +/- %.4g\n",lf->Mstar,lf->errMstar,lf->alfa,lf->erralfa);  
  
  free(y);
  if(DEBUG) printf(" MLF %g\n",_MLmax_STY_p_M_wC); 

  if(iter_amo>=MAXITER-1) return(2);
  gsl_rng_free(_random_gen_STY_gmz_p_f_M_wC);
  return(0);
}

void prepareGlobalVars_STY_p_M_wC(double *z, double *magn)
{
  size_t i;

   _Mabsn_STY_p_M_wC = vector_d(_ndata_STY_p_M_wC);
  for(i=0;i<_ndata_STY_p_M_wC;i++) _Mabsn_STY_p_M_wC[i]=Mag(z[i],magn[i],*_cosmo_STY_p_M_wC);
}

double Amoe_Funk_STY_p_M_wC_main(int n, double *x, double *y, double *p)
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

  double colori;
  double logColor;

  /* Variables for the VEGAS integral */
  double probabajo, errprobabajo;
  double xl_den[1], xu_den[1];
  size_t dim_den=1;
  size_t calls_den = 5000;
  
  gsl_monte_function G_den = 
    {&vegas_funk_denominator_STY_p_M_wC, dim_den, 0};
  

  lfamo.alfa=p[0];
  lfamo.Mstar=p[1];
  lfamo.phistar=exp(p[2]);
  _lf_STY_p_M_wC=&lfamo;

  logL=0.;

  Lstar=pow(10.,-0.4*lfamo.Mstar);

  for(i=0;i<_ndata_STY_p_M_wC;i++)
  {
    colori=_magDistn_STY_p_M_wC[i] - _magSeln_STY_p_M_wC[i];
    /* gcolor=gaussian(colori, _color_mean_STY_p_M_wC, _color_stddev_STY_p_M_wC); */
    if (_color_stddev_STY_p_M_wC==0)
    {
      logColor=0.0;
    }
    else
    {
      logColor=lngaussian(colori, _color_mean_STY_p_M_wC, _color_stddev_STY_p_M_wC);
    }
    /* Mabs=Mag(y[i],x[i],*_cosmo_STY_p_M_wC); */
    Mabs=_Mabsn_STY_p_M_wC[i];
    Mlow=Mag(y[i],_mlim_STY_p_M_wC+colori,*_cosmo_STY_p_M_wC);
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
      if(DEBUG2) printf("log_gamma_int podía petar así que GSL_LOG_DBL_MIN\n");
    }
    else
    {
      log_gamma_int=log(gsl_sf_gamma_inc(1+lfamo.alfa,Llow/Lstar));
      if(DEBUG2) printf("log_gamma_int: lfamo.alfa %g Llow %g Lstar %g\n",lfamo.alfa, Llow, Lstar);
    }
    _z_i_STY_p_M_wC=y[i];
    _magSel_i_STY_p_M_wC=_magSeln_STY_p_M_wC[i];

    xl_den[0]=_color_mean_STY_p_M_wC-6*_color_stddev_STY_p_M_wC;
    xu_den[0]=_color_mean_STY_p_M_wC+6*_color_stddev_STY_p_M_wC;
    probabajo = vegas_integrate_STY_p_M_wC
        (&G_den, xl_den, xu_den, dim_den, calls_den, &errprobabajo);
    /* log(lfamo.phistar) + log_gamma_int -> integral de la función de Schecter
    entre Llow e inf */
    if (DEBUG2) printf("logColor: %g  colori: %g\n",logColor,colori);
    if (DEBUG2) printf("logSch %g log(phistar) %g loggamm %g\n",log(Schechter_M(Mabs,lfamo)),log(lfamo.phistar),log_gamma_int);
    if (DEBUG2) printf("_color_mean: %g _color_stddev: %g\n",_color_mean_STY_p_M_wC,_color_stddev_STY_p_M_wC);
    logL-= log(Schechter_M(Mabs,lfamo)) + logColor - probabajo;
    //logL-= log(Schechter_M(Mabs,lfamo)) + logColor - log(lfamo.phistar) - log_gamma_int;
    if (DEBUG2) printf("BUCLE logL: %g\n",logL);
    if (DEBUG2) printf("#LSCH %g %g %g\n",log(Schechter_M(Mabs,lfamo)),log(lfamo.phistar),log_gamma_int);
    if (DEBUG2) printf("iobj %d Mabs %g z %g\n",i,Mabs,y[i]); 
  }

  if(DEBUG) printf(" logL %g\n",logL);

  /* Aqui viene la parte de la poissoniana de npob */
  if (_color_stddev_STY_p_M_wC==0) /* to avoid problem in Int_sch_M_wC */
  {
    /* printf("Using Int_sch_M instead of Int_sch_M_wC because color_stddev = 0\n"); */
    Ntot=Int_sch_M(lfamo,_zlow_STY_p_M_wC,_zup_STY_p_M_wC,_mlim_STY_p_M_wC,*_cosmo_STY_p_M_wC)*_strrad_STY_p_M_wC/4./M_PI;
  }
  else
  {
    Ntot=Int_sch_M_wC(lfamo,_zlow_STY_p_M_wC,_zup_STY_p_M_wC,_color_mean_STY_p_M_wC,_color_stddev_STY_p_M_wC,_mlim_STY_p_M_wC,*_cosmo_STY_p_M_wC)*_strrad_STY_p_M_wC/4./M_PI;
  }

  logL-= (_ndata_STY_p_M_wC*log(Ntot) - Ntot - gammln((double)_ndata_STY_p_M_wC+1.));

  _iter_m_STY_p_M_wC++;


  if(DEBUG) printf(" Ntot %g ndata %d rad %g zlw %f zup %f mag %f lfamophi %f\n",Ntot, _ndata_STY_p_M_wC, _strrad_STY_p_M_wC,_zlow_STY_p_M_wC,_zup_STY_p_M_wC,_mlim_STY_p_M_wC,lfamo.phistar);

  if(DEBUG) printf(" Iter %d  logL %f par0 %g par1 %g par2 %g (%.10g)\n",_iter_m_STY_p_M_wC,logL,p[0],p[1],exp(p[2]),p[2]);
  return(logL);
}
double vegas_funk_denominator_STY_p_M_wC (double *x, size_t dim, void *params)
{
  double res;
  double zreal;
  double mDistreal;
  double colorobs;
  static double dMDistdmDist=1.;
  double logfacLF;
  double logfacColor;
  double Mlow, Llow, Lstar;
  double log_gamma_int;

  if(DEBUG3) printf(" inside vegas_funk_denominator.\n");

  /*
  x[0] -> mdist - mobssel
  */

  colorobs  = x[0];
  zreal     = _z_i_STY_p_M_wC;
  mDistreal = _magSel_i_STY_p_M_wC + colorobs;

  Mlow=Mag(zreal,_mlim_STY_p_M_wC+colorobs,*_cosmo_STY_p_M_wC);
  Llow=pow(10.,-0.4*Mlow);
  Lstar=pow(10.,-0.4*_lf_STY_p_M_wC->Mstar);
  if(Llow/Lstar > 0.25 && (_lf_STY_p_M_wC->alfa*log(Llow/Lstar) - Llow/Lstar) <= GSL_LOG_DBL_MIN)
  {
    log_gamma_int=GSL_LOG_DBL_MIN;
    if(DEBUG2) printf("log_gamma_int podía petar así que GSL_LOG_DBL_MIN\n");
  }
  else
  {
    log_gamma_int=log(gsl_sf_gamma_inc(1+_lf_STY_p_M_wC->alfa,Llow/Lstar));
    if(DEBUG2) printf("log_gamma_int: lfamo.alfa %g Llow %g Lstar %g\n",_lf_STY_p_M_wC->alfa, Llow, Lstar);
  }
  
  //logfacLF = log(Schechter_M(Mabs,*_lf_STY_gmz_p_f_M_wC));
  logfacLF = log_gamma_int + log(_lf_STY_p_M_wC->phistar);
  logfacColor = lngaussian(colorobs, _color_mean_STY_p_M_wC, _color_stddev_STY_p_M_wC);
  res=exp(logfacLF+logfacColor);
  res=res*dMDistdmDist;

  return(res); 
}



/* Errores utilizando derivadas numéricas con GSL */
void NumericalHessianCovars_STY_p_M_wC(int n,double *magn,double *z,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf)
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

  struct Amoe_Funk_gsl_param_STY_p_M_wC deriv_param;
  
  deriv_param.nData = _ndata_STY_p_M_wC;
  deriv_param.magn = magn;
  deriv_param.z = z;

  gsl_multimin_function LogLFunction;
  LogLFunction.f = &Amoe_Funk_STY_p_M_wC_main_gsl_multimin;
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

double Amoe_Funk_STY_p_M_wC_main_gsl_multimin(const gsl_vector *x, void *params)
{
  int nData           = ((struct Amoe_Funk_gsl_param_STY_p_M_wC *)params)->nData;
  double* magn        = ((struct Amoe_Funk_gsl_param_STY_p_M_wC *)params)->magn;
  double* z           = ((struct Amoe_Funk_gsl_param_STY_p_M_wC *)params)->z;
  double value;
  double* lf_param = vector_d(3);

  lf_param[0] = gsl_vector_get(x,0);
  lf_param[1] = gsl_vector_get(x,1);
  lf_param[2] = gsl_vector_get(x,2);

  value=Amoe_Funk_STY_p_M_wC_main(nData, magn, z, lf_param);
  free(lf_param);
  return value;
}

double vegas_integrate_STY_p_M_wC(gsl_monte_function * f,
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

