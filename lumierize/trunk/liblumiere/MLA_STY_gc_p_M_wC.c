#include <gsl/gsl_machine.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include "mlsty.h"
#include "gsl_hessian.h"
#include "amoeba.h"
#include "alloc.h"
#include "fermisel.h"
#include "functions.h"
#include "step.h"
#include "gaussj.h"
#include "minmax.h"
#include "vvmax.h"
#include "gaussint.h"


/* VEGAS */
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte_vegas.h>

#define FTOL  1e-10
#define FTOL2 1e-6
#define FTOL3 1e-7
#define MAXITER  1000
#define MAXITER2 120
#define VERBOSE 0
#define DEBUG  1
#define DEBUG2 1
#define DEBUG3 1
#define DEBUGPLOT 0
#define TOLERR 0.07 

#define VEGAS 1

/* Este contiene STY con colores y errores en los colores */

/* #define TOLERR 0.0001 */
/* Estructura para contener los par�metros de Amo..main_gsl_multimin */
struct Amoe_Funk_gsl_param_STY_gc_p_M_wC
{
  int nData;
  double *magn;
  double *z;
};

void   PrepareGlobalVars_STY_gc_p_M_wC(double *z, double *magn);
double Amoe_Funk_STY_gc_p_M_wC_main(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_gc_p_M_wC_main_gsl_multimin(const gsl_vector *x, void *params);
double Funk2_intmag_STY_gc_p_M_wC(double fluxreal);
double Funk1_intMag_STY_gc_p_M_wC(double x);
void   NumericalHessianCovars_STY_gc_p_M_wC(int n,double *magDistn,double *errColorn, double *z, double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf);

/* VEGAS */
double Amoe_Funk_STY_gc_p_M_wC_main_VEGAS(int n, double *x, double *y, double *p);
double f1 (double *k1, size_t dim1, void *params1);
double f2 (double *k2, size_t dim2, void *params2);

struct Schlf_M *_lf_STY_gc_p_M_wC;
struct cosmo_param *_cosmo_STY_gc_p_M_wC;
double _mlim_STY_gc_p_M_wC;
double _z_i_STY_gc_p_M_wC;

int _ndata_STY_gc_p_M_wC;
int _iter_m_STY_gc_p_M_wC;
int _iter_c_STY_gc_p_M_wC;
double _MLmax_STY_gc_p_M_wC;
double _zlow_STY_gc_p_M_wC;
double _zup_STY_gc_p_M_wC;
double _strrad_STY_gc_p_M_wC;
double *_magSeln_STY_gc_p_M_wC;
double *_magDistn_STY_gc_p_M_wC;
double _magDistn_i_STY_gc_p_M_wC;
double _magSeln_i_STY_gc_p_M_wC;
double _color_mean_STY_gc_p_M_wC;
double _color_stddev_STY_gc_p_M_wC;
double *_Mabsn_STY_gc_p_M_wC;
double *_errColorn_STY_gc_p_M_wC;
double _errColorn_i_STY_gc_p_M_wC;
double _colorn_i_STY_gc_p_M_wC;

int MLA_STY_gc_p_M_wC(int n,double *magSeln, double *magDistn, double color_mean, double color_stddev, double *errColorn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf, struct MLProcessInfo *mlinfo)
{

  double par[3];
  double sigpar[3];
  int i;
  int iter_amo;
  int iter_amo_VEGAS;

  /* Variables to use for vvmax fitting */
  struct Steplf_M  lfvvmax;
  double minMabs,maxMabs;
  double *Mabs;
  double chisq;
  struct Schlf_M  lffit;

  Mabs=vector_d(n);

  /* Copias globales de estas variables para que se vean en toda la subrutina */
  
  _ndata_STY_gc_p_M_wC=n;
  _cosmo_STY_gc_p_M_wC=&cosmo;
  _mlim_STY_gc_p_M_wC=mlim;
  _strrad_STY_gc_p_M_wC=strrad;
  _zlow_STY_gc_p_M_wC=zlow;
  _zup_STY_gc_p_M_wC=zup;
  _magSeln_STY_gc_p_M_wC=magSeln;
  _magDistn_STY_gc_p_M_wC=magDistn;
  _color_mean_STY_gc_p_M_wC=color_mean;
  _color_stddev_STY_gc_p_M_wC=color_stddev;
  _errColorn_STY_gc_p_M_wC=errColorn;

  _iter_m_STY_gc_p_M_wC=0;

  if(DEBUG2)
  {
    for(i=0;i<n;i++)
    {
      printf(" Entrada x %g %g %g\n",magDistn[i],magSeln[i],errColorn[i]);
    }
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
  lfvvmax.ngalbin   =vector_i(lfvvmax.nbin);
  lfvvmax.covarlnlf =matrix_d(lfvvmax.nbin,lfvvmax.nbin);

  for(i=0;i<n;i++)   Mabs[i]=Mag(z[i],magDistn[i],cosmo);
  MinMax_d(n,Mabs,&minMabs,&maxMabs);
  for(i=0;i<=lfvvmax.nbin;i++) lfvvmax.magni[i]=minMabs+i*(maxMabs-minMabs)/lfvvmax.nbin;
  printf("Antes de VVmax\n");
  VVmax_M(n,magSeln,magDistn,z,mlim,strrad,zlow,zup,cosmo,&lfvvmax);
  /* VVmax_M(n,magDistn,magDistn,z,mlim,strrad,zlow,zup,cosmo,&lfvvmax); */
/*   cpgopen("?"); */
/*   PlotStepLF_M(lfvvmax); */
  if(DEBUG) printf(" Salida VVmax\n");
  if(DEBUG) for(i=0;i<lfvvmax.nbin;i++) printf(" Mabs %g - %g LF %g\n",lfvvmax.magni[i],lfvvmax.magni[i+1],lfvvmax.lnlf[i]/log(10));
  FitSch2StepLF_M(lfvvmax,&lffit, &chisq);
  if(DEBUG) {
    printf(" Despues ajuste MRQ\n");
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
  free(lfvvmax.ngalbin);
  free_matrix_d(lfvvmax.covarlnlf,lfvvmax.nbin,lfvvmax.nbin);

  /* dabreu */
  /* struct MLProcessInfo mlinfo2;
  i=MLA_STY_p_M(n,magn,z,mlim,strrad,zlow,zup,cosmo,&lffit, &mlinfo2);
  printf(" STY as a first solution.\n");
  printf(" STY -> Mstar: %g alpha: %g phistar: %g\n",lffit.Mstar,lffit.alfa,lffit.phistar); */

  printf(" Computing LF...\n");

  iter_amo=0;
  iter_amo_VEGAS=0;
  while(iter_amo==0) { 
    par[0]=lffit.alfa;  /* Alpha */
    par[1]=lffit.Mstar;
    par[2]=log(lffit.phistar);

    sigpar[0]=10.*lffit.erralfa;
    sigpar[1]=10.*lffit.errMstar;
    sigpar[2]=10.*lffit.errphistar/lffit.phistar;

    if(VEGAS)
    {
      printf("Entering VEGAS\n");
      iter_amo_VEGAS=Amoeba_d(n,magDistn,z,3,par,sigpar,FTOL,MAXITER,Amoe_Funk_STY_gc_p_M_wC_main_VEGAS);
    }

    iter_amo=Amoeba_d(n,magDistn,z,3,par,sigpar,FTOL,MAXITER,Amoe_Funk_STY_gc_p_M_wC_main);
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
  _MLmax_STY_gc_p_M_wC=Amoe_Funk_STY_gc_p_M_wC_main(n,magDistn,z,par);
  mlinfo->nIter = iter_amo;
  mlinfo->MLmax = _MLmax_STY_gc_p_M_wC;

  /* Meto la solucion en la salida */
 
  lf->alfa=par[0];
  lf->Mstar=par[1];
  lf->phistar=exp(par[2]);

  /* Estimacion de los errores en los parametros */

  printf(" Computing errors in LF parameters...\n");

  NumericalHessianCovars_STY_gc_p_M_wC(n,magDistn,errColorn,z,par,sigpar,mlim,cosmo,lf);
  if(DEBUG) printf(" Solucion final: Mstar %.15g +/- %.15g alpha %.4g +/- %.4g\n",lf->Mstar,lf->errMstar,lf->alfa,lf->erralfa);

  if(DEBUG) printf(" MLF %g\n",_MLmax_STY_gc_p_M_wC);

  /* Code exit status */
  if(iter_amo>=MAXITER-1) return(2);
  return(0);
}


double Amoe_Funk_STY_gc_p_M_wC_main(int n, double *x, double *y, double *p) {
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
  double offset;
  double scale;
  double x1,x2;
  int npa=21;
  int npb=21;
  /* int npc=31; */
  double probarriba;
  double probabajo;

  double intsup;
  (void)n;

  lfamo.alfa=p[0];
  lfamo.Mstar=p[1];
  lfamo.phistar=exp(p[2]);

  _lf_STY_gc_p_M_wC=&lfamo;

  logL=0.;
  logLold=0.;

  Lstar=pow(10.,-0.4*lfamo.Mstar);
  
  if(DEBUG3) printf(" Entra con  par0 %g par1 %g par2 %g\n",p[0],p[1],p[2]);

  intsup=incom(1+lfamo.alfa,100.);
  for(i=0;i<_ndata_STY_gc_p_M_wC;i++) 
  {
    _z_i_STY_gc_p_M_wC=y[i];
    _magDistn_i_STY_gc_p_M_wC=x[i];
    _magSeln_i_STY_gc_p_M_wC=_magSeln_STY_gc_p_M_wC[i];
    _colorn_i_STY_gc_p_M_wC=_magDistn_i_STY_gc_p_M_wC-_magSeln_i_STY_gc_p_M_wC;
    _errColorn_i_STY_gc_p_M_wC=_errColorn_STY_gc_p_M_wC[i];
    Mabs=Mag(y[i],x[i],*_cosmo_STY_gc_p_M_wC);
    Mlow=Mag(y[i],_mlim_STY_gc_p_M_wC+_colorn_i_STY_gc_p_M_wC,*_cosmo_STY_gc_p_M_wC);
    xmin =pow(10.,-0.4*Mlow)/Lstar;
    if(xmin>=100) /* cambiar por el GSL_LOG_DBL_MIN ? */
    {
      logL+= GSL_LOG_DBL_MAX;
    }
    else 
    {
      x1=lfamo.Mstar-10.0; /* 5 por poner un n�mero (deber�a ser -inf) */
      x2=Mag(y[i],_mlim_STY_gc_p_M_wC,*_cosmo_STY_gc_p_M_wC);
      probabajo=gaussintleg_d(Funk1_intMag_STY_gc_p_M_wC,x1,x2,npb);
      if(DEBUG2) printf(" Primer abajo %g con 1000 %g x1 %f x2 %f\n",probabajo,gaussintleg_d(Funk1_intMag_STY_gc_p_M_wC,x1-10,x2,1000),x1,x2);
/*      x2=Mag(y[i],_mlim_STY_gc_p_M_wC+6*_errColorn_i_STY_gc_p_M_wC,*_cosmo_STY_gc_p_M_wC);
      x1=Mag(y[i],_mlim_STY_gc_p_M_wC-6*_errColorn_i_STY_gc_p_M_wC,*_cosmo_STY_gc_p_M_wC);
      probabajo+=gaussintleg_d(Funk1_intMag_STY_gc_p_M_wC,x1,x2,npc);
      if(DEBUG2) printf(" Segundo abajo %g con 1000 %g x1 %f x2 %f\n",gaussintleg_d(Funk1_intMag_STY_gc_p_M_wC,x1,x2,npc),gaussintleg_d(Funk1_intMag_STY_gc_p_M_wC,x1,Mag(y[i],_mlim_STY_gc_p_M_wC+12*_errColorn_i_STY_gc_p_M_wC,*_cosmo_STY_gc_p_M_wC),1000),x1,x2); */
      /* La integral de abajo se hac�a en dos intervalos, ahora vamos a probar de un tir�n */

      if(DEBUG2) printf(" Calculo abajo %g old %g \n",probabajo,lfamo.phistar*(incom(1+lfamo.alfa,200.)-incom(1+lfamo.alfa,pow(10.,-0.4*Mlow)/Lstar)));

      /* We have to compute them again because Funk? may have overriden them */
      _z_i_STY_gc_p_M_wC=y[i];
      _magDistn_i_STY_gc_p_M_wC=x[i];
      _magSeln_i_STY_gc_p_M_wC=_magSeln_STY_gc_p_M_wC[i];
      _colorn_i_STY_gc_p_M_wC=_magDistn_i_STY_gc_p_M_wC-_magSeln_i_STY_gc_p_M_wC;
      _errColorn_i_STY_gc_p_M_wC=_errColorn_STY_gc_p_M_wC[i];

      if(_errColorn_i_STY_gc_p_M_wC==0)
      {
        probarriba=Schechter_M(Mabs,*_lf_STY_gc_p_M_wC); /* FIXME: esto seguro que es as�? */
      }
      else
      {
        offset=_magDistn_i_STY_gc_p_M_wC;
        scale=_errColorn_i_STY_gc_p_M_wC*sqrt(2.);
        /* x2=_mlim_STY_gc_p_M_wC;
        if(_mlim_STY_gc_p_M_wC>(_magDistn_i_STY_gc_p_M_wC+6*_errColorn_i_STY_gc_p_M_wC)) x2=_magDistn_i_STY_gc_p_M_wC+6*_errColorn_i_STY_gc_p_M_wC; 
        x1=_magDistn_i_STY_gc_p_M_wC-6*_errColorn_i_STY_gc_p_M_wC; 
        if(_mlim_STY_gc_p_M_wC<(_magDistn_i_STY_gc_p_M_wC-5*_errColorn_i_STY_gc_p_M_wC)) x1=x2-6*_errColorn_i_STY_gc_p_M_wC; */
        x1=_colorn_i_STY_gc_p_M_wC-6*_errColorn_i_STY_gc_p_M_wC; 
        x2=_colorn_i_STY_gc_p_M_wC+6*_errColorn_i_STY_gc_p_M_wC;
        probarriba=gaussintleg_d(Funk2_intmag_STY_gc_p_M_wC,x1,x2,npa);
        if(DEBUG2) printf(" parriba %g  con 1000 %g con her %g\n",probarriba,gaussintleg_d(Funk2_intmag_STY_gc_p_M_wC,x1,x2,1000),gaussinther_d(Funk2_intmag_STY_gc_p_M_wC,offset,scale,100));
      }

      if(DEBUG2) printf(" Calculo arriba %g old %g x1 %g  x2 % g magn %g err %g magnl %g\n",probarriba,Schechter_M(Mabs,*_lf_STY_gc_p_M_wC),x1,x2,_magDistn_i_STY_gc_p_M_wC,_errColorn_i_STY_gc_p_M_wC,_mlim_STY_gc_p_M_wC);


      if(probarriba==0 || probabajo==0) logL+=GSL_LOG_DBL_MAX; 
      logL-= log(probarriba) - log(probabajo);
      if(DEBUG3) printf(" iobj %d logL %f loglold %f      sch %g    pa %g pb %g (%g)  xmin %g x1 %g x2 %g magn %g err %g\n",i,logL,logLold,Schechter_M(Mabs,lfamo),log(probarriba),log(probabajo),probabajo,xmin,x1,x2,_magDistn_i_STY_gc_p_M_wC,_errColorn_i_STY_gc_p_M_wC);
/*       printf(" iobj %d logL %f loglold %f      sch %g int %g (%g)   pa %g pb %g (%g)  xmin %g x1 %f x2 %f\n",i,logL,logLold,Schechter_M(Lumi,lfamo),log(intsch),intsch,log(probarriba),log(probabajo),probabajo,xmin,x1,x2); */
    }
  }
  /* Aqui viene la parte de la poissoniana de npob */
  if (_color_stddev_STY_gc_p_M_wC==0.)
  {
    Ntot=Int_sch_M(lfamo,_zlow_STY_gc_p_M_wC,_zup_STY_gc_p_M_wC,_mlim_STY_gc_p_M_wC,*_cosmo_STY_gc_p_M_wC)*_strrad_STY_gc_p_M_wC/4./M_PI; 
  }
  else
  {
    Ntot=Int_sch_M_wC(lfamo,_zlow_STY_gc_p_M_wC,_zup_STY_gc_p_M_wC,_color_mean_STY_gc_p_M_wC,_color_stddev_STY_gc_p_M_wC,_mlim_STY_gc_p_M_wC,*_cosmo_STY_gc_p_M_wC)*_strrad_STY_gc_p_M_wC/4./M_PI; 
  }
  logL-=    (_ndata_STY_gc_p_M_wC*log(Ntot) - Ntot - gammln((double)_ndata_STY_gc_p_M_wC+1.)); 
/*   logLold-= (ndata*log(Ntot) - Ntot - gammln((double)ndata)+1.);  */
  
  if(DEBUG2) printf(" NTOT %f ndata*log(Ntot) %f gamm %f\n",Ntot,_ndata_STY_gc_p_M_wC*log(Ntot),gammln((double)_ndata_STY_gc_p_M_wC)+1.);

  _iter_m_STY_gc_p_M_wC++;
  if(DEBUG) printf(" Iter %d  logL %f logLold %f par0 %g par1 %g par2 %g\n",_iter_m_STY_gc_p_M_wC,logL,logLold,p[0],p[1],p[2]);
  return(logL); 
}

double Funk1_intMag_STY_gc_p_M_wC(double Mabs) 
{
  /* Lo suyo es que sea en magnitudes absolutas, ya que 
     la FL est� expresada en intervalos de Mabs. 
     Si no, tendr�a que multiplicar por dMabs/dmag */

  int npa=21;
/*   double scale,offset; */
  double x1,x2;
  double firstint;

  _magDistn_i_STY_gc_p_M_wC=mag(_z_i_STY_gc_p_M_wC,Mabs,*_cosmo_STY_gc_p_M_wC);
  /* x2=_mlim_STY_gc_p_M_wC; 
  if(_mlim_STY_gc_p_M_wC>_magDistn_i_STY_gc_p_M_wC+6*_errColorn_i_STY_gc_p_M_wC) x2=_magDistn_i_STY_gc_p_M_wC+6*_errColorn_i_STY_gc_p_M_wC; */

  x1=_colorn_i_STY_gc_p_M_wC-6*_errColorn_i_STY_gc_p_M_wC; 
  x2=_colorn_i_STY_gc_p_M_wC+6*_errColorn_i_STY_gc_p_M_wC;

  /* la integral ser� entre color-6*errcolor y color+6*errcolor */
  /* los l�mites se chequean para que la magSel de aplicar este color a la magDist no se salga de mlim */
 
  /* if(_mlim_STY_gc_p_M_wC<_magDistn_i_STY_gc_p_M_wC-5*_errColorn_i_STY_gc_p_M_wC) x1=x2-6*_errColorn_i_STY_gc_p_M_wC; */
  if(_errColorn_i_STY_gc_p_M_wC==0)
  {
    firstint=Schechter_M(Mabs,*_lf_STY_gc_p_M_wC); /* FIXME: seguro que as� no */
  }
  else
  { 
    firstint=gaussintleg_d(Funk2_intmag_STY_gc_p_M_wC,x1,x2,npa);
  }
/*   printf(" Segunda integral %g\n",firstint); */
/*   printf(" El bueno %g y el sch %g lum %g\n",log(firstint),Schechter_M(Lumi,*lf_STY_gc_p_M_wC),lumlog);  */
/*   printf(" Sale con %g\n",firstint*Lumi); */
  return(firstint); /* Este producto es para tener en cuenta que hacemos la integral en log(flux) */
}

double Funk2_intmag_STY_gc_p_M_wC(double colornreal) 
{
  double Mabs;
  double logfacLF,logfacerr, logColor;
  double magDistnreal;
  magDistnreal=colornreal+_magSeln_i_STY_gc_p_M_wC;
  if((_magDistn_i_STY_gc_p_M_wC-colornreal)>_mlim_STY_gc_p_M_wC) return(0); /* mirar si esto est� bien */
  else 
  {
    Mabs=Mag(_z_i_STY_gc_p_M_wC,magDistnreal,*_cosmo_STY_gc_p_M_wC);
    logfacLF = log(Schechter_M(Mabs,*_lf_STY_gc_p_M_wC));
    if (_color_stddev_STY_gc_p_M_wC==0.)
    {
      logColor=0.0;
    }
    else
    {
      logColor = lngaussian(colornreal, _color_mean_STY_gc_p_M_wC, _color_stddev_STY_gc_p_M_wC);
    }
    /* if (_errColorn_i_STY_gc_p_M_wC==0.) este if ya se hace en Funk1
    {
      logfacerr=0.0;
    }
    else
    { */
    logfacerr= lngaussian(colornreal,_colorn_i_STY_gc_p_M_wC,_errColorn_i_STY_gc_p_M_wC);
      /* no estoy seguro de si va primero colornreal, pero la funci�n es sim�trica respecto a los dos primeros par�metros */
    /* } */
    return(exp(logfacLF+logfacerr+logColor)); /* aqu� va +logColor ? */
  }
}

/* Errores utilizando derivadas num�ricas con GSL */

void NumericalHessianCovars_STY_gc_p_M_wC(int n,double *magDistn,double *errColorn,double *z,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf)
{
  gsl_matrix *gsl_hessian;
  gsl_vector *param;
  double **covar;
  double **bb;
  size_t i,j;

  double hessianStep = GSL_ROOT4_DBL_EPSILON;
  struct Amoe_Funk_gsl_param_STY_gc_p_M_wC deriv_param;
  gsl_multimin_function LogLFunction;
  (void)cosmo;
  (void)mlim;
  (void)n;
  (void)sigpar;
  (void)errColorn;
  /* derivStep = 0.01; */

  gsl_hessian=gsl_matrix_alloc(3,3);
  param=gsl_vector_alloc(3);
  gsl_vector_set(param, 0, par[0]);
  gsl_vector_set(param, 1, par[1]);
  gsl_vector_set(param, 2, par[2]);

  covar=matrix_d(3 ,3 );
  bb=matrix_d(3,1);


  deriv_param.nData = _ndata_STY_gc_p_M_wC;
  deriv_param.magn = magDistn;
  deriv_param.z = z;

  LogLFunction.f = &Amoe_Funk_STY_gc_p_M_wC_main_gsl_multimin;
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

double Amoe_Funk_STY_gc_p_M_wC_main_gsl_multimin(const gsl_vector *x, void *params)
{
  int nData           = ((struct Amoe_Funk_gsl_param_STY_gc_p_M_wC *)params)->nData;
  double* magn        = ((struct Amoe_Funk_gsl_param_STY_gc_p_M_wC *)params)->magn;
  double* z           = ((struct Amoe_Funk_gsl_param_STY_gc_p_M_wC *)params)->z;
  double value;
  double* lf_param = vector_d(3);

  lf_param[0] = gsl_vector_get(x,0);
  lf_param[1] = gsl_vector_get(x,1);
  lf_param[2] = gsl_vector_get(x,2);


  value=Amoe_Funk_STY_gc_p_M_wC_main(nData, magn, z, lf_param);
  free(lf_param);
  return value;
}

/* with VEGAS integration */
double Amoe_Funk_STY_gc_p_M_wC_main_VEGAS(int n, double *x, double *y, double *p)
{
  int i; 
  double logL=0.;
  double logLold=0.;
  struct Schlf_M lfamo;
  double Lstar;
  double Mabs;
  double Mlow;
  double xmin;
  double Ntot;
  double offset;
  double scale;
  double x1,x2;
  /* int npa=21;
  int npb=21; */
  double probarriba;
  double probabajo;

  /* VEGAS */
  const gsl_rng_type *T;
  gsl_rng *r1;
  gsl_rng *r2;
  double res, err;
  double xl1[2];
  double xu1[2];
  double xl2[1];
  double xu2[1];
  size_t calls=100;
  gsl_monte_function F1 = {&f1, 2, p};
  gsl_monte_function F2 = {&f2, 1, p};
  (void)n;


  gsl_rng_env_setup();
  T = gsl_rng_default;
  r1 = gsl_rng_alloc(T);
  r2 = gsl_rng_alloc(T);
  /* falta definir la f1 y la f2 -> parecidas a Funk1 y Funk2 */


  lfamo.alfa=p[0];
  lfamo.Mstar=p[1];
  lfamo.phistar=exp(p[2]);

  _lf_STY_gc_p_M_wC=&lfamo;

  logL=0.;
  logLold=0.;

  Lstar=pow(10.,-0.4*lfamo.Mstar);

  printf("Inside VEGAS\n");  
  if(DEBUG3) printf(" Entra con  par0 %g par1 %g par2 %g\n",p[0],p[1],p[2]);

  for(i=0;i<_ndata_STY_gc_p_M_wC;i++) 
  {
    _z_i_STY_gc_p_M_wC=y[i];
    _magDistn_i_STY_gc_p_M_wC=x[i];
    _magSeln_i_STY_gc_p_M_wC=_magSeln_STY_gc_p_M_wC[i];
    _colorn_i_STY_gc_p_M_wC=_magDistn_i_STY_gc_p_M_wC-_magSeln_i_STY_gc_p_M_wC;
    _errColorn_i_STY_gc_p_M_wC=_errColorn_STY_gc_p_M_wC[i];
    Mabs=Mag(y[i],x[i],*_cosmo_STY_gc_p_M_wC);
    Mlow=Mag(y[i],_mlim_STY_gc_p_M_wC+_colorn_i_STY_gc_p_M_wC,*_cosmo_STY_gc_p_M_wC);
    xmin =pow(10.,-0.4*Mlow)/Lstar;
    if(xmin>=100) 
    {
      logL+=GSL_LOG_DBL_MAX; 
    }
    else 
    {
      gsl_monte_vegas_state *s1 = gsl_monte_vegas_alloc(2);
      x1=lfamo.Mstar-10.0; /* 5 por poner un n�mero (deber�a ser -inf) */
      x2=Mag(y[i],_mlim_STY_gc_p_M_wC,*_cosmo_STY_gc_p_M_wC);
      xl1[0]=x1;
      xu1[0]=x2;
      /* falta xl1[1] y xu1[1] */
      gsl_monte_vegas_integrate (&F1, xl1, xu1, 2, calls/10., r1, s1, &res, &err); /* warm up */
      if(DEBUG2) printf("Warm up in VEGAS %g %g",res,err);
      gsl_monte_vegas_integrate (&F1, xl1, xu1, 2, calls, r1, s1, &res, &err);
      probabajo=res;

      gsl_monte_vegas_free(s1);
      gsl_rng_free(r1);

      if(DEBUG2) printf(" Primer abajo %g con 1000 %g x1 %f x2 %f\n",probabajo,gaussintleg_d(Funk1_intMag_STY_gc_p_M_wC,x1-10,x2,1000),x1,x2);

      if(DEBUG2) printf(" Calculo abajo %g old %g \n",probabajo,lfamo.phistar*(incom(1+lfamo.alfa,200.)-incom(1+lfamo.alfa,pow(10.,-0.4*Mlow)/Lstar)));

      /* We have to compute them again because Funk? may have overriden them */
      _z_i_STY_gc_p_M_wC=y[i];
      _magDistn_i_STY_gc_p_M_wC=x[i];
      _magSeln_i_STY_gc_p_M_wC=_magSeln_STY_gc_p_M_wC[i];
      _colorn_i_STY_gc_p_M_wC=_magDistn_i_STY_gc_p_M_wC-_magSeln_i_STY_gc_p_M_wC;
      _errColorn_i_STY_gc_p_M_wC=_errColorn_STY_gc_p_M_wC[i];

      if(_errColorn_i_STY_gc_p_M_wC==0)
      {
        probarriba=Schechter_M(Mabs,*_lf_STY_gc_p_M_wC); /* FIXME: esto seguro que es as�? */
      }
      else
      {
        gsl_monte_vegas_state *s2 = gsl_monte_vegas_alloc(1);
        offset=_magDistn_i_STY_gc_p_M_wC;
        scale=_errColorn_i_STY_gc_p_M_wC*sqrt(2.);
        x1=_colorn_i_STY_gc_p_M_wC-6*_errColorn_i_STY_gc_p_M_wC; 
        x2=_colorn_i_STY_gc_p_M_wC+6*_errColorn_i_STY_gc_p_M_wC;

        xl2[0]=x1;
        xu2[0]=x2;

        gsl_monte_vegas_integrate (&F2, xl2, xu2, 1, calls/10., r2, s2, &res, &err); /* warm up */
        if(DEBUG2) printf("Warm up in VEGAS %g %g",res,err);
        gsl_monte_vegas_integrate (&F2, xl2, xu2, 1, calls, r2, s2, &res, &err);

        probarriba=res;

        gsl_monte_vegas_free(s2);
        gsl_rng_free(r2);

        if(DEBUG2) printf(" parriba %g  con 1000 %g con her %g\n",probarriba,gaussintleg_d(Funk2_intmag_STY_gc_p_M_wC,x1,x2,1000),gaussinther_d(Funk2_intmag_STY_gc_p_M_wC,offset,scale,100));
      }

      if(DEBUG2) printf(" Calculo arriba %g old %g x1 %g  x2 % g magn %g err %g magnl %g\n",probarriba,Schechter_M(Mabs,*_lf_STY_gc_p_M_wC),x1,x2,_magDistn_i_STY_gc_p_M_wC,_errColorn_i_STY_gc_p_M_wC,_mlim_STY_gc_p_M_wC);


      if(probarriba==0 || probabajo==0) logL+=GSL_LOG_DBL_MAX; 
      logL-= log(probarriba) - log(probabajo);
      if(DEBUG3) printf(" iobj %d logL %f loglold %f      sch %g    pa %g pb %g (%g)  xmin %g x1 %g x2 %g magn %g err %g\n",i,logL,logLold,Schechter_M(Mabs,lfamo),log(probarriba),log(probabajo),probabajo,xmin,x1,x2,_magDistn_i_STY_gc_p_M_wC,_errColorn_i_STY_gc_p_M_wC);
    }
  }
  /* Aqui viene la parte de la poissoniana de npob */
  if (_color_stddev_STY_gc_p_M_wC==0.)
  {
    Ntot=Int_sch_M(lfamo,_zlow_STY_gc_p_M_wC,_zup_STY_gc_p_M_wC,_mlim_STY_gc_p_M_wC,*_cosmo_STY_gc_p_M_wC)*_strrad_STY_gc_p_M_wC/4./M_PI; 
  }
  else
  {
    Ntot=Int_sch_M_wC(lfamo,_zlow_STY_gc_p_M_wC,_zup_STY_gc_p_M_wC,_color_mean_STY_gc_p_M_wC,_color_stddev_STY_gc_p_M_wC,_mlim_STY_gc_p_M_wC,*_cosmo_STY_gc_p_M_wC)*_strrad_STY_gc_p_M_wC/4./M_PI; 
  }
  logL-=    (_ndata_STY_gc_p_M_wC*log(Ntot) - Ntot - gammln((double)_ndata_STY_gc_p_M_wC+1.)); 
  
  if(DEBUG2) printf(" NTOT %f ndata*log(Ntot) %f gamm %f\n",Ntot,_ndata_STY_gc_p_M_wC*log(Ntot),gammln((double)_ndata_STY_gc_p_M_wC)+1.);

  _iter_m_STY_gc_p_M_wC++;
  if(DEBUG) printf(" Iter %d  logL %f logLold %f par0 %g par1 %g par2 %g\n",_iter_m_STY_gc_p_M_wC,logL,logLold,p[0],p[1],p[2]);
  return(logL); 
}

double f1 (double *k1, size_t dim1, void *params1)
{
  double firstint=0.;
  (void)dim1;
  (void)params1;
  firstint=Funk1_intMag_STY_gc_p_M_wC(k1[0]);
  return firstint;
}
double f2 (double *k2, size_t dim2, void *params2)
{
  double firstint=0.;
  (void)dim2;
  (void)params2;
  Funk2_intmag_STY_gc_p_M_wC(k2[0]);
  return firstint;
}

