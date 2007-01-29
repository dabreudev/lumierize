#include "modulos.h"
/* dabreu */
#include <gsl/gsl_machine.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_deriv.h>

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

#ifdef TRYEMPIRICAL
#else 
#define TRYEMPIRICAL 0
#endif

#ifdef TRYNUMERICAL
#else 
#define TRYNUMERICAL 1
#endif

double probando_deriv(int n, double *x, double *y, double *p);


/* #define TOLERR 0.0001 */

/* Ahora mismo esta con el pcut definido, de modo que se aleja de valores 
   cercanos a 0. Puede afectar al calculo de las covarianzas!! */

/* Estructura para contener los parámetros */
struct num_deriv_param_STY_p_M
{
  int nData;
  int iCovar;
  int jCovar;
  int iDeriv;
  double *magn;
  double *z;
  double *lf_param;
};

double Amoe_Funk_STY_p_M_main(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_p_M_norm(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_p_M_prior(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_p_M_conf(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_p_M_main_gsl(double x, void *params);
double derivAmoe_Funk_STY_p_M_main_gsl(double x, void *params);
double Funk2_int_STY_p_M(double x);
double Funk1_norm_STY_p_M(double x);
void   EmpiricalCovars_STY_p_M(int n,double *magn,double *z,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf);
void   NumericalDerivCovars_STY_p_M(int n,double *magn,double *z,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf);

void   ComputeNorma_STY_p_M(int n, double mlim, double strrad, double zlow, double zup, struct Schlf_M *lf);

struct cosmo_param *_cosmo_STY_p_M;
double _mlim_STY_p_M;

int _ndata_STY_p_M;
int _iter_m_STY_p_M;
int _iter_c_STY_p_M;
int _nconfl_STY_p_M;
double _conflim_STY_p_M;
/* double ML[5*MAXTRIES*MAXITER]; */
double *_pp_STY_p_M;
double _MLmax_STY_p_M;
double _xf_STY_p_M,_Tf_STY_p_M;
double _sigi_STY_p_M;
double _xtmp_STY_p_M;
double _zlow_STY_p_M,_zup_STY_p_M,_strrad_STY_p_M;


int  MLA_STY_p_M(int n,double *magn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf, struct MLProcessInfo *mlinfo)
{

  double *y;
  double par[3];
  double sigpar[3];
  int i;
  int iter_amo;

  /* Variables to use for vvmax fitting */
  struct Steplf_M  lfvvmax;
  double minMabs,maxMabs;
  double *Mabs;
  struct Schlf_M  lffit;
  double chisq;


  /* Copias globales de estas variables para que se vean en toda la subrutina */
  _ndata_STY_p_M=n;
  _cosmo_STY_p_M=&cosmo;
  _mlim_STY_p_M=mlim;
  _strrad_STY_p_M=strrad;
  _zlow_STY_p_M=zlow;
  _zup_STY_p_M=zup;

  y=vector_d(n);
  _iter_m_STY_p_M=0;

  iter_amo=MAXITER+1;

  /* Los límites en z se reajustan */
  _zlow_STY_p_M = (_zlow_STY_p_M < ZMIN ? ZMIN : _zlow_STY_p_M);


  if(DEBUG3) printf(" iter_amo %d\n",iter_amo);

  /* VVmax computing as a first solution */
  lfvvmax.nbin=(int)(n/3);
  if(lfvvmax.nbin>30) lfvvmax.nbin=30;

  lfvvmax.magni   =vector_d(lfvvmax.nbin+1);
  lfvvmax.errmagni=vector_d(lfvvmax.nbin+1);
  lfvvmax.lf      =vector_d(lfvvmax.nbin);
  lfvvmax.errlf   =vector_d(lfvvmax.nbin);
  lfvvmax.covarlf =matrix_d(lfvvmax.nbin,lfvvmax.nbin);

  Mabs            =vector_d(n);
  for(i=0;i<n;i++)   Mabs[i]=Mag(z[i],magn[i],cosmo);
  MinMax_d(n,Mabs,&minMabs,&maxMabs);
  for(i=0;i<=lfvvmax.nbin;i++) lfvvmax.magni[i]=minMabs+i*(maxMabs-minMabs)/lfvvmax.nbin;
  VVmax_M(n,magn,z,mlim,strrad,zlow,zup,cosmo,&lfvvmax);

  /* Fit of the VVmax solution to a Schechter function */
  FitSch2StepLF_M(lfvvmax,&lffit, &chisq);

  /* Free vvmax related things */
  free(Mabs);
  free(lfvvmax.magni);
  free(lfvvmax.errmagni);
  free(lfvvmax.lf);
  free(lfvvmax.errlf);
  free_matrix_d(lfvvmax.covarlf,lfvvmax.nbin,lfvvmax.nbin);

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


  _MLmax_STY_p_M=Amoe_Funk_STY_p_M_main(n,magn,z,par);


  /* Meto la solucion en la salida */
 
  lf->alfa=par[0];
  lf->Mstar=par[1];
  lf->phistar=exp(par[2]);

  /* Estimacion de los errores en los parametros */



  if(TRYEMPIRICAL) {
    if(DEBUG) printf(" llamantry \n");
    _conflim_STY_p_M=exp(-.5/10.);
    EmpiricalCovars_STY_p_M(n,magn,z,par,sigpar,mlim,cosmo,lf); 
    if(DEBUG) printf(" sale \n");

    if(DEBUG) printf(" Solucion final: Mstar %.15g +/- %.15g alpha %.4g +/- %.4g\n",lf->Mstar,lf->errMstar,lf->alfa,lf->erralfa);

  }

/*  printf("Antes del trynumerical\n"); */
  if(TRYNUMERICAL) {
    if(DEBUG) printf(" llamantry \n");
    _conflim_STY_p_M=exp(-.5/10.);
/*    printf("Dentro del trynumerical\n"); */
    NumericalDerivCovars_STY_p_M(n,magn,z,par,sigpar,mlim,cosmo,lf); 
    if(DEBUG) printf(" sale \n");

    if(DEBUG) printf(" Solucion final: Mstar %.15g +/- %.15g alpha %.4g +/- %.4g\n",lf->Mstar,lf->errMstar,lf->alfa,lf->erralfa);  
  
  }
/*  printf("Después del trynumerical\n"); */

  if(CONTOURPLOT) {
/*     nemp_f++; */
/*     conflim=exp(-.5*16.);     */
/*     EmpiricalCovars_STY_p_M(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*9.);     */
/*     EmpiricalCovars_STY_p_M(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*4.);     */
/*     EmpiricalCovars_STY_p_M(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/1.);     */
/*     EmpiricalCovars_STY_p_M(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/4.);     */
/*     EmpiricalCovars_STY_p_M(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     cpgsci(1); */
/*     cpglab("P\\d1\\u","P\\d3\\u","Contornos de límites de confianza"); */
  }
  if(DEBUG) printf(" Calculo empirico\n");

  
/*   conflim=readf(conflim); */

/*   ComputeNorma_STY_p_M(n, mlim,  strrad, zlow,  zup, lf); */
  
  free(y);
  if(DEBUG) printf(" MLF %g\n",_MLmax_STY_p_M); 

  if(iter_amo>=MAXITER-1) return(2);
  return(0);
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
/*   double Mup=-30; */
/*   double Lup; */
  double intsup;
  double Ntot;

  lfamo.alfa=p[0];
  lfamo.Mstar=p[1];
  lfamo.phistar=exp(p[2]);

  logL=0.;

  Lstar=pow(10.,-0.4*lfamo.Mstar);
/*  intsup=incom(1+lfamo.alfa,200.); */
  intsup=gsl_sf_gamma(1.+lfamo.alfa);

  /* dabreu */
/*  printf("   alfa %g incom %g  gamma %g  gamma_gsl %g\n",lfamo.alfa, intsup, exp(gammln(1+lfamo.alfa)), gsl_sf_gamma(1.+lfamo.alfa)); */

  for(i=0;i<_ndata_STY_p_M;i++)
  {
    Mabs=Mag(y[i],x[i],*_cosmo_STY_p_M);
    Mlow=Mag(y[i],_mlim_STY_p_M,*_cosmo_STY_p_M);
    Llow=pow(10.,-0.4*Mlow);

    /* intsch=lfamo.phistar*(intsup-incom(1+lfamo.alfa,Llow/Lstar)); */
    /* debido a un underflow, tuvimos que poner este if */
    if(Llow/Lstar > 0.25 && (lfamo.alfa*log(Llow/Lstar) - Llow/Lstar) <= GSL_LOG_DBL_MIN)
    {
      log_gamma_int=GSL_LOG_DBL_MIN;
    }
    else
    {
      log_gamma_int=log(gsl_sf_gamma_inc(1+lfamo.alfa,Llow/Lstar));
    }
    logL-= log(Schechter_M(Mabs,lfamo)) - log(lfamo.phistar) - log_gamma_int;
    
/* dabreu    printf(" ndata %d  logL %g x %f  y %f Mabs %f Lstar %g Llow %g intsup %g  sch %g int %g\n",i,logL,x[i],y[i],Mabs,Lstar,Llow,intsup,Schechter_M(Mabs,lfamo),intsch); */
/*     printf(" LF: lfamo.Mstar %g lfamo.phistar %f lf.alfa  %f\n",lfamo.Mstar,lfamo.phistar,lfamo.alfa); */
  }

  if(DEBUG) printf(" logL %g\n",logL);
  /* Aqui viene la parte de la poissoniana de npob */
  Ntot=Int_sch_M(lfamo,_zlow_STY_p_M,_zup_STY_p_M,_mlim_STY_p_M,*_cosmo_STY_p_M)*_strrad_STY_p_M/4./M_PI;
  logL-= (_ndata_STY_p_M*log(Ntot) - Ntot - gammln((double)_ndata_STY_p_M+1.));

  if(DEBUG) printf(" Ntot %g ndata %d rad %g zlw %f zup %f mag %f lfamophi %f\n",Ntot, _ndata_STY_p_M, _strrad_STY_p_M,_zlow_STY_p_M,_zup_STY_p_M,_mlim_STY_p_M,lfamo.phistar);

  _iter_m_STY_p_M++;
  if(DEBUG) printf(" Iter %d  logL %f par0 %g par1 %g par2 %g (%.10g)\n",_iter_m_STY_p_M,logL,p[0],p[1],exp(p[2]),p[2]);
  return(logL);
}

double Amoe_Funk_STY_p_M_conf(int n, double *x, double *y, double *p) 
{
/*   printf(" MLM %g conf %f   %f %f\n",MLmax,conflim,fabs(Amoe_Funk_STY_p_M_main(n,x,y,p)-(MLmax-log(conflim))),Amoe_Funk_STY_p_M_main(n,x,y,p)); */
   
  if(DEBUG) printf(" AmoeFunkConf conf %f return = %g\n",_conflim_STY_p_M,Amoe_Funk_STY_p_M_main(n,x,y,p)-(_MLmax_STY_p_M-log(_conflim_STY_p_M)));
  return(fabs(Amoe_Funk_STY_p_M_main(n,x,y,p)-(_MLmax_STY_p_M-log(_conflim_STY_p_M))));
}



void EmpiricalCovars_STY_p_M(int n,double *magn,double *z,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf)
{

  int i,j;  
  double *parconf; 
  double *sigparconf; 
  double **invcovar;
  double **covar;
  double **parelip;
  double *y;

  double **bb;
  int nconfl,nconflini;
  double first, median, third, *distmax;

  if(DEBUG) printf(" n vale %d \n",n);
  nconfl=NCONFL;
  nconflini=NCONFL;

  bb=matrix_d(3,1);
  y=vector_d(n);
  parconf=vector_d(3);
  sigparconf=vector_d(3);
  parelip=matrix_d(3,nconflini);
  invcovar=matrix_d(3 ,3 );
  covar=matrix_d(3 ,3 );
  distmax=vector_d(3);
 
/*   printf(" antes for \n"); */

  nconfl=NCONFL;
  printf(" Error step  ");
  for(i=0;i<nconfl;i++) printf(".");
  for(i=0;i<nconfl;i++) printf("\b");
  fflush(NULL);

  for(i=0;i<nconfl;i++) 
  {
    printf("#");
    fflush(NULL);
    parconf[0]=par[0]+3*sigpar[0]*Gasdev();
    parconf[1]=par[1]+3*sigpar[1]*Gasdev();
    parconf[2]=par[2]+3*sigpar[2]*Gasdev();
    sigparconf[0]=sigpar[0];
    sigparconf[1]=sigpar[1];
    sigparconf[2]=sigpar[2];
    if(i>(int)(nconfl/2.)) 
    {
      parconf[0]=par[0]-((parelip[0])[(int)(i-nconfl/2.)+1]-par[0]);
      parconf[1]=par[1]-((parelip[1])[(int)(i-nconfl/2.)+1]-par[1]);
      parconf[2]=par[2]-((parelip[2])[(int)(i-nconfl/2.)+1]-par[2]);
    }
    _iter_c_STY_p_M=0;
    if(DEBUG) printf(" Calling amoeba for error step %d\n",i+1); 
    _iter_c_STY_p_M=Amoeba_d(n,magn,z,3,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_STY_p_M_conf);
    if(VERBOSE) printf("CONF %d SOL iter %d par0 %f    par1 %f  par2 %f\n", i,  _iter_c_STY_p_M ,parconf[0],   parconf[1], exp(parconf[2]));
    (parelip[0])[i]=parconf[0];
    (parelip[1])[i]=parconf[1];
    (parelip[2])[i]=parconf[2];
    if(_iter_c_STY_p_M==0 && Amoe_Funk_STY_p_M_conf(n,magn,z,parconf)>FTOL2 ) 
    {
      printf("\b.\b");
      fflush(NULL);
      i--;
    }
  }

  /* Supongo que el centro de la elipse es el valor que maximiza ML */
  for(i=0;i<nconfl;i++) 
  {
    (parelip[0])[i]-=par[0];
    (parelip[1])[i]-=par[1];
    (parelip[2])[i]-=par[2];
  }

  /* Detecto puntos que esten muy alejados de la supuesta elipse */
  for(j=0;j<3;j++) 
  {
    Quartil_d(nconfl,parelip[j],&first,&median,&third);
    distmax[j]=maxf(fabs(first),fabs(third));
  }  
  for(i=0;i<nconfl;i++) 
  {
    for(j=0;j<3;j++) 
    {
      if(fabs((parelip[j])[i])>2*2*distmax[j]/1.35) 
      {
	for(j=0;j<3;j++) 
	{
	  memcpy(parelip[j]+i,parelip[j]+i+1,(nconfl-i-1)*sizeof(double));
	}
	i--;
	nconfl--;
	break;
      }
    }
  }

  MCElipN_d(nconfl,3,parelip,invcovar);
  if(VERBOSE) printf(" invcovar(0,0) %g invcovar(0,1) %g invcovar(1,0) %g invcovar(1,1) %g\n",invcovar[0][0],invcovar[0][1],invcovar[1][0],invcovar[1][1]);
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      covar[i][j]=invcovar[i][j];
    }
  }
  gaussj_d(covar,3,bb,1);
  if(VERBOSE) printf("    covar(0,0) %g    covar(0,1) %g    covar(1,0) %g    covar(1,1) %g\n",   covar[0][0],   covar[0][1],   covar[1][0],   covar[1][1]);
  /* Deshago el cambio para el limite de confidencia */
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      covar[i][j]/=(-2*log(_conflim_STY_p_M));
    }
  }
  if(VERBOSE) printf("CON covar(0,0) %g    covar(0,1) %g    covar(1,0) %g    covar(1,1) %g\n",   covar[0][0],   covar[0][1],   covar[1][0],   covar[1][1]);

  lf->erralfa=sqrt(covar[0][0]);
  lf->errMstar=sqrt(covar[1][1]);
  lf->errphistar=sqrt(covar[2][2])*lf->phistar; /* Ya que p[2] esta en logaritmos */
  lf->covaralfaMstar=covar[0][1];
  lf->covaralfaphistar=covar[0][2]*lf->phistar;  /* Por la misma razon */
  lf->covarphistarMstar=covar[1][2]*lf->phistar;

  free_matrix_d(bb,3,1);
  free(y);
  free(parconf);
  free(sigparconf);

  free_matrix_d(parelip,3,nconflini);
  free_matrix_d(invcovar,3  ,3);
  free_matrix_d(   covar,3  ,3);

  printf("\n");
}


void   ComputeNorma_STY_p_M(int n, double mlim, double strrad, double zlow, double zup, struct Schlf_M *lf) {

  double Ntot[NBSNORMA];
  double pi=3.1415926535897932384;
  struct Schlf_M lfbs[NBSNORMA];
  double Ntot_mean, Ntot_sigma;
  int i;

  for(i=0;i<NBSNORMA;i++) {
    lfbs[i].Mstar=lf->Mstar+lf->errMstar*Gasdev();
    lfbs[i].alfa =lf->alfa +lf->erralfa *Gasdev();
    lfbs[i].phistar=1;
    Ntot[i]=Int_sch_M(lfbs[i],zlow,zup,mlim,*_cosmo_STY_p_M)*strrad/4./pi;
    if(VERBOSE) {
      printf(" %d Ntot %g\n",i,Ntot[i]);
      printf(" Mstar %f  (Mstar %f err %f)\n",lfbs[i].Mstar,lf->Mstar,lf->errMstar);
      printf(" alfa %f  (alfa  %f err %f)\n",lfbs[i].alfa,lf->alfa,lf->erralfa);
    }
  }
  Ntot_mean=StMedia_d(NBSNORMA, Ntot, &Ntot_sigma);

  lf->phistar=1;
  Ntot_mean=Int_sch_M(*lf,zlow,zup,mlim,*_cosmo_STY_p_M)*strrad/4./pi;

  if(VERBOSE) {
    printf(" Ntot_mean %g Ntot_sigma %g\n",Ntot_mean, Ntot_sigma);
    printf(" phistar_mean %g phistar_sigma %g\n",n/Ntot_mean, n*Ntot_sigma/Ntot_mean/Ntot_mean);
  }

  lf->phistar=(float)(n)/Ntot_mean;
  /* Sumo el error poissoniano del n de la muestra mas el error derivado 
     del calculo de Ntot debido a los errores en Mstar y alfa: Ntot_sigma */
  lf->errphistar=sqrt( sqrt(n)/Ntot_mean*sqrt(n)/Ntot_mean + (Ntot_sigma*n/Ntot_mean/Ntot_mean)*(Ntot_sigma*n/Ntot_mean/Ntot_mean) );
  lf->errphistar=sqrt(n)/Ntot_mean;
  lf->errphistar=Ntot_sigma*n/Ntot_mean/Ntot_mean;
  if(VERBOSE) printf(" Error poissoniano: %g \n",sqrt(n)/Ntot_mean);


}

/* Errores utilizando derivadas numéricas con GSL */

void NumericalDerivCovars_STY_p_M(int n,double *magn,double *z,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf)
{
  double **hessian;
  double **covar;
  double **bb;
  size_t i,j;

  double result, abserr, derivStep;
  derivStep = GSL_ROOT4_DBL_EPSILON;
  /* derivStep = 0.01; */

  hessian=matrix_d(3 ,3 );
  covar=matrix_d(3 ,3 );
  bb=matrix_d(3,1);

  struct num_deriv_param_STY_p_M deriv_param;
  
  deriv_param.lf_param = vector_d(3);
  
  deriv_param.nData = _ndata_STY_p_M;
  deriv_param.magn = magn;
  deriv_param.z = z;
  deriv_param.lf_param[0] = par[0];
  deriv_param.lf_param[1] = par[1];
  deriv_param.lf_param[2] = par[2];

  gsl_function derivLogLFunction;
  derivLogLFunction.function = &derivAmoe_Funk_STY_p_M_main_gsl;
  derivLogLFunction.params = &deriv_param;
  
  for(i = 0 ; i < 3 ; ++i)
    for(j = 0 ; j <= i ; ++j)
    {
      deriv_param.iCovar = i;
      deriv_param.jCovar  = j;
      /* x_i = par[i]; */
      gsl_deriv_central(&derivLogLFunction, par[i], derivStep, &result, &abserr);
/*      printf(" Con iCovar %d jCovar %d derivResult %g abserr %g\n",i,j,result,abserr); */
      hessian[i][j] = result;
      hessian[j][i] = result;
    }

  for(i=0;i<3;++i) 
  {
    for(j=0;j<3;++j) 
    {
      covar[i][j]=hessian[i][j];
    }
  }
  gaussj_d(covar,3,bb,1);
   

  lf->erralfa=sqrt(covar[0][0]);
  lf->errMstar=sqrt(covar[1][1]);
  lf->errphistar=sqrt(covar[2][2])*lf->phistar; /* Ya que p[2] esta en logaritmos */
  lf->covaralfaMstar=covar[0][1];
  lf->covaralfaphistar=covar[0][2]*lf->phistar;  /* Por la misma razon */
  lf->covarphistarMstar=covar[1][2]*lf->phistar;

  free_matrix_d(bb,3,1);
  free(deriv_param.lf_param);
}


double derivAmoe_Funk_STY_p_M_main_gsl(double x, void *params)
{
  int iCovar              = ((struct num_deriv_param_STY_p_M *)params)->iCovar;
  int jCovar              = ((struct num_deriv_param_STY_p_M *)params)->jCovar;
  double* lf_param        = ((struct num_deriv_param_STY_p_M *)params)->lf_param;

  double derivResult, derivAbserr, derivStep;
  derivStep = GSL_SQRT_DBL_EPSILON;
  /* derivStep = 0.01; */

  gsl_function logLfunction;
  logLfunction.function = &Amoe_Funk_STY_p_M_main_gsl;
  logLfunction.params = params;
  
  ((struct num_deriv_param_STY_p_M *)params)->iDeriv = jCovar;
  
  /* x = lf_param[jCovar] */
  lf_param[iCovar] = x;
  
  gsl_deriv_central(&logLfunction, lf_param[jCovar], derivStep, &derivResult, &derivAbserr);

/*  printf(" Con p0 %g p1 %g p2 %g  derivResult %g abserr
%g\n",lf_param[0],lf_param[1],lf_param[2],derivResult,derivAbserr); */

  return derivResult;
}

double Amoe_Funk_STY_p_M_main_gsl(double x, void *params)
{
  int nData           = ((struct num_deriv_param_STY_p_M *)params)->nData;
  int iDeriv          = ((struct num_deriv_param_STY_p_M *)params)->iDeriv;
  double* magn        = ((struct num_deriv_param_STY_p_M *)params)->magn;
  double* z           = ((struct num_deriv_param_STY_p_M *)params)->z;
  double* p           = ((struct num_deriv_param_STY_p_M *)params)->lf_param;
  double value;
  
  double* lf_param = vector_d(3);
  lf_param[0] = p[0];
  lf_param[1] = p[1];
  lf_param[2] = p[2];
  
  lf_param[iDeriv] = x;

  value=Amoe_Funk_STY_p_M_main(nData, magn, z, lf_param);
  /*  value = probando_deriv(nData, magn, z, lf_param);  */
  free(lf_param);
  return value;
}

/* Función para comprobar que la derivada se está haciendo correcatamente */
double probando_deriv(int n, double *x, double *y, double *p)
{
  double result;
  result =  4*pow(p[0],3)+5*pow(p[1],3)+6*pow(p[2],3)+3*pow(p[0],2)*pow(p[1],2)+6*pow(p[1],2)*pow(p[2],2)+9*pow(p[0],2)*pow(p[2],2);
  /*  printf(" p0 %g p1 %g p2 %g result %g\n",p[0],p[1],p[2],result); */
  return result;
}
