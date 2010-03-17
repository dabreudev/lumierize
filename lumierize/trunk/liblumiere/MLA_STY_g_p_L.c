#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_machine.h>
#include "alloc.h"
#include "cosmology.h"
#include "schechter.h"
#include "mlsty.h"
#include "amoeba.h"
#include "step.h"
#include "functions.h"
#include "random.h"
#include "gaussj.h"
#include "minmax.h"
#include "elip.h"
#include "quartil.h"
#include "stmedia.h"
#include "vvmax.h"

#define FTOL  1e-13
#define FTOL2 1e-6
#define FTOL3 1e-7
#define MAXITER  1000
#define MAXITER2 120
#define NBSNORMA  50
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
#define TRYEMPIRICAL 1
#endif

/* #define TOLERR 0.0001 */

/* Ahora mismo esta con el pcut definido, de modo que se aleja de valores 
   cercanos a 0. Puede afectar al calculo de las covarianzas!! */

double Amoe_Funk_STY_g_p_L_main(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_g_p_L_norm(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_g_p_L_prior(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_g_p_L_conf(int n, double *x, double *y, double *p);
double Funk2_int_STY_g_p_L(double x);
double Funk1_norm_STY_g_p_L(double x);
void   EmpiricalCovars_STY_g_p_L(int n,double *flux,double *errflux, double *z,double *errz, double *par, double *sigpar,double flim, struct cosmo_param cosmo,struct Schlf_L *lf);

void   ComputeNorma_STY_g_p_L(int n, double flim, double strrad, double zlow, double zup, struct Schlf_L *lf);
double *fluxsig;
double *zsig;
struct cosmo_param *co;
double fluxl;

int ndata;
int iter_m;
int iter_c;
int nconfl_STY_g_p_L;
double conflim;
/* double ML[5*MAXTRIES*MAXITER]; */
double *pp;
double MLmax;
double *pcut;
double xf,Tf;
double sigi;
double xtmp;
double zlow_g,zup_g,strrad_g;


int  MLA_STY_g_p_L(int n,double *flux,double *errflux,double *z,double *errz,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf) {

  double *y;
  double par[3];
  double sigpar[3];
  int i;
  int iter_amo;

  struct Steplf_L  lfvvmax;
  double minlum,maxlum;
  double *lumi;
  double chisq;
  struct Schlf_L  lffit;

  lumi=vector_d(n);

  /* Copias globales de estas variables para que se vean en toda la subrutina */
  zsig=errz;
  fluxsig=errflux;
  ndata=n;
  co=&cosmo;
  fluxl=flim;
  strrad_g=strrad;
  zlow_g=zlow;
  zup_g=zup;

  y=vector_d(n);
  iter_m=0;

  if(DEBUG3) {
    for(i=0;i<n;i++) printf(" Entrada x %g\n",flux[i]);
  }

  lfvvmax.nbin=(int)(n/3);
  if(lfvvmax.nbin>30) lfvvmax.nbin=30;

  lfvvmax.lumi      =vector_d(lfvvmax.nbin+1);
  lfvvmax.errlumi   =vector_d(lfvvmax.nbin+1);
  lfvvmax.lnlf      =vector_d(lfvvmax.nbin);
  lfvvmax.errlnlf   =vector_d(lfvvmax.nbin);
  lfvvmax.covarlnlf =matrix_d(lfvvmax.nbin,lfvvmax.nbin);

  for(i=0;i<n;i++)   lumi[i]=log(Lum(z[i],flux[i],cosmo));
  MinMax_d(n,lumi,&minlum,&maxlum);
  for(i=0;i<=lfvvmax.nbin;i++) lfvvmax.lumi[i]=minlum+i*(maxlum-minlum)/lfvvmax.nbin;
  VVmax_L(n,flux,flux,z,flim,strrad,zlow,zup,cosmo,&lfvvmax);
/*   cpgopen("?"); */
/*   PlotStepLF_L(lfvvmax); */
  if(DEBUG) printf(" Salida VVmax\n");
  if(DEBUG) for(i=0;i<lfvvmax.nbin;i++) printf(" Lum %g - %g LF %g\n",lfvvmax.lumi[i]/log(10),lfvvmax.lumi[i+1]/log(10),lfvvmax.lnlf[i]/log(10));
  FitSch2StepLF_L(lfvvmax,&lffit, &chisq);
  if(DEBUG) {
    printf(" Desuess ajuste MRQ\n");
    printf(" Schechter FIT\n");
    printf(" log(Lstar (W)):  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",log10(lffit.Lstar),lffit.alfa,lffit.phistar,log10(lffit.phistar));
    printf(" E_log(Lstar):    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lffit.errLstar/lffit.Lstar/log(10.),lffit.erralfa,lffit.errphistar,lffit.errphistar/lffit.phistar/log(10.));
    printf(" Covar(Lstar,alpha): %g Covar(log(Lstar),alpha): %g\n",lffit.covaralfaLstar,lffit.covaralfaLstar/lffit.Lstar/log(10.));
    printf(" Covar(Lstar,Phistar): %g Covar(log(Lstar),log(Phistar)): %g\n",lffit.covarphistarLstar,lffit.covarphistarLstar/lffit.Lstar/log(10.)/lffit.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lffit.covaralfaphistar,lffit.covaralfaphistar/lffit.phistar/log(10.));
  }

/*   PlotStepSchLF_L(lfvvmax,lffit); */

/*   cpgclos(); */
  
  
  printf(" Computing LF...\n");

  iter_amo=0;
  while(iter_amo==0) { 
    par[0]=lffit.alfa;  /* Alpha */
    par[1]=log10(lffit.Lstar);
    par[2]=log(lffit.phistar);
    sigpar[0]=5.*lffit.erralfa;
    sigpar[1]=3.*lffit.errLstar/lffit.Lstar/log(10.);
    sigpar[2]=3.*lffit.errphistar/lffit.phistar;
    iter_amo=Amoeba_d(n,flux,z,3,par,sigpar,FTOL,MAXITER,Amoe_Funk_STY_g_p_L_main);  
    if(DEBUG) printf(" iteramo %d\n",iter_amo);
    lf->alfa=par[0];
    lf->Lstar=pow(10.,par[1]);
    lf->phistar=exp(par[2]);
    if(DEBUG) {
      printf(" Solucion MALA\n");
      printf(" log(Lstar (W)):  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",log10(lf->Lstar),lf->alfa,lf->phistar,log10(lf->phistar));
      printf(" E_log(Lstar):    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lf->errLstar/lf->Lstar/log(10.),lf->erralfa,lf->errphistar,lf->errphistar/lf->phistar/log(10.));
    } 
  }
  MLmax=Amoe_Funk_STY_g_p_L_main(n,flux,z,par);


  /* Meto la solucion en la salida */
 
  lf->alfa=par[0];
  lf->Lstar=pow(10.,par[1]);
  lf->phistar=exp(par[2]);

  /* Estimacion de los errores en los parametros */



  if(TRYEMPIRICAL) {
    if(DEBUG) printf(" llamantry \n");
    conflim=exp(-.5/10.);
    EmpiricalCovars_STY_g_p_L(n,flux,errflux,z,errz,par,sigpar,flim,cosmo,lf); 
    if(DEBUG) printf(" sale \n");

    if(DEBUG) printf(" Solucion final: Lstar %.15g +/- %.15g alpha %.4g +/- %.4g\n",lf->Lstar,lf->errLstar,lf->alfa,lf->erralfa);

  }
  if(CONTOURPLOT) {
/*     nemp_f++; */
/*     conflim=exp(-.5*16.);     */
/*     EmpiricalCovars_STY_g_p_L(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*9.);     */
/*     EmpiricalCovars_STY_g_p_L(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*4.);     */
/*     EmpiricalCovars_STY_g_p_L(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/1.);     */
/*     EmpiricalCovars_STY_g_p_L(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/4.);     */
/*     EmpiricalCovars_STY_g_p_L(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     cpgsci(1); */
/*     cpglab("P\\d1\\u","P\\d3\\u","Contornos de l�mites de confianza"); */
  }
  if(DEBUG) printf(" Calculo empirico\n");

  
/*   conflim=readf(conflim); */

/*   ComputeNorma_STY_g_p_L(n, flim,  strrad, zlow,  zup, lf); */
  
  free(y);
  free(lumi);
  free(lfvvmax.lumi);
  free(lfvvmax.errlumi);
  free(lfvvmax.lnlf);
  free(lfvvmax.errlnlf);
  free_matrix_d(lfvvmax.covarlnlf,lfvvmax.nbin,lfvvmax.nbin);
  return(iter_amo);
}


double Amoe_Funk_STY_g_p_L_main(int n, double *x, double *y, double *p) {

  int i; 
  double logL=0.;
  struct Schlf_L lfamo;
  struct Schlf_M lfamo_M;
  double Lumi;
  double Mabs;
  double intsch;
  double Llow;
  double xmin;
/*   double Mup=-30; */
/*   double Lup; */
  double intsup;
  double Ntot;
  (void)n;


  lfamo.alfa=p[0];
  lfamo.Lstar=pow(10.,p[1]);
  lfamo.phistar=exp(p[2]);

  lfamo_M.alfa=p[0];
  lfamo_M.Mstar=-2.5*p[1];
  lfamo_M.phistar=exp(p[2]);

  logL=0.;

  intsup=incom(1+lfamo.alfa,200.);
  for(i=0;i<ndata;i++) {
    Lumi=Lum(y[i],x[i],*co);
    Mabs=-2.5*log10(Lumi);
    Llow=Lum(y[i],fluxl,*co);
    xmin=Llow/lfamo.Lstar;
    if(xmin>=50) {
      logL+=GSL_LOG_DBL_MAX; 
    }
    else {
      intsch=lfamo.phistar*(intsup-incom(1+lfamo.alfa,xmin));
      /* Asi es haciendolo con luminosidades 
	 logL-= Schechter_L(Lumi,lfamo) -log(intsch); 
	 Pero prefiero usar la funcion de Sch en magni */
      logL-= log(Schechter_M(Mabs,lfamo_M)/intsch);
      /* Funciona perfectamente en magnitudes */
    }
/*     printf(" ndata %d  logL %g x %g  y %g Lumi %g Mabs %f alfa %f  Lstar %g Mstar %f Llow %g intsup %g  sch %g int %g xmin %g inco %g sch-int%g\n",i,logL,x[i],y[i],Lumi,Mabs,lfamo.alfa,lfamo.Lstar,lfamo_M.Mstar,Llow,intsup,Schechter_L(Lumi,lfamo),intsch,Llow/lfamo.Lstar,incom(1+lfamo.alfa,Llow/lfamo.Lstar),Schechter_L(Lumi,lfamo) -log(intsch)); */
/*     printf(" LF: lfamo.Lstar %g lfamo.phistar %f lf.alfa  %f\n",lfamo.Lstar,lfamo.phistar,lfamo.alfa); */
  }
  /* Aqui viene la parte de la poissoniana de npob */
  Ntot=Int_sch_L(lfamo,zlow_g,zup_g,fluxl,*co)*strrad_g/4./M_PI;
  logL-= (ndata*log(Ntot) - Ntot - gammln((double)ndata+1.));
  
  if(DEBUG) printf(" NTOT %f ndata*log(Ntot) %f gamm %f\n",Ntot,ndata*log(Ntot),gammln((double)ndata)+1.);

  iter_m++;
  if(DEBUG) printf(" Iter %d  logL %f par0 %g par1 %g\n",iter_m,logL,p[0],p[1]);
  return(logL);
}


double Amoe_Funk_STY_g_p_L_conf(int n, double *x, double *y, double *p) {

/*   printf(" MLM %g conf %f   %f %f\n",MLmax,conflim,fabs(Amoe_Funk_STY_g_p_L_main(n,x,y,p)-(MLmax-log(conflim))),Amoe_Funk_STY_g_p_L_main(n,x,y,p)); */
   
  return(fabs(Amoe_Funk_STY_g_p_L_main(n,x,y,p)-(MLmax-log(conflim)))); 
}



void   EmpiricalCovars_STY_g_p_L(int n,double *flux,double *errflux,double *z,double *errz,double *par, double *sigpar,double flim, struct cosmo_param cosmo,struct Schlf_L *lf) {

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
  (void)cosmo;
  (void)errz;
  (void)errflux;
  (void)flim;



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

  for(i=0;i<nconfl;i++) {
    printf("#");
    fflush(NULL);
    parconf[0]=par[0]+3*sigpar[0]*Gasdev();
    parconf[1]=par[1]+3*sigpar[1]*Gasdev();
    parconf[2]=par[2]+3*sigpar[2]*Gasdev();
    sigparconf[0]=sigpar[0];
    sigparconf[1]=sigpar[1];
    sigparconf[2]=sigpar[2];
    if(i>(int)(nconfl/2.)) {
      parconf[0]=par[0]-((parelip[0])[(int)(i-nconfl/2.)+1]-par[0]);
      parconf[1]=par[1]-((parelip[1])[(int)(i-nconfl/2.)+1]-par[1]);
      parconf[2]=par[2]-((parelip[2])[(int)(i-nconfl/2.)+1]-par[2]);
    }
    iter_c=0;
    if(DEBUG) printf(" antes a nmo\n"); 
    iter_c=Amoeba_d(n,flux,z,3,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_STY_g_p_L_conf);
    if(DEBUG) printf(" %d SOL iter %d par0 %f    par1 %f\n", i,  iter_c ,parconf[0],   parconf[1]);
    (parelip[0])[i]=parconf[0];
    (parelip[1])[i]=parconf[1];
    (parelip[2])[i]=parconf[2];
    if(iter_c==0 && Amoe_Funk_STY_g_p_L_conf(n,flux,z,parconf)>FTOL2 ) {
      printf("\b.\b");
      fflush(NULL);
      i--;
    }
  }

  /* Supongo que el centro de la elipse es el valor que maximiza ML */
  for(i=0;i<nconfl;i++) {
    (parelip[0])[i]-=par[0];
    (parelip[1])[i]-=par[1];
    (parelip[2])[i]-=par[2];
  }

  /* Detecto puntos que esten muy alejados de la supuesta elipse */
  for(j=0;j<3;j++) {
    Quartil_d(nconfl,parelip[j],&first,&median,&third);
    distmax[j]=maxf(fabs(first),fabs(third));
  }  
  for(i=0;i<nconfl;i++) {
    for(j=0;j<3;j++) {
      if(fabs((parelip[j])[i])>2*2*distmax[j]/1.35) {
	for(j=0;j<3;j++) {
	  memcpy(parelip[j]+i,parelip[j]+i+1,(nconfl-i-1)*sizeof(double));
	}
	i--;
	nconfl--;
	break;
      }
    }
  }

  MCElipN_d(nconfl,3,parelip,invcovar);
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      covar[i][j]=invcovar[i][j];
    }
  }
  gaussj_d(covar,3,bb,1);
  /* Deshago el cambio para el limite de confidencia */
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      covar[i][j]/=(-2*log(conflim));
    }
  }

  lf->erralfa=sqrt(covar[0][0]);
  lf->errLstar=sqrt(covar[1][1])*lf->Lstar*log(10.);
  lf->errphistar=sqrt(covar[2][2])*lf->phistar; /* Ya que p[2] esta en logaritmos */
  lf->covaralfaLstar=covar[0][1]*lf->Lstar*log(10.);
  lf->covaralfaphistar=covar[0][2]*lf->phistar;  /* Por la misma razon */
  lf->covarphistarLstar=covar[1][2]*lf->phistar*lf->Lstar*log(10.);

  free_matrix_d(bb,3,1);
  free(y);
  free(parconf);
  free(sigparconf);

  free_matrix_d(parelip,3,nconflini);
  free_matrix_d(invcovar,3  ,3);
  free_matrix_d(   covar,3  ,3);

  printf("\n");
}





void   ComputeNorma_STY_g_p_L(int n, double flim, double strrad, double zlow, double zup, struct Schlf_L *lf) {

  double Ntot[NBSNORMA];
  double pi=3.1415926535897932384;
  struct Schlf_L lfbs[NBSNORMA];
  double Ntot_mean, Ntot_sigma;
  int i;

  for(i=0;i<NBSNORMA;i++) {
    lfbs[i].Lstar=lf->Lstar+lf->errLstar*Gasdev();
    lfbs[i].alfa =lf->alfa +lf->erralfa *Gasdev();
    lfbs[i].phistar=1;
    Ntot[i]=Int_sch_L(lfbs[i],zlow,zup,flim,*co)*strrad/4./pi;
    if(VERBOSE) {
      printf(" %d Ntot %g\n",i,Ntot[i]);
      printf(" Lstar %g  (Lstar %g err %g)\n",lfbs[i].Lstar,lf->Lstar,lf->errLstar);
      printf(" alfa %f  (alfa  %f err %f)\n",lfbs[i].alfa,lf->alfa,lf->erralfa);
    }
  }
  Ntot_mean=StMedia_d(NBSNORMA, Ntot, &Ntot_sigma);

  lf->phistar=1;
  Ntot_mean=Int_sch_L(*lf,zlow,zup,flim,*co)*strrad/4./pi;

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
