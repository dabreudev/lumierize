#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_machine.h>
#include "alloc.h"
#include "cosmology.h"
#include "schechter.h"
#include "mlsty.h"
#include "vvmax.h"
#include "minmax.h"
#include "amoeba.h"
#include "functions.h"
#include "random.h"
#include "quartil.h"


#define FTOL  1e-11
#define FTOL2 1e-6
#define FTOL3 1e-7
#define MAXITER  1000
#define MAXITER2 120
#define NBSNORMA  50
#define MAXTRIES   5
#define VERBOSE 0
#define DEBUG  1
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


double Amoe_Funk_STY_p_f_L_main(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_p_f_L_norm(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_p_f_L_prior(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_p_f_L_conf(int n, double *x, double *y, double *p);
double Funk2_int_STY_p_f_L(double fluxreal);
double Funk1_intlumlog_STY_p_f_L(double x);
void   EmpiricalCovars_STY_p_f_L(int n,double *flux, double *z, double *par, double *sigpar,struct fermifsel_L fsel, struct cosmo_param cosmo,struct Schlf_L *lf);

void   ComputeNorma_STY_p_f_L(int n, struct fermifsel_L fsel, double strrad, double zlow, double zup, struct Schlf_L *lf);

struct Schlf_L *lf_STY_p_f_L;
struct cosmo_param *co_STY_p_f_L;
double fluxl_STY_p_f_L;
double z_STY_p_f_L;
double flux_STY_p_f_L;
struct fermifsel_L fsel_STY_p_f_L;

int ndata;
int iter_m;
int iter_c;
int nconfl_p_f_L;
double conflim;
/* double ML[5*MAXTRIES*MAXITER]; */
double *pp;
double MLmax;
double zlow_STY_p_f_L,zup_STY_p_f_L,strrad_STY_p_f_L;


int  MLA_STY_p_f_L(int n,double *flux,double *z,struct fermifsel_L fsel, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_L *lf) {

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
  
  ndata=n;
  co_STY_p_f_L=&cosmo;
  fsel_STY_p_f_L=fsel;
  strrad_STY_p_f_L=strrad;
  zlow_STY_p_f_L=zlow;
  zup_STY_p_f_L=zup;

  y=vector_d(n);
  iter_m=0;

  if(DEBUG3) {
    for(i=0;i<n;i++) printf(" Entrada x %g\n",flux[i]);
  }

  lfvvmax.nbin=(int)(n/3);
  if(lfvvmax.nbin>30) lfvvmax.nbin=30;

  lfvvmax.lumi       =vector_d(lfvvmax.nbin+1);
  lfvvmax.errlumi    =vector_d(lfvvmax.nbin+1);
  lfvvmax.lnlf       =vector_d(lfvvmax.nbin);
  lfvvmax.errlnlf    =vector_d(lfvvmax.nbin);
  lfvvmax.covarlnlf  =matrix_d(lfvvmax.nbin,lfvvmax.nbin);

  for(i=0;i<n;i++)   lumi[i]=log(Lum(z[i],flux[i],cosmo));
  MinMax_d(n,lumi,&minlum,&maxlum);
  for(i=0;i<=lfvvmax.nbin;i++) lfvvmax.lumi[i]=minlum+i*(maxlum-minlum)/lfvvmax.nbin;
  VVmax_L(n,flux,flux,z,fsel_STY_p_f_L.fluxcut,strrad,zlow,zup,cosmo,&lfvvmax);
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
/*     par[0]=-0.9; */
    par[1]=log10(lffit.Lstar);
    par[2]=log(lffit.phistar);
    sigpar[0]=6.*lffit.erralfa;
    sigpar[1]=6.*lffit.errLstar/lffit.Lstar/log(10.);
    sigpar[2]=6.*lffit.errphistar/lffit.phistar;
    iter_amo=Amoeba_d(n,flux,z,3,par,sigpar,FTOL,MAXITER,Amoe_Funk_STY_p_f_L_main);  
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
  MLmax=Amoe_Funk_STY_p_f_L_main(n,flux,z,par);


  /* Meto la solucion en la salida */
 
  lf->alfa=par[0];
  lf->Lstar=pow(10.,par[1]);
  lf->phistar=exp(par[2]);

  /* Estimacion de los errores en los parametros */



  if(TRYEMPIRICAL) {
    if(DEBUG) printf(" llamantry \n");
    conflim=exp(-.5/10.);
    EmpiricalCovars_STY_p_f_L(n,flux,z,par,sigpar,fsel,cosmo,lf); 
    if(DEBUG) printf(" sale \n");

    if(DEBUG) printf(" Solucion final: Lstar %.15g +/- %.15g alpha %.4g +/- %.4g\n",lf->Lstar,lf->errLstar,lf->alfa,lf->erralfa);

  }

/*   ComputeNorma_STY_p_f_L(n, flim,  strrad, zlow,  zup, lf); */
  
  free(y);
  free(lumi);
  free(lfvvmax.lumi);
  free(lfvvmax.errlumi);
  free(lfvvmax.lnlf);
  free(lfvvmax.errlnlf);
  free_matrix_d(lfvvmax.covarlnlf,lfvvmax.nbin,lfvvmax.nbin);
  return(iter_amo);
}


double Amoe_Funk_STY_p_f_L_main(int n, double *x, double *y, double *p) {

  int i; 
  double logL=0.;
  double logLold=0.;
  struct Schlf_L lfamo;
  double Lumi;
  double L;
  double fluxtmp,fluxmin;
  double Lleft,Lright;
  double xleft,xright;
/*   double Mup=-30; */
/*   double Lup; */
  double Ntot;
  double intabajo;
  double logprobarriba;
  double probabajo;

  double intsup;
  int j;
  int nL_fs=200;

  lfamo.alfa=p[0];
  lfamo.Lstar=pow(10.,p[1]);
  lfamo.phistar=exp(p[2]);

  lf_STY_p_f_L=&lfamo;

  logL=0.;
  logLold=0.;

  if(DEBUG3) printf(" Entra con  par0 %g par1 %g par2 %g\n",p[0],p[1],p[2]);

  intsup=incom(1+lfamo.alfa,100.);
  for(i=0;i<ndata;i++) {
    z_STY_p_f_L=y[i];
    flux_STY_p_f_L=x[i];
    Lumi=Lum(y[i],x[i],*co_STY_p_f_L);
    fluxmin=fsel_STY_p_f_L.fluxcut-5*fsel_STY_p_f_L.deltaflux;
    if(fluxmin<=0) fluxmin=fsel_STY_p_f_L.fluxcut/20;
    Lleft=Lum(z_STY_p_f_L,fluxmin,*co_STY_p_f_L);
    if(Lleft<0) Lleft=0;
    Lright=Lum(z_STY_p_f_L,fsel_STY_p_f_L.fluxcut+5*fsel_STY_p_f_L.deltaflux,*co_STY_p_f_L);
    xright=Lright/lfamo.Lstar;
    xleft=Lleft/lfamo.Lstar;
/*     Llow=Lum(y[i],fluxl_STY_p_f_L,*co_STY_p_f_L); */
/*     xmin=Llow/lfamo.Lstar; */
    if(xleft>=150) {
      logL+=GSL_LOG_DBL_MAX; 
    }
    else {
      intabajo=0;
      /* Primero la integral a pelo de Lleft a Lright */
      for(j=0;j<nL_fs;j++) {
        L=Lleft+j*(Lright-Lleft)/(nL_fs-1.);
        fluxtmp=Flux(z_STY_p_f_L,L,*co_STY_p_f_L);
        intabajo=intabajo+exp(Schechter_L(L,lfamo))*Fermi(fluxtmp,fsel_STY_p_f_L.fluxcut,-fsel_STY_p_f_L.deltaflux);
      }
      intabajo=intabajo/nL_fs*(Lright-Lleft);  
      /* Ahora la integral desde Lright a inf */
      intabajo+=lfamo.phistar*(intsup-incom(1+lfamo.alfa,xright));
      probabajo=intabajo;

      logprobarriba=Schechter_L(Lumi,lfamo) + log(Fermi(flux_STY_p_f_L,fsel_STY_p_f_L.fluxcut,-fsel_STY_p_f_L.deltaflux));


      if(Fermi(flux_STY_p_f_L,fsel_STY_p_f_L.fluxcut,-fsel_STY_p_f_L.deltaflux)==0 || probabajo==0) logL+=GSL_LOG_DBL_MAX;
      else logL-= logprobarriba  - log(probabajo);   /* Perfectamente testado */
      /*       if(DEBUG3) printf(" iobj %d logL %f loglold %f      sch %g    pa %g pb %g (%g)  xmin %g x1 %g x2 %g flux %g err %g\n",i,logL,logLold,Schechter_L(Lumi,lfamo),log(probarriba),log(probabajo),probabajo,xmin,xleft,xright,flux_STY_p_f_L,errflux_STY_p_f_L); */
      /*       printf(" iobj %d logL %f loglold %f      sch %g int %g (%g)   pa %g pb %g (%g)  xmin %g x1 %f x2 %f\n",i,logL,logLold,Schechter_L(Lumi,lfamo),log(intsch),intsch,log(probarriba),log(probabajo),probabajo,xmin,x1,x2); */
    }
  }
  /* Aqui viene la parte de la poissoniana de npob */
  Ntot=Int_sch_f_L(lfamo,zlow_STY_p_f_L,zup_STY_p_f_L,fsel_STY_p_f_L,*co_STY_p_f_L)*strrad_STY_p_f_L/4./M_PI; 
  logL-=    (ndata*log(Ntot) - Ntot - gammln((double)ndata+1.)); 
/*   logLold-= (ndata*log(Ntot) - Ntot - gammln((double)ndata)+1.);  */
  
  if(DEBUG2) printf(" NTOT %f ndata*log(Ntot) %f gamm %f\n",Ntot,ndata*log(Ntot),gammln((double)ndata)+1.);

  iter_m++;
  if(DEBUG) printf(" Iter %d  logL %f logLold %f par0 %g par1 %g par2 %g\n",iter_m,logL,logLold,p[0],p[1],p[2]);
  return(logL); 
}


double Amoe_Funk_STY_p_f_L_conf(int n, double *x, double *y, double *p) {

/*   printf(" MLM %g conf %f   %f %f\n",MLmax,conflim,fabs(Amoe_Funk_STY_p_f_L_main(n,x,y,p)-(MLmax-log(conflim))),Amoe_Funk_STY_p_f_L_main(n,x,y,p)); */
   
  return(fabs(Amoe_Funk_STY_p_f_L_main(n,x,y,p)-(MLmax-log(conflim)))); 
}



void   EmpiricalCovars_STY_p_f_L(int n,double *flux,double *z,double *par, double *sigpar,struct fermifsel_L fsel, struct cosmo_param cosmo,struct Schlf_L *lf) {

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

  for(i=0;i<nconfl;i++) {
    printf("#");
    fflush(NULL);
    parconf[0]=par[0]+2*sigpar[0]*Gasdev();
    parconf[1]=par[1]+2*sigpar[1]*Gasdev();
    parconf[2]=par[2]+2*sigpar[2]*Gasdev();
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
    iter_c=Amoeba_d(n,flux,z,3,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_STY_p_f_L_conf);
    if(DEBUG) printf(" %d SOL iter %d par0 %f    par1 %f\n", i,  iter_c ,parconf[0],   parconf[1]);
    (parelip[0])[i]=parconf[0];
    (parelip[1])[i]=parconf[1];
    (parelip[2])[i]=parconf[2];
    if(iter_c==0 && Amoe_Funk_STY_p_f_L_conf(n,flux,z,parconf)>FTOL2 ) {
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


