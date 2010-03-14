#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_machine.h>
#include "alloc.h"
#include "cosmology.h"
#include "schechter.h"
#include "minmax.h"
#include "mlhist.h"

#define FTOL  1e-12
#define FTOL2 5e-6
#define MAXITER  300
#define MAXITER2 60
#define DEBUG 0
#define DEBUG2 0
#define DEBUG3 0
#define DEBUGPLOT 0
#define NCONFL 20


double Amoe_Funk_g_n_f_d_main(int n, double *x, double *y, double *p);
double Amoe_Funk_g_n_f_d_conf(int n, double *x, double *y, double *p);
double Funk_norm_g_n_f_d(double x);
void   EmpiricalCovars_g_n_f_d(int n,double *x,double *par,double *sigpar,double xfermi, double Tfermi, double *mean,double *sigma,double *errmean,double *errsigma,double *covarmeansigma);



int ndata;
int iter_m;
int iter_c;
double MLmax;
double conflim;

double xf;
double Tf;

double p0,p1;

int MLA_g_n_f_d(int n,double *x,double xfermi, double Tfermi, double *mean,double *sigma,double *errmean,double *errsigma,double *covarmeansigma) {

  double par[2];
  double sigpar[2];
/*   double parconf[2]; */
/*   double sigparconf[2]; */

  double *y;

  double wmean,wsigma;


  int iter_amo;

/*   printf(" Entra aqui\n"); */

  ndata=n;
  xf=xfermi;
  Tf=Tfermi;



/*   return(); */

  if((y=malloc(n*sizeof(double)))==NULL) { printf("I cannot dimension y of %d elements \n",n);exit(1);}

  iter_m=0;

  iter_amo=MAXITER+1;
  if(DEBUG)   printf(" iter_amo %d\n",iter_amo);
  
 

  iter_amo=0;
  while(iter_amo==0) {
    wmean=StMedia_d(n,x,&wsigma); 
    par[0]=wmean;
    par[1]=wsigma;
    sigpar[0]=par[1]/sqrt(n);
    if(maxf(fabs(xf-Tf-par[0]),fabs(xf+Tf-par[0]))<5) sigpar[0]=maxf(fabs(xf-Tf-par[0]),fabs(xf+Tf-par[0]))/3.;
    sigpar[1]=sqrt(par[1]);
    if(DEBUG2) printf(" INI    par0 %f    par1 %f\n",   par[0],   par[1]); 
    if(DEBUG2) printf(" INI sigpar0 %f sigpar1 %f\n",sigpar[0],sigpar[1]); 
    iter_amo=Amoeba_d(n,x,y,2,par,sigpar,FTOL,MAXITER,Amoe_Funk_g_n_f_d_main);
  }
  MLmax=Amoe_Funk_g_n_f_d_main(n,x,y,par);
  if(DEBUG) printf(" logL Final: %20f iter %d\n",MLmax,iter_amo);
  if(DEBUG) printf("Amoeba  Mean %14.10f Sigma %14.10f\n",par[0],par[1]); 


  *mean=par[0];
  *sigma=par[1];

  conflim=exp(-.5/10.);
  EmpiricalCovars_g_n_f_d( n, x, par,sigpar, xfermi,Tfermi, mean,sigma,errmean,errsigma,covarmeansigma); 



  free(y);

  if(iter_amo>=MAXITER) return(2);
  return(0);
  
  
}


double Amoe_Funk_g_n_f_d_main(int n, double *x, double *y, double *p) {
  /* p[0] es h media */
  /* p[1] es sigma_h */

  int i,j;
  double logL=0.,logLb;
  double conv=0;
  double norm;
  int nstepnorm=1000;
  double xminnorm,xmaxnorm;
  double xintnorm;  
  double scale,offset;
  double gauss;
  
  int nul=0;
  int nul2=0;
  (void)y;/* To avoid warning */
  (void)n;/* To avoid warning */


  if((p[1])<0) p[1]=-p[1];  

  
  if(p[1]/Tf>3) {
    xminnorm=p[0]-4*p[1];
    xmaxnorm=p[0]+4*p[1];
    nstepnorm=maxi(10*p[1]/Tf,nstepnorm);
    if(nstepnorm>10000) nstepnorm=10000;
  }
  else if(Tf/p[1]>3) {
    xminnorm=p[0]-4*p[1];
    xmaxnorm=p[0]+4*p[1];
  }
  else {
    xminnorm=minf(p[0]-4*p[1],xf-2*Tf);
    xmaxnorm=maxf(p[0]+4*p[1],xf+2*Tf);
    nstepnorm=maxi(100*p[1]/Tf,nstepnorm);
  }
  norm=0;
  for(j=0;j<nstepnorm;j++) {
    xintnorm=xminnorm+(xmaxnorm-xminnorm)*j/(nstepnorm-1);
    norm+=gaussian(xintnorm,p[0],p[1])*Fermi(xintnorm,xf,Tf);
  }
  norm*=(xmaxnorm-xminnorm)/(nstepnorm-1.);

  p0=p[0];
  p1=p[1];
  scale=p[1];
  offset=p[0];
  gauss=gaussinther_d(Funk_norm_g_n_f_d,offset,scale,30);
  /* Cuando demuestre que gaussinther va mejor, la uso */

  

/*   if((p[1])<=meanerrx/500.) p[1]=wsigma;  */
  logL=0.;
  logLb=0.;
  for(i=0;i<ndata;i++) {
    if(1./Fermi(x[i],xf,Tf)>1e5) {
      conv+=0;
      nul++;
      if(DEBUG) printf(" caso 1 \n");
/*       exit(1);  */
    }
    else         {
      conv=gaussian(x[i],p[0],p[1])*Fermi(x[i],xf,Tf)/norm;
      if(DEBUG3) printf(" obj %d conv %f f %f norm %f cc %f xi %f\n",i,conv,(0.5*(erf((xf-p[0])/sqrt(2)/p[1])+1)),norm,pow(conv,1./Fermi(x[i],xf,Tf)),x[i]);
      /*       conv=pow(conv,1./Fermi(x[i],xf,Tf));  */
    }
    if(conv==0) {
      logL-=32;
      nul2++;
      if(DEBUG) printf(" caso 2 \n");
    }
    else            logL-=log(conv);
  }
  if(DEBUG)  printf(" %d logL %f loglb %f m %f s %f conv %f nul %d nul2 %d norm %f\n",iter_m,logL,logLb,p[0],p[1],conv,nul,nul2,norm);

  iter_m++;
  return(logL);
}

double Amoe_Funk_g_n_f_d_conf(int n, double *x, double *y, double *p) {
  /* p[0] es h media */
  /* p[1] es sigma_h */

  int i,j;
  double logL=0.,logLb;
  double conv;
  double norm;
  int nstepnorm=1000;
  double xminnorm,xmaxnorm;
  double xintnorm;  
  
  int nul=0;
  int nul2=0;
  (void)y;/* To avoid warning */
  (void)n;/* To avoid warning */

  if((p[1])<0) p[1]=-p[1];  

  if(p[1]/Tf>3) {
    xminnorm=p[0]-4*p[1];
    xmaxnorm=p[0]+4*p[1];
    nstepnorm=maxi(10*p[1]/Tf,nstepnorm);
    if(nstepnorm>10000) nstepnorm=10000;
  }
  else if(Tf/p[1]>3) {
    xminnorm=p[0]-4*p[1];
    xmaxnorm=p[0]+4*p[1];
  }
  else {
    xminnorm=minf(p[0]-4*p[1],xf-2*Tf);
    xmaxnorm=maxf(p[0]+4*p[1],xf+2*Tf);
    nstepnorm=maxi(100*p[1]/Tf,nstepnorm);
  }
  norm=0;
  for(j=0;j<nstepnorm;j++) {
    xintnorm=xminnorm+(xmaxnorm-xminnorm)*j/(nstepnorm-1);
    norm+=gaussian(xintnorm,p[0],p[1])*Fermi(xintnorm,xf,Tf);
  }
  norm*=(xmaxnorm-xminnorm)/(nstepnorm-1.);

  
  

/*   if((p[1])<=meanerrx/500.) p[1]=wsigma;  */
  logL=0.;
  logLb=0.;
  for(i=0;i<ndata;i++) {
    if(1./Fermi(x[i],xf,Tf)>1e5) {
      conv+=0;
      nul++;
    }
    else         {
      conv=gaussian(x[i],p[0],p[1])*Fermi(x[i],xf,Tf)/norm;
      if(DEBUG3) printf(" obj %d conv %f f %f norm %f cc %f xi %f\n",i,conv,(0.5*(erf((xf-p[0])/sqrt(2)/p[1])+1)),norm,pow(conv,1./Fermi(x[i],xf,Tf)),x[i]);
    }
    if(conv==0) {
      logL-=32;
      nul2++;
    }
    else            logL-=log(conv);
  }

  return(fabs(logL-(MLmax-log(conflim))));
}




void   EmpiricalCovars_g_n_f_d(int n,double *x,double *par,double *sigpar,double xfermi, double Tfermi, double *mean,double *sigma,double *errmean,double *errsigma,double *covarmeansigma) {

  double *parconf;
  double *sigparconf;
  double **invcovar;
  double **covar;
  double **parelip;
  double *y;

  double **bb;
  int nconfl,nconflini;
  int i,j;
  double first, median, third, *distmax;

  nconfl=NCONFL;
  nconflini=nconfl;

  /* Dimensiono las matrices */

  bb=matrix_d(2,1);
  y=vector_d(n);
  parconf=vector_d(2);
  sigparconf=vector_d(2);
  parelip=matrix_d(2,nconflini);
  invcovar=matrix_d(2 ,2 );
  covar=matrix_d(2 ,2 );
  distmax=vector_d(2);


  for(i=0;i<nconfl;i++) {
    parconf[0]=par[0]+3*sigpar[0]*Gasdev();
    parconf[1]=par[1]+3*sigpar[1]*Gasdev();
    sigparconf[0]=sigpar[0]; 
    sigparconf[1]=sigpar[1]; 
    if(i>(int)(nconfl/2.)) {
      parconf[0]=par[0]-((parelip[0])[(int)(i-nconfl/2.)+1]-par[0]);
      parconf[1]=par[1]-((parelip[1])[(int)(i-nconfl/2.)+1]-par[1]);
    }
    iter_c=0;
    iter_c=Amoeba_d(n,x,y,2,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_g_n_f_d_conf);
    if(DEBUG) printf(" %d SOL par0 %f    par1 %f\n", i,  parconf[0],   parconf[1]); 
    (parelip[0])[i]=parconf[0];
    (parelip[1])[i]=parconf[1];
    if(iter_c==0 && Amoe_Funk_g_n_f_d_conf(n,x,y,parconf)>FTOL2 ) i--;
  }
  
  /* Supongo que el centro de la elipse es el valor que maximiza ML */
  for(i=0;i<nconfl;i++) {
    (parelip[0])[i]-=par[0];
    (parelip[1])[i]-=par[1];
  }

  /* Detecto puntos que esten muy alejados de la supuesta elipse */
  for(j=0;j<2;j++) {
    Quartil_d(nconfl,parelip[j],&first,&median,&third);
    distmax[j]=maxf(fabs(first),fabs(third));
  }  
  for(i=0;i<nconfl;i++) {
    for(j=0;j<2;j++) {
      if(fabs((parelip[j])[i])>2*2*distmax[j]/1.35) {
	for(j=0;j<2;j++) {
	  memcpy(parelip[j]+i,parelip[j]+i+1,(nconfl-i-1)*sizeof(double));
	}
	i--;
	nconfl--;
	break;
      }
    }
  }

  MCElipN_d(nconfl,2,parelip,invcovar);
  for(i=0;i<2;i++) {
    for(j=0;j<2;j++) {
      covar[i][j]=invcovar[i][j];
    }
  } 
  gaussj_d(covar,2,bb,1);
  /* Deshago el cambio para el limite de confidencia */
  for(i=0;i<2;i++) {
    for(j=0;j<2;j++) {
      covar[i][j]/=(-2*log(conflim));
    }
  } 

  *errmean=sqrt(covar[0][0]);
  *errsigma=sqrt(covar[1][1]);
  *covarmeansigma=covar[0][1];
  
  free_matrix_d(bb,2,1);
  free(y);
  free(parconf);
  free(sigparconf);

  free_matrix_d(parelip,2,nconflini);
  free_matrix_d(invcovar,2  ,2);
  free_matrix_d(   covar,2  ,2);
 
}


double Funk_norm_g_n_f_d(double x) {
  return(gaussian(x,p0,p1)*Fermi(x,xf,Tf));
}



