#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_machine.h>
#include "alloc.h"
#include "cosmology.h"
#include "schechter.h"
#include "minmax.h"

/* #define FTOL 1e-14  */
#define FTOL 1e-17
#define MAXITER 500
#define DEBUG 0
#define DEBUG2 0
#define NBOOT 1000
#define NCONFL 5000

int ML_g_g_corr(int n,float *x,float *errx,float *mean,float *sigma,float *errmean,float *errsigma,float *covarmeansigma);
 
 
float Amoe_Funk_d1(int n, float *x, float *y, float *p);
float Amoe_Funk_d2(int n, float *x, float *y, float *p);
float Amoe_Funk_d3(int n, float *x, float *y, float *p);
float Amoe_Funk_conf(int n, float *x, float *y, float *p);
float Amoe_Funk(int n, float *x, float *y, float *p);
void  Covars_g_g_corr(int n,float *x,float *errx,float mean,float sigma,float covar[2][2]);
float *sig;
int ndata;
int iter;
int amoeba_use;
float ML[NBOOT];
float  m[NBOOT];
float  s[NBOOT];
float conflim;
float parconf0;
float parconf1;


int ML_g_g_corr(int n,float *x,float *errx,float *mean,float *sigma,float *errmean,float *errsigma,float *covarmeansigma) {

  float par[2];
  float sigpar[2];
  float parconf[2];
  
  float par0_conf_min;
  float par0_conf_max;
  float par1_conf_min;
  float par1_conf_max;

  float *y;
  int i;
  float covar[2][2];

/*   printf(" Entra aqui\n"); */


  if((y=malloc(n*sizeof(float)))==NULL) printf("I cannot dimension y of %d elements \n",n);



  sig=errx;
  ndata=n;
  
/*   printf(" Entra sin problemas \n"); */

/*   printf(" par0 %f par1 %f\n",par[0],par[1]); */

  par[0]=StWeightMedia(n,x,errx,par+1); 
/*   printf(" Aqui\n"); */
  sigpar[0]=par[1]/sqrt(n);
  sigpar[1]=par[1]/sqrt(n);

/*   printf(" INI par0 %f par1 %f\n",par[0],par[1]); */
/*   printf(" INI sigpar0 %f sigpar1 %f\n",sigpar[0],sigpar[1]); */

/*   printf(" Lallmo ameoab\n");  */
  amoeba_use=1;
/*   printf(" Auiqn\n"); */
  iter=0;
  Amoeba(n,x,y,2,par,sigpar,FTOL,MAXITER,Amoe_Funk);
  if(DEBUG) printf("Amoeba  Mean %14.10f Sigma %14.10f\n",par[0],par[1]); 
 

  Covars_g_g_corr(n,x,errx,par[0],par[1],covar);


  if(DEBUG) printf(" Covarianzas:     Var2(mean)= %f        Var(mean,sigma)=%f\n",covar[0][0],covar[0][1]);
  if(DEBUG) printf("              Var(mean,sigma)=%f            Var2(sigma)=%f\n",covar[1][0],covar[1][1]);


  if(DEBUG) printf(" alla\n");


  
/*   cpgpage(); */
/*   cpgsch(3); */
/*   cpgsci(1); */
/*   cpgswin(par[0]-1.5*sigpar[0],par[0]+1.5*sigpar[0],par[1]-1.5*sigpar[1],par[1]+1.5*sigpar[1]); */
/*   cpgbox("BCTNS",0,0,"BCTNS",0,0); */
/*   for(i=0;i<iter;i++) {  */
/*      printf(" tam %f\n",3.*exp((-ML[i]+ML[iter-1])*3.));  */
/*      cpgsch(3.*exp((-ML[i]+ML[iter-1])*3.));  */
/*      cpgpt1(m[i],s[i],2);  */
/*    }  */

/*   cpgpt1(m[iter-1],s[iter-1],2); */
  
/*   cpgsch(.5); */
  amoeba_use=2;
  conflim=exp(-.5);             //Los puntos corresponderan a 1 sigma de desviacion para dist. normal en par
  parconf[0]=par[0]-5.;
  parconf1=par[1];
  Amoeba(n,x,y,1,parconf,sigpar,FTOL,MAXITER,Amoe_Funk);
  if(DEBUG) printf("Amoeba  Mean-sigma(mean) %14.10f \n",parconf[0]);
/*   cpgpt1(parconf[0],parconf1,3);  */

  if(parconf[0]<par[0])  {
    par0_conf_min=parconf[0];
    parconf[0]=par[0]+(par[0]-par0_conf_min);
  }
  else {
    par0_conf_max=parconf[0];
    parconf[0]=par[0]+(par[0]-par0_conf_max);
  }
  Amoeba(n,x,y,1,parconf,sigpar,FTOL,MAXITER,Amoe_Funk);
  if(DEBUG) printf("Amoeba  Mean+sigma(mean) %14.10f \n",parconf[0]);
/*   cpgpt1(parconf[0],parconf1,3);  */
  if(parconf[0]<par[0])  {
    par0_conf_min=parconf[0];
  }
  else {
    par0_conf_max=parconf[0];
  }
  if(DEBUG) printf(" Sigma(mean))= %f   sigma(mean)_teor %14.10f\n",(par0_conf_max-par0_conf_min)/2.,sqrt(covar[0][0]));


  amoeba_use=3;
  parconf0=par[0];
  parconf[0]=0.9*par[1];
  Amoeba(n,x,y,2,parconf,sigpar,FTOL,MAXITER,Amoe_Funk);
  if(DEBUG) printf("Amoeba  Sigma+sigma(sigma) %14.10f\n",parconf[0]);
/*   cpgpt1(parconf0,parconf[0],3);  */
  if(parconf[0]<par[1])  {
    par1_conf_min=parconf[0];
    parconf[0]=par[1]+(par[1]-par1_conf_min);
  }
  else {
    par1_conf_max=parconf[0];
    parconf[0]=par[1]+(par[1]-par1_conf_max);
  }
  Amoeba(n,x,y,1,parconf,sigpar,FTOL,MAXITER,Amoe_Funk);
  if(parconf[0]<par[1])  {
    par1_conf_min=parconf[0];
  }
  else {
    par1_conf_max=parconf[0];
  }
  if(DEBUG) printf("Amoeba  Sigma+sigma(sigma) %14.10f \n",parconf[0]);
/*   cpgpt1(parconf0,parconf[0],3);  */
  if(DEBUG) printf(" Sigma(sigma))= %f sigma(sigma)_teor %14.10f\n",(par1_conf_max-par1_conf_min)/2.,sqrt(covar[1][1]));

  *mean=par[0];
  *sigma=par[1];


  return 0;

  for(i=0;i<NCONFL;i++) {
    parconf[0]=par[1];
    parconf0=par0_conf_min+(par0_conf_max-par0_conf_min)*i/(NCONFL-1.);
    Amoeba(n,x,y,2,parconf,sigpar,FTOL,MAXITER,Amoe_Funk);
/*     cpgpt1(parconf0,parconf[0],3);  */
 
  }





  amoeba_use=2;
  conflim=0.95;
  parconf[0]=par[0]-5.;
  parconf1=par[1];
  Amoeba(n,x,y,1,parconf,sigpar,FTOL,MAXITER,Amoe_Funk);
  if(DEBUG) printf("Amoeba  Mean %14.10f Sigma %14.10f\n",parconf[0],parconf1);
/*   cpgpt1(parconf[0],parconf1,3);  */

  if(parconf[0]<par[0])  {
    par0_conf_min=parconf[0];
    parconf[0]=par[0]+(par[0]-par0_conf_min);
  }
  else {
    par0_conf_max=parconf[0];
    parconf[0]=par[0]-(par[0]-par0_conf_min);
  }
  Amoeba(n,x,y,1,parconf,sigpar,FTOL,MAXITER,Amoe_Funk);
  if(DEBUG) printf("Amoeba  Mean %14.10f Sigma %14.10f\n",parconf[0],parconf1);
/*   cpgpt1(parconf[0],parconf1,3);  */
  if(parconf[0]<par[0])  {
    par0_conf_min=parconf[0];
  }
  else {
    par0_conf_max=parconf[0];
  }
  
  amoeba_use=3;
  for(i=0;i<NCONFL;i++) {
    parconf[0]=par[1];
    parconf0=par0_conf_min+(par0_conf_max-par0_conf_min)*i/(NCONFL-1.);
    Amoeba(n,x,y,2,parconf,sigpar,FTOL,MAXITER,Amoe_Funk);
/*     cpgpt1(parconf0,parconf[0],3);  */
 
  }

  amoeba_use=2;
  conflim=0.99;
  parconf[0]=par[0]-5.;
  parconf1=par[1];
  Amoeba(n,x,y,1,parconf,sigpar,FTOL,MAXITER,Amoe_Funk);
  if(DEBUG) printf("Amoeba  Mean %14.10f Sigma %14.10f\n",parconf[0],parconf1);
/*   cpgpt1(parconf[0],parconf1,3);  */

  if(parconf[0]<par[0])  {
    par0_conf_min=parconf[0];
    parconf[0]=par[0]+(par[0]-par0_conf_min);
  }
  else {
    par0_conf_max=parconf[0];
    parconf[0]=par[0]-(par[0]-par0_conf_min);
  }
  Amoeba(n,x,y,1,parconf,sigpar,FTOL,MAXITER,Amoe_Funk);
  if(DEBUG) printf("Amoeba  Mean %14.10f Sigma %14.10f\n",parconf[0],parconf1);
/*   cpgpt1(parconf[0],parconf1,3);  */
  if(parconf[0]<par[0])  {
    par0_conf_min=parconf[0];
  }
  else {
    par0_conf_max=parconf[0];
  }
  
  amoeba_use=3;
  for(i=0;i<NCONFL;i++) {
    parconf[0]=par[1];
    parconf0=par0_conf_min+(par0_conf_max-par0_conf_min)*i/(NCONFL-1.);
    Amoeba(n,x,y,2,parconf,sigpar,FTOL,MAXITER,Amoe_Funk);
/*     cpgpt1(parconf0,parconf[0],3);  */
 
  }
  

/*   cpgopen("/xserve"); */
/*   cpgswin(par[0]-sigpar[0],par[0]+sigpar[0],0.,1.); */
/*   cpgbox("BCTNS",0,0,"BCTNS",0,0); */
/*   for(i=0;i<iter;i++) { */
/*     cpgpt1(m[i],exp(-ML[i]+ML[iter-1]),2); */
/*   } */
/*   cpgopen("/xserve"); */
/*   cpgswin(par[1]-sigpar[1],par[1]+sigpar[1],0.,1.); */
/*   cpgbox("BCTNS",0,0,"BCTNS",0,0); */
/*   for(i=0;i<iter;i++) { */
/*     cpgpt1(s[i],exp(-ML[i]+ML[iter-1]),2); */
/*   } */


/*   iter=Mrq(x,y,sigy,nbin,par,ipar,3,covarpar,&chisq,funmaglim);    */
  



  printf("Amoeba  Mean %14.10f Sigma %14.10f\n",par[0],par[1]);
  *mean=par[0];
  *sigma=par[1];

  free(y);
  return 0;
}


float Amoe_Funk(int n, float *x, float *y, float *p) { 

  int i;
  double logL=0.;

  if(amoeba_use==1) {
  /* p[0] es h media */
  /* p[1] es sigma_h */
    
    if(p[1]<=0) p[1]=-p[1]; 
    logL=0.;
    for(i=0;i<ndata;i++) {
      logL-=-0.5*log(sig[i]*sig[i]+p[1]*p[1]);
      logL-=-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+p[1]*p[1]);
    }
    if(DEBUG2) printf(" logL %f p[0] %f p[1] %f\n",logL,p[0],p[1]);
/*     printf(" logL %f p[0] %f p[1] %f\n",logL,p[0],p[1]);  */
    ML[iter]=logL;
    m[iter]=p[0];
    s[iter]=p[1];
    iter++;
    return(logL);

  }

  if(amoeba_use==2) {
    /* p[0] es h media */
    logL=0.;
    for(i=0;i<ndata;i++) {
      logL-=-0.5*log(sig[i]*sig[i]+parconf1*parconf1);
      logL-=-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+parconf1*parconf1);
    }
    if(DEBUG2) printf(" logL %f logL-logLmax %f p[0] %f parconf1 %f\n",logL,fabs(logL-(ML[iter-1]-log(conflim))),p[0],parconf1);
    return(fabs(logL-(ML[iter-1]-log(conflim))));
  }

  if(amoeba_use==3) {
    /* p[0] es sigma media */
    if(p[0]<=0) p[0]=-p[0]; 
    logL=0.;
    for(i=0;i<ndata;i++) {
      logL-=-0.5*log(sig[i]*sig[i]+p[0]*p[0]);
      logL-=-(x[i]-parconf0)*(x[i]-parconf0)/2./(sig[i]*sig[i]+p[0]*p[0]);
    }
    if(DEBUG2) printf(" logL %f logL-logLmax %f p[0] %f parconf0 %f\n",logL,fabs(logL-(ML[iter-1]-log(conflim))),p[0],parconf0);
    return(fabs(logL-(ML[iter-1]-log(conflim))));
  }
  return 0;
}

float Amoe_Funk_d1(int n, float *x, float *y, float *p) { 
  /* p[0] es h media */
  /* p[1] es sigma_h */

  int i;
  double logL=0.;
  
  if(p[1]<=0) p[1]=-p[1]; 
  logL=0.;
  for(i=0;i<ndata;i++) {
    logL-=-0.5*log(sig[i]*sig[i]+p[1]*p[1]);
    logL-=-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+p[1]*p[1]);
  }
  ML[iter]=logL;
  m[iter]=p[0];
  s[iter]=p[1];
  iter++;
  return(logL);
}


float Amoe_Funk_d2(int n, float *x, float *y, float *p) {
  /* p[0] es h media */
  /* p[1] es sigma_h */

  float pi=4*atan(1.); 
  static int iter=0;
  int i,j;
  double logL=0.,logLb=0.;
  int nstep=1000;
  double conv,xint,xmin,xmax;
  (void)y;
  (void)n;

  if(p[1]<=0) p[1]=-p[1]; 
  logL=0.;
  logLb=0.;
  for(i=0;i<ndata;i++) {
    conv=0;
    xmin=minf(x[i]-4*sig[i],p[0]-4*p[1]);
    xmax=maxf(x[i]+4*sig[i],p[0]+4*p[1]);
    for(j=0;j<nstep;j++) {
      xint=xmin+(xmax-xmin)*j/(nstep-1);
      conv+=gaussian(x[i],xint,sig[i])*gaussian(xint,p[0],p[1]); 
      if(DEBUG2)      printf(" j %d  x %f e %f h %f  conv %g g1 %g  g3 %g\n",j,x[i],sig[i],xint,conv,gaussian(x[i],xint,sig[i]),gaussian(xint,p[0],p[1]));  
    }
    conv*=(xmax-xmin)/(nstep-1.);
    logLb-=-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+p[1]*p[1]); 
    logLb-=-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+p[1]*p[1]); 
    logL-=log(conv);
    if(DEBUG2) printf("i %4d sig %10f x[i] %10f p[0] %10f p[1] %10f anal %10f int %10f diff %10g x1 %f x2 %f\n",i,sig[i],x[i],p[0],p[1],-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+p[1]*p[1])-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+p[1]*p[1]),log(conv),-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+p[1]*p[1])-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+p[1]*p[1])-log(conv),xmin,xmax);   
  }
  if(DEBUG) printf(" iter %d logL1 %f logL2 %f p0 %f p1 %f\n",iter,logL,logLb,p[0],p[1]);
  ML[iter]=logL;
  m[iter]=p[0];
  s[iter]=p[1];
  iter++;
  return(logL);
}

float Amoe_Funk_d3(int n, float *x, float *y, float *p) {
  /* p[0] es h media */
  /* p[1] es sigma_h */
 
  float pi=4*atan(1.); 
  int i,j,k;
  double logL=0.,logLb=0.;
  int nstep=100;
  double conv,xint,xmin,xmax;
  double sigint,sigmin,sigmax;
  (void)y;/* To avoid warning */
  (void)n;/* To avoid warning */

  if(p[1]<=0) p[1]=-p[1]; 
  logL=0.;
  logLb=0.;
  for(i=0;i<ndata;i++) {
    conv=0;
    xmin=minf(x[i]-4*sig[i],p[0]-4*p[1]);
    xmax=maxf(x[i]+4*sig[i],p[0]+4*p[1]);
    sigmin=4+xint-3*fabs(xint);
    sigmax=4+xint+3*fabs(xint);
    sigmin=4-3*5.;
    sigmax=4+3*5.;
    for(j=0;j<nstep;j++) {
      for(k=0;k<nstep;k++) {
	xint=xmin+(xmax-xmin)*j/(nstep-1);
	sigint=sigmin+(sigmax-sigmin)*k/(nstep-1);
/* 	conv+=gaussian(x[i],xint,sigint)*gaussian(sigint,4+xint,fabs(xint))*gaussian(xint,p[0],p[1]);  */
	conv+=gaussian(x[i],xint,sigint)*gaussian(sigint,4.,10.)*gaussian(xint,p[0],p[1]); 
	
	/*       conv+=gaussian(x[i],xint,sigint)*gaussian(xint,p[0],p[1]);  */
/*  	printf(" j %d k %d x %f e %f h %f s %f conv %g g1 %g g2 %g g3 %g\n",j,k,x[i],sig[i],xint,sigint,conv,gaussian(x[i],xint,sigint),gaussian(sigint,4+xint,fabs(xint)),gaussian(xint,p[0],p[1]));  */
      }
      
    }
/*     printf(" conv %f\n",conv); */
    conv*=(xmax-xmin)/(nstep-1.)*(sigmax-sigmin)/(nstep-1.);
/*     printf(" conv... %f\n",conv); */
    logLb-=-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+p[1]*p[1]); 
    logLb-=-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+p[1]*p[1]); 
    logL-=log(conv);
/*     printf("i %4d sig %10f x[i] %10f p[0] %10f p[1] %10f anal %10f int %10f diff %10g x1 %f x2 %f\n",i,sig[i],x[i],p[0],p[1],-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+p[1]*p[1])-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+p[1]*p[1]),log(conv),-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+p[1]*p[1])-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+p[1]*p[1])-log(conv),xmin,xmax);   */
  }
  printf(" iter %d logL1 %f logL2 %f p0 %f p1 %f\n",iter,logL,logLb,p[0],p[1]);
  
  iter++;
  return(logL);
}


void  Covars_g_g_corr(int n,float *x,float *errx,float mean,float sigma,float covar[2][2]) {
  
  int i;
  double dum1,dum2,dum3;
  dum1=0;
  dum2=0;
  dum3=0;
  for(i=0;i<n;i++) {
    dum1+=1/(sigma*sigma+errx[i]+errx[i]);
/*     dum2+=2.*sigma*sigma/(sigma*sigma+errx[i]*errx[i])/(sigma*sigma+errx[i]*errx[i]); */
/*     dum2+=-4.*sigma*sigma*(x[i]-mean)*(x[i]-mean)/(sigma*sigma+errx[i]*errx[i])/(sigma*sigma+errx[i]*errx[i])/(sigma*sigma+errx[i]*errx[i]); */
    dum2+=2.*sigma*sigma*(-sigma*sigma+errx[i]*errx[i]+2*(x[i]-mean)*(x[i]-mean))/(sigma*sigma+errx[i]*errx[i])/(sigma*sigma+errx[i]*errx[i])/(sigma*sigma+errx[i]*errx[i]);
    dum3+=-1./(sigma*sigma+errx[i]*errx[i])+1.*(x[i]-mean)*(x[i]-mean)/(sigma*sigma+errx[i]*errx[i])/(sigma*sigma+errx[i]*errx[i]);
  }

/*   if(DEBUG) printf(" dum3 %g dum2  %g sigma %f mean %f\n",dum3,dum2,sigma,mean); */

  covar[0][0]=1./dum1;
  covar[1][1]=1./dum2;
  covar[0][1]=0.;
  covar[1][0]=0.;
}

