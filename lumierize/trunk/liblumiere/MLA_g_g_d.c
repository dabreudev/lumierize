#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include <gsl/gsl_machine.h>
#include "alloc.h"
#include "mlhist.h"
#include "stmedia.h"
#include "amoeba.h"
#include "random.h"
#include "cpgplot.h"
#include "elip.h"

#define FTOL  1e-12
#define FTOL2 5e-6
#define MAXITER  300
#define MAXITER2 40
#define DEBUG 1
#define DEBUG2 0
#define DEBUG3 0
#define DEBUGPLOT 0
#define NCONFL 20


double Amoe_Funk_g_g_d_main(int n, double *x, double *y, double *p);
double Amoe_Funk_g_g_d_confp0(int n, double *x, double *y, double *p);
double Amoe_Funk_g_g_d_confp1(int n, double *x, double *y, double *p);
void   DebugPlot(int n,double *x,double *errx); 


double *sig_g_g;
int ndata_g_g;
int iter_m_g_g;
int iter_c_g_g;
int amoeba_use_g_g;
double MLVal_g_g[10*MAXITER];
double  m_g_g[10*MAXITER];
double  s_g_g[10*MAXITER];
double conflim_g_g;
double parconf0_g_g;
double parconf1_g_g;
double meanerrx_g_g,stderrx_g_g;
double wmean_g_g,wsigma_g_g;


int MLA_g_g_d(int n,double *x,double *errx,double *mean,double *sigma,double *errmean,double *errsigma,double *covarmeansigma) {

  double par[2];
  double sigpar[2];
  double parconf[2];
  double sigparconf[2];
/*   double par[2]; */
/*   double sigpar[2]; */
/*   double parconf[2]; */

  double pargood[2];
  
  double par0_conf_min = 0;
  double par0_conf_max = 0;
  double par1_conf_min = 0;
  double par1_conf_max = 0;

  double *y;
  int i;
  double covar[2][2]; 


  double par0elip[NCONFL];
  double par1elip[NCONFL];

  double a,b,c,d,f,e;
  double xcelip,ycelip,elipa,elipb,elipt; 
  int ctrl1,ctrl2;

  double rho;

  int iter_amo;

/*   printf(" Entra aqui\n"); */

  sig_g_g=errx;
  ndata_g_g=n;



  if(DEBUGPLOT) DebugPlot(n,x,errx);

/*   return(); */

  if((y=malloc(n*sizeof(double)))==NULL) { printf("I cannot dimension y of %d elements \n",n);exit(1);}






  iter_m_g_g=0;

  meanerrx_g_g=StMedia_d(n,errx,&stderrx_g_g);

  if(DEBUG) printf(" meanerrx %g /500 %g\n",meanerrx_g_g,meanerrx_g_g/500.);

  if(DEBUG3) {
    for(i=0;i<n;i++) printf(" Entrada x %g errx %g\n",x[i],errx[i]);
  }


  iter_amo=MAXITER+1;
  printf(" iter_amo %d\n",iter_amo);


  while(iter_amo>= MAXITER) {
    wmean_g_g=StErrWeightMedia_d(n,x,errx,&wsigma_g_g);
    par[0]=wmean_g_g;
    par[1]=wsigma_g_g;
    sigpar[0]=par[1]/sqrt(n);
    sigpar[1]=par[1]/sqrt(n)/10.;
    if(DEBUG2) printf(" INI    par0 %f    par1 %f\n",   par[0],   par[1]); 
    if(DEBUG2) printf(" INI sigpar0 %f sigpar1 %f\n",sigpar[0],sigpar[1]); 
    iter_amo=Amoeba_d(n,x,y,2,par,sigpar,FTOL,MAXITER,Amoe_Funk_g_g_d_main);
  }
  pargood[0]=par[0];  pargood[1]=par[1];
/*   pargood[2]=par[2];  pargood[3]=par[3]; */
  if(DEBUG) printf(" logL Final: %20f\n",Amoe_Funk_g_g_d_main(n,x,y,pargood));
  if(DEBUG) printf("Amoeba  Mean %14.10f Sigma %14.10f\n",par[0],par[1]); 


  *mean=par[0];
  *sigma=par[1];
  pargood[0]=16.5;
  pargood[1]=3.;
  
  conflim_g_g=exp(-.5/100.);             //Los puntos corresponderan a 1 sigma entre 10 de desviacion para dist. normal en par

  parconf[0]=par[0]-5*par[0]/sqrt(n);
  sigparconf[0]=par[0]/sqrt(n);
  parconf1_g_g=par[1];
  iter_c_g_g=0;
  if(DEBUG) printf(" ml[iter_m-1] %.8g\n",MLVal_g_g[iter_m_g_g-1]);
  Amoeba_d(n,x,y,1,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_g_g_d_confp0);
  par0elip[0]=parconf[0];
  par1elip[0]=parconf1_g_g;
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
  if(DEBUG) printf(" par_conf_min %g  par_conf_max %g parconf[0] %g \n",par0_conf_min,par0_conf_max,parconf[0]);
  iter_c_g_g=0;
  Amoeba_d(n,x,y,1,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_g_g_d_confp0);
  par0elip[1]=parconf[0];
  par1elip[1]=parconf1_g_g;
  if(DEBUG) printf("Amoeba  Mean+sigma(mean) %14.10f \n",parconf[0]);
/*   cpgpt1(parconf[0],parconf1,3);  */
  if(parconf[0]<par[0])  {
    par0_conf_min=parconf[0];
  }
  else {
    par0_conf_max=parconf[0];
  }
  if(DEBUG) printf(" par_conf_min %g  par_conf_max %g parconf[0] %g \n",par0_conf_min,par0_conf_max,parconf[0]);

  if(DEBUG2) printf(" Sigma(mean))= %f   sigma(mean)_teor %14.10f\n",(par0_conf_max-par0_conf_min)/2.,sqrt(covar[0][0]));


  parconf[1]=par[1]-5*par[1]/sqrt(n);
  sigparconf[1]=par[1]/sqrt(n);
  parconf0_g_g=par[0];
  iter_c_g_g=0;
  if(DEBUG) printf(" ml[iter_m-1] %.8g\n",MLVal_g_g[iter_m_g_g-1]);
  Amoeba_d(n,x,y,1,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_g_g_d_confp1);
  par0elip[2]=parconf0_g_g;
  par1elip[2]=parconf[0];
  if(DEBUG) printf("Amoeba  Sgigma-sigma(sigma) %14.10f \n",parconf[0]);
/*   cpgpt1(parconf[0],parconf1,3);  */

  if(parconf[0]<par[1])  {
    par1_conf_min=parconf[0];
    parconf[0]=par[1]+(par[1]-par1_conf_min);
  }
  else {
    par1_conf_max=parconf[0];
    parconf[0]=par[1]+(par[1]-par1_conf_max);
  }
  if(DEBUG) printf(" par_conf_min %g  par_conf_max %g parconf[1] %g \n",par1_conf_min,par1_conf_max,parconf[0]);
  iter_c_g_g=0;
  Amoeba_d(n,x,y,1,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_g_g_d_confp1);
  par0elip[3]=parconf0_g_g;
  par1elip[3]=parconf[0];
  if(DEBUG) printf("Amoeba  Sigma+sigma(sigma) %14.10f \n",parconf[0]);
/*   cpgpt1(parconf[0],parconf1,3);  */
  if(parconf[0]<par[1])  {
    par1_conf_min=parconf[0];
  }
  else {
    par1_conf_max=parconf[0];
  }
  if(DEBUG) printf(" par_conf_min %g  par_conf_max %g parconf[1] %g \n",par1_conf_min,par1_conf_max,parconf[0]);

  if(DEBUG2) printf(" Sigma(sigma))= %f   sigma(sigma)_teor %14.10f\n",(par1_conf_max-par1_conf_min)/2.,sqrt(covar[1][1]));

  for(i=4;i<NCONFL;i++) {
    parconf[0]=par[0]+par[1]*Gasdev();
    parconf1_g_g=par1_conf_min+(par1_conf_max-par1_conf_min)*i/(NCONFL-1.);
    iter_c_g_g=0;
    iter_amo=Amoeba_d(n,x,y,1,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_g_g_d_confp0);
    par0elip[i]=parconf[0];
    par1elip[i]=parconf1_g_g;
    if(DEBUG) printf(" i %d par0 %g par1 %g\n",i,par0elip[i],par1elip[i]);
    if(iter_amo==MAXITER2) i--;
    /*     cpgpt1(parconf0,parconf[0],3);  */
  }

  for(i=0;i<NCONFL;i++) {
    par0elip[i]-=par[0];
    par1elip[i]-=par[1];
  }

  /* Para deshacer el asunto de que haya escogido un limite de confidencia diferente */
  for(i=0;i<NCONFL;i++) {
    par0elip[i]/=sqrt(-2*log(conflim_g_g));
    par1elip[i]/=sqrt(-2*log(conflim_g_g));
  }




  ctrl1=MCElip_d(NCONFL,par0elip,par1elip,&a,&b,&c,&d,&f,&e);
  if(ctrl1) {
    if(DEBUG) printf(" El ajuste tuvo exito\n");
    ctrl2=ElipPar_d(a,b,c,d,f,e,&xcelip,&ycelip,&elipa,&elipb,&elipt); 
    if(ctrl2) {if(DEBUG) printf(" Parece que ha funiona el otro ademas\n");}
    else      {if(DEBUG) printf(" El otro ha fallado.\n");}
  }
  else {
    if(DEBUG) printf(" No hubo exito\n");
  }
  


  if(DEBUG) printf(" a %g b %g c %g d %g e %g f %g\n",a,b,c,d,e,f);

  if(DEBUG) printf(" xcelip %g ycelip %g elipa %g elipb %g elipt %g\n",xcelip,ycelip,elipa,elipb,elipt); 

  if(DEBUG) printf("Amoeba  Mean %14.10f Sigma %14.10f\n",par[0],par[1]);
  *mean=par[0];
  *sigma=par[1];
  if(ctrl1) {

    rho=-f/(2.*sqrt(a*b));
    if(DEBUG) printf(" rho %f\n",rho);

    *errmean=1./sqrt(a*(1-rho*rho));
    *errsigma=1./sqrt(b*(1-rho*rho));
    *covarmeansigma=rho**errmean**errsigma;
  }
  else {
    *errmean=0.;
    *errsigma=0.;
    *covarmeansigma=0.;
  }

  if(DEBUG) {
    cpgopen("/xserve");
    cpgswin(par[0]-2**errmean,par[0]+2**errmean,par[1]-2**errsigma,par[1]+2**errsigma); 
    cpgbox("BCTNS",0,0,"BCTNS",0,0); 
    for(i=0;i<NCONFL;i++) {
      if(DEBUG2) printf(" Elip %g %g \n",par0elip[i],par1elip[i]);
      cpgpt1(par0elip[i]+par[0],par1elip[i]+par[1],2);
    }
    
    cpgmove(elipb*sin(-elipt/180.*M_PI)+par[0],elipa*cos(-elipt/180.*M_PI)+par[1]); 
    for(i=0;i<100;i++) {
      cpgdraw(elipa*sin(i/99.*2.*M_PI)*cos(elipt/180.*M_PI)-elipb*cos(i/99.*2.*M_PI)*sin(elipt/180.*M_PI)+par[0]+xcelip,+elipa*sin(i/99.*2.*M_PI)*sin(elipt/180.*M_PI)+elipb*cos(i/99.*2.*M_PI)*cos(elipt/180.*M_PI)+par[1]+ycelip); 
    }
    cpgclos();
  }
  
  free(y);

  if(iter_amo>=MAXITER) return(2);
  return(0);
  
  
}


double Amoe_Funk_g_g_d_main(int n, double *x, double *y, double *p) {
  /* p[0] es h media */
  /* p[1] es sigma_h */

  double pi=4*atan(1.); 
  int i;
  double logL=0.;
  (void)y;/* To avoid warning */
  (void)n;/* To avoid warning */


/*   if((p[1])<=meanerrx/500.) p[1]=wsigma;  */
  logL=0.;
  for(i=0;i<ndata_g_g;i++) {
    if(DEBUG3) printf(" a ame_main le llega x %g errx %g\n",x[i],sig_g_g[i]);
    logL-=-log(sqrt(2.*pi))-0.5*log(sig_g_g[i]*sig_g_g[i]+p[1]*p[1]);
    logL-=-(x[i]-p[0])*(x[i]-p[0])/2./(sig_g_g[i]*sig_g_g[i]+p[1]*p[1]);
  }
  if(DEBUG2) printf(" mlF iter %3d logL %15g  p0 %f p1 %f ndata %d\n",iter_m_g_g,logL,p[0],p[1],ndata_g_g);

  MLVal_g_g[iter_m_g_g]=logL;
  m_g_g[iter_m_g_g]=p[0];
  s_g_g[iter_m_g_g]=p[1];
  iter_m_g_g++;
  return(logL);
}



double Amoe_Funk_g_g_d_confp0(int n, double *x, double *y, double *p) {
  /* p[0] es h media */
  /* parconf1 es sigma_h */

  double pi=4*atan(1.); 
  int i;
  double logL=0.;
  (void)y;/* To avoid warning */
  (void)n;/* To avoid warning */

  
  logL=0.;
  for(i=0;i<ndata_g_g;i++) {
    logL-=-log(sqrt(2.*pi))-0.5*log(sig_g_g[i]*sig_g_g[i]+parconf1_g_g*parconf1_g_g);
    logL-=-(x[i]-p[0])*(x[i]-p[0])/2./(sig_g_g[i]*sig_g_g[i]+parconf1_g_g*parconf1_g_g);
  }
  if(DEBUG2) printf(" iter %d logL %g log(conlim_g_g) %g p0 %f p1 %f\n",iter_c_g_g,fabs(logL-(MLVal_g_g[iter_m_g_g-1]-log(conflim_g_g))),log(conflim_g_g),p[0],parconf1_g_g);
  if(DEBUG2 && log(conflim_g_g)==0 ) printf(" mlF2 iter %3d logL %15g  p0 %f p1 %f\n",iter_c_g_g,logL,p[0],parconf1_g_g);
  iter_c_g_g++;
  return(fabs(logL-(MLVal_g_g[iter_m_g_g-1]-log(conflim_g_g))));
/*   return(fabs(logL-log(conflim))); */
}

double Amoe_Funk_g_g_d_confp1(int n, double *x, double *y, double *p) {
  /* parconf0 es h media */
  /* p[0] es sigma_h */

  double pi=4*atan(1.); 
  int i;
  double logL=0.;
  (void)y;/* To avoid warning */
  (void)n;/* To avoid warning */


  if(p[0]<=0) p[0]=-p[0]; 
  logL=0.;
  for(i=0;i<ndata_g_g;i++) {
    logL-=-log(sqrt(2.*pi))-0.5*log(sig_g_g[i]*sig_g_g[i]+p[0]*p[0]);
    logL-=-(x[i]-parconf0_g_g)*(x[i]-parconf0_g_g)/2./(sig_g_g[i]*sig_g_g[i]+p[0]*p[0]);
  }
  if(DEBUG2) printf(" iter %d logL %g p0 %f p1 %f\n",iter_c_g_g,fabs(logL-(MLVal_g_g[iter_m_g_g-1]-log(conflim_g_g))),parconf0_g_g,p[0]);
  if(DEBUG2) printf(" mlF3 iter %3d logL %15g  p0 %f p1 %f\n",iter_c_g_g,logL,p[0],p[1]);
  iter_c_g_g++;
  return(fabs(logL-(MLVal_g_g[iter_m_g_g-1]-log(conflim_g_g))));
/*   return(fabs(logL-log(conflim))); */
  
}


void DebugPlot(int n,double *x,double *errx) {

  int i,j;
  

  int ni=250;
  int nj=250;
  
  float minm=10.0,maxm=24.;
  float mins=-4,maxs=6.;

  float m,s;

  double ppar[2];
  double *y;

  static int nf=0;

  char fname[100]="surplot.dat";
  
  FILE *fp;
  (void)errx; //Avoid warning
  if((y=malloc(n*sizeof(double)))==NULL) { printf("I cannot dimension y of %d elements \n",n);exit(1);}

  sprintf(fname,"surplot.dat%d",nf);

  fp=fopen(fname,"w");

  for(i=0;i<ni;i++) {
    for(j=0;j<nj;j++) {
      m=minm+(maxm-minm)*i/(ni-1.);
      s=mins+(maxs-mins)*j/(nj-1.);
      ppar[0]=m;
      ppar[1]=s;
      fprintf(fp," %.20g %g %g\n",Amoe_Funk_g_g_d_main(n,x,y,ppar),ppar[0],ppar[1]);
      iter_m_g_g=0;
    }
  }
  
  fclose(fp);
  free(y);

  nf++;

}
