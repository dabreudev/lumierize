#include "modulos.h"
#define FTOL  1e-12
#define FTOL2 5e-6
#define MAXITER  300
#define MAXITER2 40
#define DEBUG 0
#define DEBUG2 0
#define DEBUG3 0
#define DEBUGPLOT 0
#define NCONFL 20


float Amoe_Funk_main(int n, float *x, float *y, float *p);
float Amoe_Funk_confp0(int n, float *x, float *y, float *p);
float Amoe_Funk_confp1(int n, float *x, float *y, float *p);
void   DebugPlot(int n,float *x,float *errx); 


float *sig;
int ndata;
int iter_m;
int iter_c;
int amoeba_use; 
float ML[10*MAXITER];
float  m[10*MAXITER];
float  s[10*MAXITER];
float conflim;
float parconf0;
float parconf1;
float meanerrx,stderrx;
float wmean,wsigma;


int MLA_g_g(int n,float *x,float *errx,float *mean,float *sigma,float *errmean,float *errsigma,float *covarmeansigma) {

  float par[2];
  float sigpar[2];
  float parconf[2];
  float sigparconf[2];
/*   float par[2]; */
/*   float sigpar[2]; */
/*   float parconf[2]; */

  float pargood[2];
  
  float par0_conf_min;
  float par0_conf_max;
  float par1_conf_min;
  float par1_conf_max;

  float *y;
  int i;
  float covar[2][2]; 


  float par0elip[NCONFL];
  float par1elip[NCONFL];

  float a,b,c,d,f,e;
  float xcelip,ycelip,elipa,elipb,elipt; 
  int ctrl1,ctrl2;

  float rho;

  int iter_amo;

/*   printf(" Entra aqui\n"); */

  sig=errx;
  ndata=n;



  if(DEBUGPLOT) DebugPlot(n,x,errx);

/*   return(); */

  if((y=malloc(n*sizeof(float)))==NULL) { printf("I cannot dimension y of %d elements \n",n);exit(1);}






  iter_m=0;

  meanerrx=StMedia(n,errx,&stderrx);

  if(DEBUG) printf(" meanerrx %g /500 %g\n",meanerrx,meanerrx/500.);

  if(DEBUG3) {
    for(i=0;i<n;i++) printf(" Entrada x %g errx %g\n",x[i],errx[i]);
  }


  iter_amo=MAXITER+1;
  if(DEBUG) printf(" iter_amo %d\n",iter_amo);


  while(iter_amo>= MAXITER) {
    wmean=StErrWeightMedia(n,x,errx,&wsigma); 
    par[0]=wmean;
    par[1]=wsigma;
    sigpar[0]=par[1]/sqrt(n);
    sigpar[1]=par[1]/sqrt(n)/10.;
    if(DEBUG2) printf(" INI    par0 %f    par1 %f\n",   par[0],   par[1]); 
    if(DEBUG2) printf(" INI sigpar0 %f sigpar1 %f\n",sigpar[0],sigpar[1]); 
    iter_amo=Amoeba(n,x,y,2,par,sigpar,FTOL,MAXITER,Amoe_Funk_main);
  }
  pargood[0]=par[0];  pargood[1]=par[1];
/*   pargood[2]=par[2];  pargood[3]=par[3]; */
  if(DEBUG) printf(" logL Final: %20f\n",Amoe_Funk_main(n,x,y,pargood));
  if(DEBUG) printf("Amoeba  Mean %14.10f Sigma %14.10f\n",par[0],par[1]); 


  *mean=par[0];
  *sigma=par[1];
  pargood[0]=16.5;
  pargood[1]=3.;
  
  conflim=exp(-.5/100.);             //Los puntos corresponderan a 1 sigma entre 10 de desviacion para dist. normal en par

  parconf[0]=par[0]-5*par[0]/sqrt(n);
  sigparconf[0]=par[0]/sqrt(n);
  parconf1=par[1];
  iter_c=0;
  if(DEBUG) printf(" ml[iter_m-1] %.8g\n",ML[iter_m-1]);
  Amoeba(n,x,y,1,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_confp0);
  par0elip[0]=parconf[0];
  par1elip[0]=parconf1;
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
  iter_c=0;
  Amoeba(n,x,y,1,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_confp0);
  par0elip[1]=parconf[0];
  par1elip[1]=parconf1;
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
  parconf0=par[0];
  iter_c=0;
  if(DEBUG) printf(" ml[iter_m-1] %.8g\n",ML[iter_m-1]);
  Amoeba(n,x,y,1,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_confp1);
  par0elip[2]=parconf0;
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
  iter_c=0;
  Amoeba(n,x,y,1,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_confp1);
  par0elip[3]=parconf0;
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
    parconf1=par1_conf_min+(par1_conf_max-par1_conf_min)*i/(NCONFL-1.);
    iter_c=0;
    iter_amo=Amoeba(n,x,y,1,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_confp0);
    par0elip[i]=parconf[0];
    par1elip[i]=parconf1;
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
    par0elip[i]/=sqrt(-2*log(conflim));
    par1elip[i]/=sqrt(-2*log(conflim));
  }




  ctrl1=MCElip(NCONFL,par0elip,par1elip,&a,&b,&c,&d,&f,&e);
  if(ctrl1) {
    if(DEBUG) printf(" El ajuste tuvo exito\n");
    ctrl2=ElipPar(a,b,c,d,f,e,&xcelip,&ycelip,&elipa,&elipb,&elipt); 
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


float Amoe_Funk_main(int n, float *x, float *y, float *p) {
  /* p[0] es h media */
  /* p[1] es sigma_h */

  float pi=4*atan(1.); 
  int i;
  double logL=0.;

/*   if((p[1])<=meanerrx/500.) p[1]=wsigma;  */
  logL=0.;
  for(i=0;i<ndata;i++) {
    if(DEBUG3) printf(" a ame_main le llega x %g errx %g\n",x[i],sig[i]);
    logL-=-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+p[1]*p[1]);
    logL-=-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+p[1]*p[1]);
  }
  if(DEBUG2) printf(" mlF iter %3d logL %15g  p0 %f p1 %f ndata %d\n",iter_m,logL,p[0],p[1],ndata);

  ML[iter_m]=logL; 
  m[iter_m]=p[0]; 
  s[iter_m]=p[1]; 
  iter_m++;
  return((float)logL);
}



float Amoe_Funk_confp0(int n, float *x, float *y, float *p) {
  /* p[0] es h media */
  /* parconf1 es sigma_h */

  float pi=4*atan(1.); 
  int i;
  double logL=0.;
  
  logL=0.;
  for(i=0;i<ndata;i++) {
    logL-=-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+parconf1*parconf1);
    logL-=-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+parconf1*parconf1);
  }
  if(DEBUG2) printf(" iter %d logL %g log(conlim) %g p0 %f p1 %f\n",iter_c,fabs(logL-(ML[iter_m-1]-log(conflim))),log(conflim),p[0],parconf1);
  if(DEBUG2 && log(conflim)==0 ) printf(" mlF2 iter %3d logL %15g  p0 %f p1 %f\n",iter_c,logL,p[0],parconf1);
  iter_c++;
  return(fabs(logL-(ML[iter_m-1]-log(conflim)))); 
/*   return(fabs(logL-log(conflim))); */
}

float Amoe_Funk_confp1(int n, float *x, float *y, float *p) {
  /* parconf0 es h media */
  /* p[0] es sigma_h */

  float pi=4*atan(1.); 
  int i;
  double logL=0.;

  if(p[0]<=0) p[0]=-p[0]; 
  logL=0.;
  for(i=0;i<ndata;i++) {
    logL-=-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+p[0]*p[0]); 
    logL-=-(x[i]-parconf0)*(x[i]-parconf0)/2./(sig[i]*sig[i]+p[0]*p[0]); 
  }
  if(DEBUG2) printf(" iter %d logL %g p0 %f p1 %f\n",iter_c,fabs(logL-(ML[iter_m-1]-log(conflim))),parconf0,p[0]);
  if(DEBUG2) printf(" mlF3 iter %3d logL %15g  p0 %f p1 %f\n",iter_c,logL,p[0],p[1]);
  iter_c++;
  return(fabs(logL-(ML[iter_m-1]-log(conflim)))); 
/*   return(fabs(logL-log(conflim))); */
  
}


void DebugPlot(int n,float *x,float *errx) {

  int i,j;
  

  int ni=250;
  int nj=250;
  
  float minm=10.0,maxm=24.;
  float mins=-4,maxs=6.;

  float m,s;

  float ppar[2];
  float *y;

  static int nf=0;

  char fname[100]="surplot.dat";
  
  FILE *fp;

  if((y=malloc(n*sizeof(float)))==NULL) { printf("I cannot dimension y of %d elements \n",n);exit(1);}

  sprintf(fname,"surplot.dat%d",nf);

  fp=fopen(fname,"w");

  for(i=0;i<ni;i++) {
    for(j=0;j<nj;j++) {
      m=minm+(maxm-minm)*i/(ni-1.);
      s=mins+(maxs-mins)*j/(nj-1.);
      ppar[0]=m;
      ppar[1]=s;
      fprintf(fp," %.20g %g %g\n",Amoe_Funk_main(n,x,y,ppar),ppar[0],ppar[1]);
      iter_m=0;
    }
  }
  
  fclose(fp);
  free(y);

  nf++;

}
