#include "modulos.h"
#define FTOL 1e-12
#define MAXITER 300
#define DEBUG 1
#define DEBUG2 0
#define NBOOT 1000
#define NCONFL 5000
 
int ML_g_g_corr_g(int n,float *x,float *errx,float *mean,float *sigma,float *errmean,float *errsigma,float *covarmeansigma);
 
 
float Amoe_Funk_dg1(int n, float *x, float *y, float *p);
float Amoe_Funk_dg2(int n, float *x, float *y, float *p);
float Amoe_Funk_dg3(int n, float *x, float *y, float *p);
float Amoe_Funk_gconf(int n, float *x, float *y, float *p);
float Amoe_Funk_main(int n, float *x, float *y, float *p);
float Amoe_Funk_confp0(int n, float *x, float *y, float *p);
float Amoe_Funk_confp1(int n, float *x, float *y, float *p);
void  Covars_g_g_corr_g(int n,float *x,float *errx,float mean,float sigma,float covar[2][2]);
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


int ML_g_g_corr_g(int n,float *x,float *errx,float *mean,float *sigma,float *errmean,float *errsigma,float *covarmeansigma) {

  float par[4];
  float sigpar[4];
  float parconf[4];
/*   float par[2]; */
/*   float sigpar[2]; */
/*   float parconf[2]; */

  float pargood[4];
  
  float par0_conf_min;
  float par0_conf_max;
  float par1_conf_min;
  float par1_conf_max;

  float *y;
  int i;
  float covar[4][4]; 
/*   float covar[2][2]; */

/*   float fdum; */

/*   printf(" Entra aqui\n"); */


  if((y=malloc(n*sizeof(float)))==NULL) printf("I cannot dimension y of %d elements \n",n);



  sig=errx;
  ndata=n;
  
/*   printf(" Entra sin problemas \n"); */

/*   printf(" par0 %f par1 %f\n",par[0],par[1]); */

  par[0]=StWeightMedia(n,x,errx,par+1); 

  par[2]=StMedia(n,errx,sigpar+2); 
  par[3]=par[1]/sigpar[2];

/*   covar=StCovar(n,x,errx,&sigmah,&sigmasigma);  */
      
/*   printf(" Aqui\n"); */
  sigpar[0]=par[1]/sqrt(n);
  sigpar[1]=par[1]/sqrt(n);
  sigpar[2]=sigpar[2]/sqrt(n);
  sigpar[3]=par[3]/2.;

  /* De prueba */

  par[0]=25.;par[1]=5;
  sigpar[0]=0.1;sigpar[1]=0.1;


  if(DEBUG) printf(" INI par0 %f par1 %f par2 %f par3 %f\n",par[0],par[1],par[2],par[3]); 
  if(DEBUG) printf(" INI sigpar0 %f sigpar1 %f sigpar2 %f sigpar3 %f\n",sigpar[0],sigpar[1],sigpar[2],sigpar[3]); 

/*   printf(" Lallmo ameoab\n");  */
  amoeba_use=1;
/*   printf(" Auiqn\n"); */
  iter=0;
/*   par[2]=4.; */
/*   par[3]=15.*15.; */

  pargood[0]=15.;
  pargood[1]=5.;
  pargood[2]=4.;
  pargood[3]=25.*25.;
  
  printf(" logL Buena: %20f\n",Amoe_Funk_main(n,x,y,pargood));




  Amoeba(n,x,y,2,par,sigpar,FTOL,MAXITER,Amoe_Funk_main);
  pargood[0]=par[0];  pargood[1]=par[1];  pargood[2]=par[2];  pargood[3]=par[3];
  printf(" logL Final: %20f\n",Amoe_Funk_main(n,x,y,pargood));
  if(DEBUG) printf("Amoeba  Mean %14.10f Sigma %14.10f\n",par[0],par[1]); 


  *mean=par[0];
  *sigma=par[1];
  pargood[0]=15.;
  pargood[1]=5.;
  pargood[2]=4.;
  pargood[3]=15.*15.;
  
  printf(" logL Buena: 20%f\n",Amoe_Funk_main(n,x,y,pargood));


  return 0;

/*   Covars_g_g_corr_g(n,x,errx,par[0],par[1],covar); */


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
  Amoeba(n,x,y,1,parconf,sigpar,FTOL,MAXITER,Amoe_Funk_confp0);
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
  Amoeba(n,x,y,1,parconf,sigpar,FTOL,MAXITER,Amoe_Funk_confp0);
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
  Amoeba(n,x,y,2,parconf,sigpar,FTOL,MAXITER,Amoe_Funk_confp1);
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
  Amoeba(n,x,y,1,parconf,sigpar,FTOL,MAXITER,Amoe_Funk_confp1);
  if(parconf[0]<par[1])  {
    par1_conf_min=parconf[0];
  }
  else {
    par1_conf_max=parconf[0];
  }
  if(DEBUG) printf("Amoeba  Sigma+sigma(sigma) %14.10f \n",parconf[0]);
/*   cpgpt1(parconf0,parconf[0],3);  */
  if(DEBUG) printf(" Sigma(sigma))= %f sigma(sigma)_teor %14.10f\n",(par1_conf_max-par1_conf_min)/2.,sqrt(covar[1][1]));

  return 0;

  for(i=0;i<NCONFL;i++) {
    parconf[0]=par[1];
    parconf0=par0_conf_min+(par0_conf_max-par0_conf_min)*i/(NCONFL-1.);
    Amoeba(n,x,y,2,parconf,sigpar,FTOL,MAXITER,Amoe_Funk_confp1);
/*     cpgpt1(parconf0,parconf[0],3);  */
 
  }





  amoeba_use=2;
  conflim=0.95;
  parconf[0]=par[0]-5.;
  parconf1=par[1];
  Amoeba(n,x,y,1,parconf,sigpar,FTOL,MAXITER,Amoe_Funk_confp0);
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
  Amoeba(n,x,y,1,parconf,sigpar,FTOL,MAXITER,Amoe_Funk_confp0);
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
    Amoeba(n,x,y,2,parconf,sigpar,FTOL,MAXITER,Amoe_Funk_confp1);
/*     cpgpt1(parconf0,parconf[0],3);  */
 
  }

  amoeba_use=2;
  conflim=0.99;
  parconf[0]=par[0]-5.;
  parconf1=par[1];
  Amoeba(n,x,y,1,parconf,sigpar,FTOL,MAXITER,Amoe_Funk_confp0);
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
  Amoeba(n,x,y,1,parconf,sigpar,FTOL,MAXITER,Amoe_Funk_confp0);
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
    Amoeba(n,x,y,4,parconf,sigpar,FTOL,MAXITER,Amoe_Funk_confp1);
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


float Amoe_Funk_conf(int n, float *x, float *y, float *p) { 

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

float Amoe_Funkg1(int n, float *x, float *y, float *p) { 
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


float Amoe_Funk_main(int n, float *x, float *y, float *p) {
  /* p[0] es h media */
  /* p[1] es sigma_h */

  float pi=4*atan(1.); 
  static int iter=0;
  int i,j;
  double logL=0.,logLb=0.;
  int nstep=1000;
  double conv,xint,xmin,xmax;

  if(p[1]<=0) p[1]=-p[1]; 
  logL=0.;
  logLb=0.;
  for(i=0;i<ndata;i++) {
    conv=0;
    xmin=minf(x[i]-4*sig[i],p[0]-4*p[1]);
    xmax=maxf(x[i]+4*sig[i],p[0]+4*p[1]);
    for(j=0;j<nstep;j++) {
      xint=xmin+(xmax-xmin)*j/(nstep-1);
/*       conv+=gaussian(x[i],xint,sig[i])*gaussian(sig[i],p[2]+xint*xint/p[3],3.)*gaussian(xint,p[0],p[1]);  */
      conv+=gaussian(x[i],xint,sig[i])*gaussian(sig[i],6.+xint*xint/25./25.,1.)*gaussian(xint,p[0],p[1]);     
/*       conv+=gaussian(x[i],xint,sig[i])*gaussian(sig[i],4.,4.)*gaussian(xint,p[0],p[1]);  */
      if(DEBUG2)      printf(" j %d  x %f e %f h %f  conv %g g1 %g  g3 %g\n",j,x[i],sig[i],xint,conv,gaussian(x[i],xint,sig[i]),gaussian(xint,p[0],p[1]));
    }
    conv*=(xmax-xmin)/(nstep-1.);
    logLb-=-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+p[1]*p[1]);
    logLb-=-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+p[1]*p[1]);
    if(conv==0) {
      logL-=-32.;
/*       printf(" DEberia exit(1)\n");  */
      
    }
    else        logL-=log(conv);
    if(DEBUG2) printf("i %4d sig %10f x[i] %10f p[0] %10f p[1] %10f p[2] %10f p[3] %10f anal %10f int %10f diff %10g x1 %f x2 %f\n",i,sig[i],x[i],p[0],p[1],p[2],p[3],-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+p[1]*p[1])-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+p[1]*p[1]),log(conv),-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+p[1]*p[1])-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+p[1]*p[1])-log(conv),xmin,xmax);
  }
  if(DEBUG) printf(" iter %3d logL1 %15f logL2 %15f L1-L2 %15f p0 %f p1 %f p[2] %10f p[3] %10f\n",iter,logL,logLb,logL-logLb,p[0],p[1],p[2],p[3]);
/*   ML[iter]=logL; */
/*   m[iter]=p[0]; */
/*   s[iter]=p[1]; */
  iter++;
  return(logL);
}

float Amoe_Funk_confp0(int n, float *x, float *y, float *p) {
  /* p[0] es h media */
  /* p[1] es sigma_h */

  float pi=4*atan(1.); 
  static int iter=0;
  int i,j;
  double logL=0.,logLb=0.;
  int nstep=1000;
  double conv,xint,xmin,xmax;
  
  logL=0.;
  logLb=0.;
  for(i=0;i<ndata;i++) {
    conv=0;
    xmin=minf(x[i]-4*sig[i],p[0]-4*parconf1);
    xmax=maxf(x[i]+4*sig[i],p[0]+4*parconf1);
    for(j=0;j<nstep;j++) {
      xint=xmin+(xmax-xmin)*j/(nstep-1);
      conv+=gaussian(x[i],xint,sig[i])*gaussian(xint,p[0],parconf1); 
      if(DEBUG2)      printf(" j %d  x %f e %f h %f  conv %g g1 %g  g3 %g\n",j,x[i],sig[i],xint,conv,gaussian(x[i],xint,sig[i]),gaussian(xint,p[0],parconf1));  
    }
    conv*=(xmax-xmin)/(nstep-1.);
    logLb-=-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+parconf1*parconf1); 
    logLb-=-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+parconf1*parconf1); 
    logL-=log(conv);
    if(DEBUG2) printf("i %4d sig %10f x[i] %10f p[0] %10f parconf1 %10f anal %10f int %10f diff %10g x1 %f x2 %f\n",i,sig[i],x[i],p[0],parconf1,-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+parconf1*parconf1)-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+parconf1*parconf1),log(conv),-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+parconf1*parconf1)-(x[i]-p[0])*(x[i]-p[0])/2./(sig[i]*sig[i]+parconf1*parconf1)-log(conv),xmin,xmax);   
  }
  if(DEBUG) printf(" iter %d logL1 %f logL2 %f p0 %f p1 %f\n",iter,logL,logLb,p[0],parconf1);
  return(fabs(logL-(ML[iter-1]-log(conflim))));
}

float Amoe_Funk_confp1(int n, float *x, float *y, float *p) {
  /* p[0] es h media */
  /* p[1] es sigma_h */

  float pi=4*atan(1.); 
  static int iter=0;
  int i,j;
  double logL=0.,logLb=0.;
  int nstep=1000;
  double conv,xint,xmin,xmax;


  if(p[0]<=0) p[0]=-p[0]; 
  logL=0.;
  logLb=0.;
  for(i=0;i<ndata;i++) {
    conv=0;
    xmin=minf(x[i]-4*sig[i],parconf0-4*p[0]);
    xmax=maxf(x[i]+4*sig[i],parconf0+4*p[0]);
    for(j=0;j<nstep;j++) {
      xint=xmin+(xmax-xmin)*j/(nstep-1);
      conv+=gaussian(x[i],xint,sig[i])*gaussian(xint,parconf0,p[0]); 
      if(DEBUG2)      printf(" j %d  x %f e %f h %f  conv %g g1 %g  g3 %g\n",j,x[i],sig[i],xint,conv,gaussian(x[i],xint,sig[i]),gaussian(xint,parconf0,p[0]));  
    }
    conv*=(xmax-xmin)/(nstep-1.);
    logLb-=-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+p[0]*p[0]); 
    logLb-=-(x[i]-parconf0)*(x[i]-parconf0)/2./(sig[i]*sig[i]+p[0]*p[0]); 
    logL-=log(conv);
    if(DEBUG2) printf("i %4d sig %10f x[i] %10f parconf0 %10f p[0] %10f anal %10f int %10f diff %10g x1 %f x2 %f\n",i,sig[i],x[i],parconf0,p[0],-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+p[0]*p[0])-(x[i]-parconf0)*(x[i]-parconf0)/2./(sig[i]*sig[i]+p[0]*p[0]),log(conv),-log(sqrt(2.*pi))-0.5*log(sig[i]*sig[i]+p[0]*p[0])-(x[i]-parconf0)*(x[i]-parconf0)/2./(sig[i]*sig[i]+p[0]*p[0])-log(conv),xmin,xmax);   
  }
  if(DEBUG) printf(" iter %d logL1 %f logL2 %f p0 %f p1 %f\n",iter,logL,logLb,parconf0,p[0]);
  return(fabs(logL-(ML[iter-1]-log(conflim))));
  
}

float Amoe_Funk_g3(int n, float *x, float *y, float *p) {
  /* p[0] es h media */
  /* p[1] es sigma_h */

  float pi=4*atan(1.); 
  int i,j,k;
  double logL=0.,logLb=0.;
  int nstep=100;
  double conv,xint,xmin,xmax;
  double sigint,sigmin,sigmax;

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


void  Covars_g_g_corr_g(int n,float *x,float *errx,float mean,float sigma,float covar[2][2]) {
  
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

