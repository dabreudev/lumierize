#include "modulos.h"
#include <gsl/gsl_machine.h>
#define FTOL  1e-12
#define FTOL2 1e-6
#define FTOL3 1e-7
#define MAXITER  200
#define MAXITER2 70
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
#define TRYEMPIRICAL 0
#endif

/* #define TOLERR 0.0001 */

/* Ahora mismo esta con el pcut definido, de modo que se aleja de valores 
   cercanos a 0. Puede afectar al calculo de las covarianzas!! */

double Amoe_Funk_STY_M2_main(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_M2_norm(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_M2_prior(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_M2_conf(int n, double *x, double *y, double *p);
double Funk2_int_STY_M2(double x);
double Funk1_norm_STY_M2(double x);
void   EmpiricalCovars_STY_M2(int n,double *magn,double *z,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf);

void   ComputeNorma_STY_M2(int n, double mlim, double strrad, double zlow, double zup, struct Schlf_M *lf);

struct cosmo_param *co_STY_M2;
double magl_STY_M2;

int ndata;
int iter_m;
int iter_c;
int nconfl;
double conflim;
/* double ML[5*MAXTRIES*MAXITER]; */
double *pp;
double MLmax;
double *pcut;
double xf,Tf;
double sigi;
double xtmp;
double zlow_STY_M2,zup_STY_M2;


int  MLA_STY_M2(int n,double *magn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf) {

  double *y;
  double par[2];
  double sigpar[2];
  int i;
  int iter_amo;

/*   int j; */
/*   double m1=-31,m2=-13,a1=-2.9,a2=1.0; */
/*   double m,a; */
/*   int nplot=10;         */
/*   double LL[100]; */
/*   double lmin,lmax; */
/*   double cont[10]; */
/*   double tr[6]={0,1,0,0,0,1}; */
/*   char cnul='D'; */


  if(DEBUG) printf(" Estoy auiq\n");

  /* Copias globales de estas variables para que se vean en toda la subrutina */

  ndata=n;
  co_STY_M2=&cosmo;
  magl_STY_M2=mlim;
  zlow_STY_M2=zlow;
  zup_STY_M2=zup;

  y=vector_d(n);
  iter_m=0;

  if(DEBUG3) {
    for(i=0;i<n;i++) printf(" Entrada x %g\n",magn[i]);
  }


  iter_amo=MAXITER+1;
  if(DEBUG3) printf(" iter_amo %d\n",iter_amo);

  
  par[0]=-1.3;  /* Alpha */
  par[1]=-20.4;
  sigpar[0]=.1;
  sigpar[1]=.1;
  
  printf(" Computing LF...\n");

  iter_amo=Amoeba_d(n,magn,z,2,par,sigpar,FTOL,MAXITER,Amoe_Funk_STY_M2_main);  
/*   printf(" FINAL par0 %.15g par1 %.15g \n",par[0],par[1]); */


  MLmax=Amoe_Funk_STY_M2_main(n,magn,z,par);


  /* Meto la solucion en la salida */
 
  lf->alfa=par[0];
  lf->Mstar=par[1];

  /* Estimacion de los errores en los parametros */



  if(TRYEMPIRICAL) {
    if(DEBUG) printf(" llamantry \n");
    conflim=exp(-.5/10.);
    EmpiricalCovars_STY_M2(n,magn,z,par,sigpar,mlim,cosmo,lf); 
    if(DEBUG) printf(" sale \n");

    if(DEBUG) printf(" Solucion final: Mstar %.15g +/- %.15g alpha %.4g +/- %.4g\n",lf->Mstar,lf->errMstar,lf->alfa,lf->erralfa);

  }
  if(CONTOURPLOT) {
/*     nemp_f++; */
/*     conflim=exp(-.5*16.);     */
/*     EmpiricalCovars_STY_M2(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*9.);     */
/*     EmpiricalCovars_STY_M2(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*4.);     */
/*     EmpiricalCovars_STY_M2(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/1.);     */
/*     EmpiricalCovars_STY_M2(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/4.);     */
/*     EmpiricalCovars_STY_M2(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     cpgsci(1); */
/*     cpglab("P\\d1\\u","P\\d3\\u","Contornos de l�mites de confianza"); */
  }
  if(DEBUG) printf(" Calculo empirico\n");

  
/*   conflim=readf(conflim); */

  ComputeNorma_STY_M2(n, mlim,  strrad, zlow,  zup, lf);
  
  free(y);
  if(DEBUG3) printf(" otro mas\n");
  if(DEBUG3) printf(" Este ya\n");
  if(DEBUG3) printf(" el par ta\n");
  if(DEBUG3) printf(" tres mas\n");
/*   free(histo);  */
  if(DEBUG3) printf(" Queda uno \n");

  if(DEBUG3) printf(" He Tu que salgo \n");
  if(DEBUG) printf(" MLF %g\n",MLmax); 

  if(iter_amo>=MAXITER-1) return(2);
  return(0);
}


double Amoe_Funk_STY_M2_main(int n, double *x, double *y, double *p) {

  int i; 
  double logL=0.;
  struct Schlf_M lfamo;
  double Mabs;
  double intsch;
  double Mlow;
  double Lstar;
  double Llow;
/*   double Mup=-30; */
/*   double Lup; */
  double intsup;

  lfamo.alfa=p[0];
  lfamo.phistar=1;
  lfamo.Mstar=p[1];

  logL=0.;

  Lstar=pow(10.,-0.4*lfamo.Mstar);
  intsup=incom(1+lfamo.alfa,100.);
  intsch=Int_sch_M(lfamo,zlow_STY_M2,zup_STY_M2,magl_STY_M2,*co_STY_M2);
  for(i=0;i<ndata;i++) {
    Mabs=Mag(y[i],x[i],*co_STY_M2);
    Mlow=Mag(y[i],magl_STY_M2,*co_STY_M2);
    Llow=pow(10.,-0.4*Mlow);
    if(Llow/Lstar>50) {
      logL+=GSL_LOG_DBL_MAX;
    }
    else {
/*       intsch=(intsup-incom(1+lfamo.alfa,Llow/Lstar));    */
      logL-= log(dVdz(y[i],*co_STY_M2)*Schechter_M(Mabs,lfamo)/intsch); 
    }
/*      printf(" ndata %d  logL %g x %f  y %f Mabs %f Lstar %g Llow %g intsup %g  sch %g int %g\n",i,logL,x[i],y[i],Mabs,Lstar,Llow,intsup,Schechter_M(Mabs,lfamo),intsch);  */
/*     printf(" LF: lfamo.Mstar %g lfamo.phistar %f lf.alfa  %f\n",lfamo.Mstar,lfamo.phistar,lfamo.alfa); */
  }

  iter_m++;
  if(DEBUG) printf(" Iter %d  logL %f par0 %g par1 %g\n",iter_m,logL,p[0],p[1]);
  return(logL);
}


double Amoe_Funk_STY_M2_conf(int n, double *x, double *y, double *p) {

/*   printf(" MLM %g conf %f   %f %f\n",MLmax,conflim,fabs(Amoe_Funk_STY_M2_main(n,x,y,p)-(MLmax-log(conflim))),Amoe_Funk_STY_M2_main(n,x,y,p)); */
   
  return(fabs(Amoe_Funk_STY_M2_main(n,x,y,p)-(MLmax-log(conflim)))); 
}



void   EmpiricalCovars_STY_M2(int n,double *magn,double *z,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf) {

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

  bb=matrix_d(2,1);
  y=vector_d(n);
  parconf=vector_d(2);
  sigparconf=vector_d(2);
  parelip=matrix_d(2,nconflini);
  invcovar=matrix_d(2 ,2 );
  covar=matrix_d(2 ,2 );
  distmax=vector_d(2);

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
    sigparconf[0]=sigpar[0];
    sigparconf[1]=sigpar[1];
    if(i>(int)(nconfl/2.)) {
      parconf[0]=par[0]-((parelip[0])[(int)(i-nconfl/2.)+1]-par[0]);
      parconf[1]=par[1]-((parelip[1])[(int)(i-nconfl/2.)+1]-par[1]);
    }
    iter_c=0;
    if(DEBUG) printf(" antes a nmo\n"); 
    iter_c=Amoeba_d(n,magn,z,2,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_STY_M2_conf);
    if(VERBOSE) printf("CONF %d SOL iter %d par0 %f    par1 %f\n", i,  iter_c ,parconf[0],   parconf[1]);
    (parelip[0])[i]=parconf[0];
    (parelip[1])[i]=parconf[1];
    if(iter_c==0 && Amoe_Funk_STY_M2_conf(n,magn,z,parconf)>FTOL2 ) {
      printf("\b.\b");
      fflush(NULL);
      i--;
    }
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
  if(VERBOSE) printf(" invcovar(0,0) %g invcovar(0,1) %g invcovar(1,0) %g invcovar(1,1) %g\n",invcovar[0][0],invcovar[0][1],invcovar[1][0],invcovar[1][1]);
  for(i=0;i<2;i++) {
    for(j=0;j<2;j++) {
      covar[i][j]=invcovar[i][j];
    }
  }
  gaussj_d(covar,2,bb,1);
  if(VERBOSE) printf("    covar(0,0) %g    covar(0,1) %g    covar(1,0) %g    covar(1,1) %g\n",   covar[0][0],   covar[0][1],   covar[1][0],   covar[1][1]);
  /* Deshago el cambio para el limite de confidencia */
  for(i=0;i<2;i++) {
    for(j=0;j<2;j++) {
      covar[i][j]/=(-2*log(conflim));
    }
  }
  if(VERBOSE) printf("CON covar(0,0) %g    covar(0,1) %g    covar(1,0) %g    covar(1,1) %g\n",   covar[0][0],   covar[0][1],   covar[1][0],   covar[1][1]);

  lf->erralfa=sqrt(covar[0][0]);
  lf->errMstar=sqrt(covar[1][1]);
  lf->covaralfaMstar=covar[0][1];

  free_matrix_d(bb,2,1);
  free(y);
  free(parconf);
  free(sigparconf);

  free_matrix_d(parelip,2,nconflini);
  free_matrix_d(invcovar,2  ,2);
  free_matrix_d(   covar,2  ,2);

  printf("\n");
}





void   ComputeNorma_STY_M2(int n, double mlim, double strrad, double zlow, double zup, struct Schlf_M *lf) {

  double Ntot[NBSNORMA];
  double pi=3.1415926535897932384;
  struct Schlf_M lfbs[NBSNORMA];
  double Ntot_mean, Ntot_sigma;
  int i;

  for(i=0;i<NBSNORMA;i++) {
    lfbs[i].Mstar=lf->Mstar+lf->errMstar*Gasdev();
    lfbs[i].alfa =lf->alfa +lf->erralfa *Gasdev();
    lfbs[i].phistar=1;
    Ntot[i]=Int_sch_M(lfbs[i],zlow,zup,mlim,*co_STY_M2)*strrad/4./pi;
    if(VERBOSE) {
      printf(" %d Ntot %g\n",i,Ntot[i]);
      printf(" Mstar %f  (Mstar %f err %f)\n",lfbs[i].Mstar,lf->Mstar,lf->errMstar);
      printf(" alfa %f  (alfa  %f err %f)\n",lfbs[i].alfa,lf->alfa,lf->erralfa);
    }
  }
  Ntot_mean=StMedia_d(NBSNORMA, Ntot, &Ntot_sigma);

  lf->phistar=1;
  Ntot_mean=Int_sch_M(*lf,zlow,zup,mlim,*co_STY_M2)*strrad/4./pi;

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
