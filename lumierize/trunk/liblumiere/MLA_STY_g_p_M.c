#include "modulos.h"
#include <gsl/gsl_machine.h>
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
#define TRYEMPIRICAL 1
#endif

/* #define TOLERR 0.0001 */

/* Ahora mismo esta con el pcut definido, de modo que se aleja de valores 
   cercanos a 0. Puede afectar al calculo de las covarianzas!! */

double Amoe_Funk_STY_g_p_M_main(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_g_p_M_norm(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_g_p_M_prior(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_g_p_M_conf(int n, double *x, double *y, double *p);
double Funk2_int_STY_g_p_M(double x);
double Funk1_norm_STY_g_p_M(double x);
void   EmpiricalCovars_STY_g_p_M(int n,double *magn,double *errmag,double *z, double *errz,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf);

void   ComputeNorma_STY_g_p_M(int n, double mlim, double strrad, double zlow, double zup, struct Schlf_M *lf);

double *magsig;
double *zsig;
struct cosmo_param *co;
double magl;

int ndata;
int iter_m;
int iter_c;
int nconfl;
double conflim;
/* double ML[5*MAXTRIES*MAXITER]; */
double *pp;
double MLmax;
double xf,Tf;
double sigi;
double xtmp;
double zlow_g,zup_g,strrad_g;


int  MLA_STY_g_p_M(int n,double *magn,double *errmag,double *z, double *errz,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf) {

  double *y;
  double par[3];
  double sigpar[3];
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
  magsig=errmag;
  zsig=errz;
  ndata=n;
  co=&cosmo;
  magl=mlim;
  strrad_g=strrad;
  zlow_g=zlow;
  zup_g=zup;

  y=vector_d(n);
  iter_m=0;


  if(DEBUG3) {
    for(i=0;i<n;i++) printf(" Entrada x %g errx %g\n",magn[i],errmag[i]);
  }


  iter_amo=MAXITER+1;
  if(DEBUG3) printf(" iter_amo %d\n",iter_amo);

  
  par[0]=-0.7;  /* Alpha */
  par[1]=-21.5;
  par[2]=log(0.043);
  sigpar[0]=.1;
  sigpar[1]=.1;
  sigpar[2]=.1;
  
  printf(" Computing LF...\n");

  iter_amo=Amoeba_d(n,magn,z,3,par,sigpar,FTOL,MAXITER,Amoe_Funk_STY_g_p_M_main);  
/*   printf(" FINAL par0 %.15g par1 %.15g \n",par[0],par[1]); */


  MLmax=Amoe_Funk_STY_g_p_M_main(n,magn,z,par);


  /* Meto la solucion en la salida */
 
  lf->alfa=par[0];
  lf->Mstar=par[1];
  lf->phistar=exp(par[2]);

  /* Estimacion de los errores en los parametros */



  if(TRYEMPIRICAL) {
    if(DEBUG) printf(" llamantry \n");
    conflim=exp(-.5/10.);
    EmpiricalCovars_STY_g_p_M(n,magn,errmag,z,errz,par,sigpar,mlim,cosmo,lf); 
    if(DEBUG) printf(" sale \n");

    if(DEBUG) printf(" Solucion final: Mstar %.15g +/- %.15g alpha %.4g +/- %.4g\n",lf->Mstar,lf->errMstar,lf->alfa,lf->erralfa);

  }
  if(CONTOURPLOT) {
/*     nemp_f++; */
/*     conflim=exp(-.5*16.);     */
/*     EmpiricalCovars_STY_g_p_M(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*9.);     */
/*     EmpiricalCovars_STY_g_p_M(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*4.);     */
/*     EmpiricalCovars_STY_g_p_M(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/1.);     */
/*     EmpiricalCovars_STY_g_p_M(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/4.);     */
/*     EmpiricalCovars_STY_g_p_M(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     cpgsci(1); */
/*     cpglab("P\\d1\\u","P\\d3\\u","Contornos de límites de confianza"); */
  }
  if(DEBUG) printf(" Calculo empirico\n");

  
/*   conflim=readf(conflim); */

/*   ComputeNorma_STY_g_p_M(n, mlim,  strrad, zlow,  zup, lf); */
  
  free(y);
  if(DEBUG) printf(" MLF %g\n",MLmax); 

  if(iter_amo>=MAXITER-1) return(2);
  return(0);
}


double Amoe_Funk_STY_g_p_M_main(int n, double *x, double *y, double *p) {

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
  double Ntot;

  lfamo.alfa=p[0];
  lfamo.Mstar=p[1];
  lfamo.phistar=exp(p[2]);

  logL=0.;

  Lstar=pow(10.,-0.4*lfamo.Mstar);
  intsup=incom(1+lfamo.alfa,200.);
  for(i=0;i<ndata;i++) {
    Mabs=Mag(y[i],x[i],*co);
    Mlow=Mag(y[i],magl,*co);
    Llow=pow(10.,-0.4*Mlow);
    if(Llow/Lstar>50) {
      logL+=GSL_LOG_DBL_MAX;
    }
    else {
      intsch=lfamo.phistar*(intsup-incom(1+lfamo.alfa,Llow/Lstar));
      logL-= log(Schechter_M(Mabs,lfamo)/intsch);
    }
/*      printf(" ndata %d  logL %g x %f  y %f Mabs %f Lstar %g Llow %g intsup %g  sch %g int %g\n",i,logL,x[i],y[i],Mabs,Lstar,Llow,intsup,Schechter_M(Mabs,lfamo),intsch);  */
/*     printf(" LF: lfamo.Mstar %g lfamo.phistar %f lf.alfa  %f\n",lfamo.Mstar,lfamo.phistar,lfamo.alfa); */
  }

  if(DEBUG) printf(" logL %g\n",logL);
  /* Aqui viene la parte de la poissoniana de npob */
  Ntot=Int_sch_M(lfamo,zlow_g,zup_g,magl,*co)*strrad_g/4./M_PI;
  logL-= (ndata*log(Ntot) - Ntot - gammln((double)ndata+1.));

  if(DEBUG) printf(" Ntot %g ndata %d rad %g zlw %f zup %f mag %f lfamophi %f\n",Ntot, ndata, strrad_g,zlow_g,zup_g,magl,lfamo.phistar);

  iter_m++;
  if(DEBUG) printf(" Iter %d  logL %f par0 %g par1 %g par2 %g (%.10g)\n",iter_m,logL,p[0],p[1],exp(p[2]),p[2]);
  return(logL);
}


double Amoe_Funk_STY_g_p_M_conf(int n, double *x, double *y, double *p) {

/*   printf(" MLM %g conf %f   %f %f\n",MLmax,conflim,fabs(Amoe_Funk_STY_g_p_M_main(n,x,y,p)-(MLmax-log(conflim))),Amoe_Funk_STY_g_p_M_main(n,x,y,p)); */
   
  return(fabs(Amoe_Funk_STY_g_p_M_main(n,x,y,p)-(MLmax-log(conflim)))); 
}



void   EmpiricalCovars_STY_g_p_M(int n,double *magn,double *errmag,double *z, double *errz,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf) {

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
    iter_c=Amoeba_d(n,magn,z,3,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_STY_g_p_M_conf);
    if(VERBOSE) printf("CONF %d SOL iter %d par0 %f    par1 %f  par2 %f\n", i,  iter_c ,parconf[0],   parconf[1], exp(parconf[2]));
    (parelip[0])[i]=parconf[0];
    (parelip[1])[i]=parconf[1];
    (parelip[2])[i]=parconf[2];
    if(iter_c==0 && Amoe_Funk_STY_g_p_M_conf(n,magn,z,parconf)>FTOL2 ) {
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
      covar[i][j]/=(-2*log(conflim));
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





void   ComputeNorma_STY_g_p_M(int n, double mlim, double strrad, double zlow, double zup, struct Schlf_M *lf) {

  double Ntot[NBSNORMA];
  double pi=3.1415926535897932384;
  struct Schlf_M lfbs[NBSNORMA];
  double Ntot_mean, Ntot_sigma;
  int i;

  for(i=0;i<NBSNORMA;i++) {
    lfbs[i].Mstar=lf->Mstar+lf->errMstar*Gasdev();
    lfbs[i].alfa =lf->alfa +lf->erralfa *Gasdev();
    lfbs[i].phistar=1;
    Ntot[i]=Int_sch_M(lfbs[i],zlow,zup,mlim,*co)*strrad/4./pi;
    if(VERBOSE) {
      printf(" %d Ntot %g\n",i,Ntot[i]);
      printf(" Mstar %f  (Mstar %f err %f)\n",lfbs[i].Mstar,lf->Mstar,lf->errMstar);
      printf(" alfa %f  (alfa  %f err %f)\n",lfbs[i].alfa,lf->alfa,lf->erralfa);
    }
  }
  Ntot_mean=StMedia_d(NBSNORMA, Ntot, &Ntot_sigma);

  lf->phistar=1;
  Ntot_mean=Int_sch_M(*lf,zlow,zup,mlim,*co)*strrad/4./pi;

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
