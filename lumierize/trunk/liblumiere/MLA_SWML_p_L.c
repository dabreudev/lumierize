#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_machine.h>
#include <cpgplot.h>
#include "cosmology.h"
#include "schechter.h"
#include "alloc.h"
#include "mlswml.h"
#include "amoeba.h"

#define FTOL  1e-12
#define FTOL2 1e-6
#define FTOL3 1e-7
#define MAXITER  4000
#define MAXITER2 400
#define MAXTRIES   5
#define DEBUG  1
#define DEBUG2 0
#define DEBUG3 0
#define DEBUG4 0
#define DEBUGPLOT 0
#define CONTOURPLOT 0
#define NCONFL 5
#define TOLERR 0.07 

#ifdef TRYEMPIRICAL
#else 
#define TRYEMPIRICAL 1
#endif

/* #define TOLERR 0.0001 */

double Amoe_Funk_SWML_p_L_main(int n, double *x, double *y, double *p);
double Amoe_Funk_SWML_p_L_norm(int n, double *x, double *y, double *p);
double Amoe_Funk_SWML_p_L_prior(int n, double *x, double *y, double *p);
double Amoe_Funk_SWML_p_L_conf(int n, double *x, double *y, double *p);
double Funk2_int_SWML_p_L(double x);
double Funk1_norm_SWML_p_L(double x);
void   EmpiricalCovars_SWML_p_L(int n,double *flux,double *z,double *par, double *sigpar,double flim, struct cosmo_param cosmo,struct Steplf_L *lf);

void   ComputeNorma_SWML_p_L(int n, double fluxlim, double strrad, double zlow, double zup, struct Steplf_L *lf);

struct cosmo_param *co;
double fluxl;

int ndata;
int iter_m;
int iter_c;
int nconfl_SWML_p_L;
double conflim;
double *pp;
double MLmax;
double xf,Tf;
double sigi;
double xtmp;
double *lumbin;
int nbin;
double zlow_g,zup_g,strrad_g;
int *nullbin;

int  MLA_SWML_p_L(int n,double *flux,double *z,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf) {

  double *y;
  double *par;
  double *sigpar;
  int j;
  int iter_amo;
  int arezeros;
 
  if(DEBUG) printf(" Estoy auiq SWML_p_L\n");

  /* Copias globales de estas variables para que se vean en toda la subrutina */
  ndata=n;
  co=&cosmo;
  fluxl=flim;
  nbin=lf->nbin;
  lumbin=lf->lumi;
  zlow_g=zlow;
  zup_g=zup;
  strrad_g=strrad;

  if(DEBUG) printf(" N Vale %d\n",n);
  if(DEBUG) printf(" nbin %d\n",lf->nbin);

  y=vector_d(n);
  par=vector_d(lf->nbin);
  sigpar=vector_d(lf->nbin);
  nullbin=vector_i(lf->nbin);
  for(j=0;j<nbin;j++) nullbin[j]=0;
  iter_m=0;
  iter_amo=0;


  VVmax_L(n,flux,flux,z,flim,strrad,zlow,zup,cosmo,lf);
  if(DEBUG) for(j=0;j<nbin;j++) printf(" VVMAX Lum %g - %g LF %g\n",lf->lumi[j]/log(10),lf->lumi[j+1]/log(10),lf->lnlf[j]/log(10));
  arezeros=1;
  while(arezeros) {
    arezeros=0;
    for(j=0;j<nbin;j++) {
      if(lf->lnlf[j]==-1/0.) {
	nullbin[j]=1;
	arezeros=1;
	if(j==0) {
	  lf->lnlf[0]=lf->lnlf[1];
	  lf->errlnlf[0]=lf->errlnlf[1];
	}
	else if(j==nbin-1) {
	  lf->lnlf[nbin-1]=lf->lnlf[nbin-2];
	  lf->errlnlf[nbin-1]=lf->errlnlf[nbin-2];
	}
	else {
	  lf->lnlf[j]=(lf->lnlf[j+1]+lf->lnlf[j-1])/2.;
	  lf->errlnlf[j]=(lf->errlnlf[j+1]+lf->errlnlf[j-1])/2.;
	}
      }
    }
  }
  if(DEBUG) for(j=0;j<nbin;j++) printf(" VVMAX Lum %g - %g LF %g\n",lf->lumi[j]/log(10),lf->lumi[j+1]/log(10),lf->lnlf[j]/log(10));
  printf(" Computing LF...\n");
  while(iter_amo==0) 
  { 
    for(j=0;j<nbin;j++) 
    {
      par[j]=lf->lnlf[j];
      sigpar[j]=4*lf->errlnlf[j];
    }
    iter_amo=Amoeba_d(n,flux,z,lf->nbin,par,sigpar,FTOL,MAXITER,Amoe_Funk_SWML_p_L_main);
  }


  MLmax=Amoe_Funk_SWML_p_L_main(n,flux,z,par);


  /* Meto la solucion en la salida */


  for(j=0;j<nbin;j++) lf->lnlf[j]=par[j];  

  for(j=0;j<nbin;j++) {
    if(nullbin[j]) {
      if(j==0) {
	lf->lnlf[0]=lf->lnlf[1];
	lf->errlnlf[0]=lf->errlnlf[1];
      }
      else if(j==nbin-1) {
	lf->lnlf[nbin-1]=lf->lnlf[nbin-2];
	lf->errlnlf[nbin-1]=lf->errlnlf[nbin-2];
      }
      else {
	lf->lnlf[j]=(lf->lnlf[j+1]+lf->lnlf[j-1])/2.;
	lf->errlnlf[j]=(lf->errlnlf[j+1]+lf->errlnlf[j-1])/2.;
      }
    }
  }

  /* Estimacion de los errores en los parametros */



  if(TRYEMPIRICAL) {
    if(DEBUG) printf(" llamantry \n");
    conflim=exp(-.5/10.);
    EmpiricalCovars_SWML_p_L(n,flux,z,par,sigpar,flim,cosmo,lf);  
    if(DEBUG) printf(" sale \n");


  }
  if(CONTOURPLOT) {
/*     nemp_f++; */
/*     conflim=exp(-.5*16.);     */
/*     EmpiricalCovars_SWML_p_L(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*9.);     */
/*     EmpiricalCovars_SWML_p_L(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*4.);     */
/*     EmpiricalCovars_SWML_p_L(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/1.);     */
/*     EmpiricalCovars_SWML_p_L(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/4.);     */
/*     EmpiricalCovars_SWML_p_L(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     cpgsci(1); */
/*     cpglab("P\\d1\\u","P\\d3\\u","Contornos de l�mites de confianza"); */
  }
  if(DEBUG) printf(" Calculo empirico\n");

/*   ComputeNorma_SWML_p_L(n, flim,  strrad, zlow,  zup, lf); */

  
/*   conflim=readf(conflim); */

  free(y);
  free(par);
  free(sigpar);
  free(nullbin);
  return(iter_amo);
}


double Amoe_Funk_SWML_p_L_main(int n, double *x, double *y, double *p) {

  int i,j; 
  int jini = 0;
  double logL=0.;
  double Lumabs;
  double intstep;
  double Llow;
  double funl = 0;
  double Ntot;
  struct Steplf_L lfamo;

/*   for(j=0;j<nbin;j++) printf(" EN %f \n ",p[j]); */
  lfamo.lumi =vector_d(nbin+1);
  lfamo.lnlf   =vector_d(nbin);
  lfamo.nbin =nbin;

  for(j=0;j<nbin;j++) {
    lfamo.lnlf[j]=p[j];
    lfamo.lumi[j]=lumbin[j];
  }
  lfamo.lumi[nbin]=lumbin[nbin];


  if(DEBUG) {
    cpgsci(2);
    cpgsch(2.);
    for(j=0;j<nbin;j++) cpgpt1((float)((lumbin[j+1]+lumbin[j])/2./log(10)),(float)(log10(exp(p[j]))),5);
    if(DEBUG) {
      printf("FL:\n");
      for(j=0;j<nbin;j++) printf("  %f %f %f ",(lumbin[j+1]+lumbin[j])/2.,exp(p[j]),log10(exp(p[j])));
      printf("\n");
    }
  }
  logL=0.;
  for(i=0;i<ndata;i++) {
    Lumabs=Lum(y[i],x[i],*co);
    Llow=Lum(y[i],fluxl,*co);
    intstep=0;
    if(log(Llow)<lumbin[0])  jini=-1;
    if(log(Llow)>lumbin[nbin])   jini=nbin;
    for(j=0;j<nbin;j++) if(log(Llow)>lumbin[j] && log(Llow)<lumbin[j+1]) jini=j;
    if(jini!=-1   && jini!=nbin)   {
      intstep+=exp(p[jini]+lumbin[jini+1])-exp(p[jini]+log(Llow));
      if(DEBUG3) printf(" intstep %f p %g l[ini] %g low %g\n",intstep,p[jini],lumbin[jini],log(Llow));
    }
    for(j=jini+1;j<nbin;j++) 	{
      intstep+=exp(p[j]+lumbin[j+1])-exp(p[j]+lumbin[j]);
      if(DEBUG3) printf(" %d intstep %f p %g l[j+1] %g l[j] %g\n",j,intstep,p[j],lumbin[j+1],lumbin[j]);
    }
    if(log(Lumabs)>lumbin[nbin] || log(Lumabs)<lumbin[0]);
    else {
      for(j=0;j<nbin;j++) {
	if(DEBUG3) printf(" Comp j %d  %f  %f  %f \n",j,lumbin[j],log(Lumabs),lumbin[j+1]);
	if(lumbin[j]<log(Lumabs) && log(Lumabs)<lumbin[j+1]) {
	  funl=exp(p[j]);
	  if(DEBUG3) printf(" Es sta %f %g    %f %f %f \n",p[j],funl,lumbin[j],log(Lumabs),lumbin[j+1]);
	  break;
	}
      }
      logL-= log(funl/intstep); 
/*       logL-= log(funl); */
    }
    if(DEBUG3) printf(" ndata %d  logL %g x %f  y %f Lumabs %g  Llow %g funl %g int %g\n",i,logL,log10(x[i]),y[i],Lumabs,Llow,funl,intstep);
  }

  /* Aqui viene la parte de la poissoniana de npob */
  Ntot=Int_step_L(lfamo,zlow_g,zup_g,fluxl,*co)*strrad_g/4./M_PI;
  logL-= (ndata*log(Ntot) - Ntot - gammln((double)ndata+1.));


  if(DEBUG) {
    cpgsci(0);
    cpgsch(2.);
    for(j=0;j<nbin;j++) cpgpt1((float)((lumbin[j+1]+lumbin[j])/2./log(10)),(float)(log10(exp(p[j]))),5);
  }
  if(DEBUG)  printf(" iter %d logL %f\n",iter_m,logL);
  if(DEBUG2) {
    printf("FL:\n");
    for(j=0;j<nbin;j++) printf("  L %f  LF  %f ",(lumbin[j+1]+lumbin[j])/2./log(10),log10(exp(p[j])));
    printf("\n");
  }

  free(lfamo.lumi);
  free(lfamo.lnlf);

  iter_m++;
  return(logL);
}


double Amoe_Funk_SWML_p_L_conf(int n, double *x, double *y, double *p) {

/*   printf(" MLM %g conf %f   %f %f\n",MLmax,conflim,fabs(Amoe_Funk_SWML_p_L_main(n,x,y,p)-(MLmax-log(conflim))),Amoe_Funk_SWML_p_L_main(n,x,y,p)); */
   
  return(fabs(Amoe_Funk_SWML_p_L_main(n,x,y,p)-(MLmax-log(conflim)))); 
}



void   EmpiricalCovars_SWML_p_L(int n,double *flux,double *z, double *par, double *sigpar,double flim, struct cosmo_param cosmo,struct Steplf_L *lf) {


  int i,j;  
  double *parconf; 
  double *sigparconf; 
  double **invcovar;
  double **covar;
  double **parelip;
  double *y;

  double **bb;

  


  double first, median, third, *distmax;

  int nconfl,nconflini;

  if(DEBUG) printf(" n vale %d \n",n);
  nconfl=NCONFL*lf->nbin;
  nconfl=2*lf->nbin*(lf->nbin+1)/2;    /* Esto es el numero de elementos de la matriz de covar multiplicado por 2 */
  nconflini=nconfl;

  if(DEBUG4) printf(" nconfl %d\n",nconfl);

  bb=matrix_d(lf->nbin,1);
  printf(" bb %p\n",bb);
  for(i=0;i<lf->nbin;i++) (bb[i])[0]=0;
  y=vector_d(n);
  parconf=vector_d(lf->nbin);
  sigparconf=vector_d(lf->nbin);
  parelip=matrix_d(lf->nbin,  nconflini);

  invcovar=matrix_d(lf->nbin  ,lf->nbin  );
  covar   =matrix_d(lf->nbin  ,lf->nbin  );


  distmax =vector_d(lf->nbin);


  printf(" Error step  ");
  for(i=0;i<nconfl;i++) printf(".");
  for(i=0;i<nconfl;i++) printf("\b");
  fflush(NULL);

  for(i=0;i<nconfl;i++) {
    printf("#");
    fflush(NULL);
    for(j=0;j<lf->nbin;j++) {
      parconf[j]=par[j]+sigpar[j]*Gasdev();
      sigparconf[j]=sigpar[j];
    }
    if(i>(int)(nconfl/2.)) {
      for(j=0;j<lf->nbin;j++) {
	parconf[j]=par[j]-((parelip[j])[(int)(i-nconfl/2.)+1]-par[j]);
	sigparconf[j]=((parelip[j])[(int)(i-nconfl/2.)+1]-par[j])/2.; 
      }
    }
    iter_c=0;
    if(DEBUG) printf(" antes a nmo\n"); 
    iter_c=Amoeba_d(n,flux,z,lf->nbin,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_SWML_p_L_conf);
    if(DEBUG)  {
      printf(" SOL %d/%diter %d\n",i,nconfl,iter_c);
      for(j=0;j<nbin;j++) printf(" SOL par%d %.15g\n",j,parconf[j]);
    }
    for(j=0;j<lf->nbin;j++)  (parelip[j])[i]=parconf[j];
    if(DEBUG) printf(" Despues parelip =\n");
    if(iter_c==0 && Amoe_Funk_SWML_p_L_conf(n,flux,z,parconf)>FTOL2 ) {
      printf("\b.\b");
      i--;
      fflush(NULL);
    }
    if(DEBUG) printf(" Esto es despues check\n");

    /* Si hay algun bin vacio, le doy los valores de los de al lado */
    for(j=0;j<lf->nbin;j++) {
      if(nullbin[j]) {
	if(j==0) parconf[0]=parconf[1];
	else if(j==nbin-1)  parconf[nbin-1]=parconf[nbin-2];
	else                parconf[j]=(parconf[j+1]+parconf[j-1])/2.;
      }
    }
  }

  if(DEBUG) printf(" Esto es despues de todo\n");

  /* Supongo que el centro de la elipse es el valor que maximiza ML */
  for(i=0;i<nconfl;i++) {
    for(j=0;j<lf->nbin;j++)       (parelip[j])[i]-=par[j];
  }


  if(DEBUG4) {
    printf(" Las soluciones primeras son : \n");
    for(i=0;i<nconfl;i++) {
      printf(" CONFINI %d ",i);
      for(j=0;j<lf->nbin;j++)  printf(" %f ",(parelip[j])[i]);
      printf(" \n");
    }
  }
  
  /* Detecto puntos que esten muy alejados de la supuesta elipse */
  for(j=0;j<lf->nbin;j++) {
    Quartil_d(nconfl,parelip[j],&first,&median,&third);
    distmax[j]=maxf(fabs(first),fabs(third));
    if(DEBUG4) printf(" Los quartiles de %d son: fir %f med %f th %f dist %f\n",j,first,median,third,distmax[j]);
  }
  for(i=0;i<nconfl;i++) {
    for(j=0;j<lf->nbin;j++) {
      if(fabs((parelip[j])[i])>4*2*distmax[j]/1.35) {
	if(DEBUG4) printf(" Elimino %d por el %d : %f \n",i,j,(parelip[j])[i]);
	for(j=0;j<lf->nbin;j++) {
	  memcpy(parelip[j]+i,parelip[j]+i+1,(nconfl-i-1)*sizeof(double));
	}
	i--;
	nconfl--;
	break;
      }
    }
  }

  if(DEBUG4) {
    printf(" Uso estas soluciones: \n");
    for(i=0;i<nconfl;i++) {
      printf(" CONF %d ",i);
      for(j=0;j<lf->nbin;j++)  printf(" %f ",(parelip[j])[i]);
      printf(" \n");
    }
  }
	    




  MCElipN_d(nconfl,lf->nbin,parelip,invcovar);

  if(DEBUG4) printf(" Ya he calculado las elip\n");
  


  for(i=0;i<lf->nbin;i++) {
    for(j=0;j<lf->nbin;j++) {
      covar[i][j]=invcovar[i][j];
      if(DEBUG4) printf(" cov %d %d %g\n",i,j,covar[i][j]);
      printf(" cov %d %d %g\n",i,j,covar[i][j]);
    }
  } 

  printf(" Esto es justo antes de gauus bb %p\n",bb);

  gaussj_d(covar,lf->nbin,bb,1);

  printf(" Y esto depsues\n");
  
  /* Deshago el cambio para el limite de confidencia */
  for(i=0;i<lf->nbin;i++) {
    for(j=0;j<lf->nbin;j++) {
      covar[i][j]/=(-2*log(conflim));
      lf->covarlnlf[i][j]=covar[i][j];
      if(DEBUG4) printf(" invcov %d %d %g\n",i,j,covar[i][j]);
   }
  }

  /* Pongo las covar a 0 de los que tienen un bin vacio */
  for(j=0;j<lf->nbin;j++) {
    if(nullbin[j]) {
      lf->lnlf[j]=-1/0.;
      for(i=0;i<lf->nbin;i++) {
	lf->covarlnlf[i][j]=0;
	lf->covarlnlf[j][i]=0;
      }
    }
  }



  printf(" Aqui estoy acabando\n");

  for(i=0;i<lf->nbin;i++) 
    lf->errlnlf[i]=sqrt(lf->covarlnlf[i][i]);

  free_matrix_d(bb,lf->nbin,1);
  free(y);
  free(sigparconf);
  free(parconf);
  free(distmax);

  free_matrix_d(parelip,lf->nbin  ,nconflini);

  free_matrix_d(invcovar,lf->nbin  ,lf->nbin);
  free_matrix_d(   covar,lf->nbin  ,lf->nbin);


  printf("\n");
}





void   ComputeNorma_SWML_p_L(int n, double flim, double strrad, double zlow, double zup, struct Steplf_L *lf) {

  double Ntot;
  double pi=3.1415926535897932384;
  int j;

  Ntot=Int_step_L(*lf,zlow,zup,flim,*co)*strrad/4./pi; 
  

  for(j=0;j<lf->nbin;j++) {
    if(DEBUG) printf(" From %f ",lf->lnlf[j]);
    lf->lnlf[j]=lf->lnlf[j]+log((float)(n)/Ntot);
    if(DEBUG) printf(" to %f \n",lf->lnlf[j]);
    /* El error, como es logaritmico, sigue intacto */
    lf->errlnlf[j]=lf->errlnlf[j];
  }

}
