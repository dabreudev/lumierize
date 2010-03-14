#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_machine.h>
#include "alloc.h"
#include "mlhist.h"
#include "sthisto.h"

#define FTOL  1e-12
#define FTOL2 1e-6
#define FTOL3 1e-7
#define MAXITER  2000
#define MAXITER2 200
#define MAXITER3 150
#define MAXTRIES   5
#define DEBUG  1
#define DEBUG2 0
#define DEBUG3 0
#define DEBUGPLOT 0
#define CONTOURPLOT 0
#define NCONFL 10
#define TOLERR 0.07 

#ifdef TRYEMPIRICAL
#else 
#define TRYEMPIRICAL 1
#endif

#define USEPCUT 1
/* #define TOLERR 0.0001 */

/* Ahora mismo esta con el pcut definido, de modo que se aleja de valores 
   cercanos a 0. Puede afectar al calculo de las covarianzas!! */

double Amoe_Funk_h_g_f_d_main(int n, double *x, double *y, double *p);
double Amoe_Funk_h_g_f_d_norm(int n, double *x, double *y, double *p);
double Amoe_Funk_h_g_f_d_prior(int n, double *x, double *y, double *p);
double Amoe_Funk_h_g_f_d_conf(int n, double *x, double *y, double *p);
double Funk2_int_h_g_f_d(double x);
double Funk1_norm_h_g_f_d(double x);
void   EmpiricalCovars_h_g_f_d(int n,double *x,double *errx,int k,double *xk,double *Pk,double *sigPk,double **covar);


double *sig;
double *xklim;
int ndata;
int iter_m;
int iter_c;
int kpar;
int nconfl;
double conflim;
/* double ML[5*MAXTRIES*MAXITER]; */
double *pp;
double MLmax;
double *pcut;
double xf,Tf;
double sigi;
double xtmp;

int nemp_f=2;

int  MLA_h_g_f_d(int n,double *x,double *errx,double xfermi, double Tfermi, int k, double *xk, double *Pk, double **covPk) {

  double *y;
  double *par;
  double *sigpar;
  double norm;
  int *histo;
  int i,j;
  double **covar; 
  int iter_amo;
  int tries=0;

  /* Copias globales de estas variables para que se vean en toda la subrutina */
  sig=errx;
  ndata=n;
  kpar=k;
  xklim=xk;
  xf=xfermi;
  Tf=Tfermi;

  printf(" n= %d k =%d \n",n,k);

  y=vector_d(n);
  par=vector_d(k+1);
  sigpar=vector_d(k+1);
  pcut=vector_d(k);
  covar=matrix_d(k,k);

  iter_m=0;

  if(DEBUG3) {
    for(i=0;i<n;i++) printf(" Entrada x %g errx %g\n",x[i],errx[i]);
  }


  iter_amo=MAXITER+1;
  if(DEBUG3) printf(" iter_amo %d\n",iter_amo);

  if(DEBUG3) {
    printf(" Y aqhira\n");
    histo=StHisto2_d(n,x,k,xk,xk+k);
    printf(" YESOA OO A xk %f xk+n %f\n",xk[0],xk[k]);
    free(histo);
    printf(" YAA \n");
  }
  
  
  while(iter_amo>= MAXITER) {
    if((histo=StHisto2_d(n,x,k,xk,xk+k))==NULL) { printf("I cannot dimension histo of %d elements \n",k);exit(1);}

    for(j=0;j<kpar;j++) {
      if(histo[j]==0)  par[j]=   1.           /(xk[j+1]-xk[j]);
      else             par[j]=(double)histo[j]/(xk[j+1]-xk[j]);
    }
    norm=0;
    

    for(j=0;j<kpar;j++) norm+=par[j]*(xk[j+1]-xk[j]);
    for(j=0;j<kpar;j++) { 
      par[j]=par[j]/norm;
      if(histo[j]==0) sigpar[j]=sqrt(par[j]*(xk[j+1]-xk[j])*(1-par[j]*(xk[j+1]-xk[j]))/1.);
      else            sigpar[j]=sqrt(par[j]*(xk[j+1]-xk[j])*(1-par[j]*(xk[j+1]-xk[j]))/histo[j]);
      if(par[j]/sigpar[j]>1) pcut[j]=par[j]/sigpar[j]/10000./(xk[j+1]-xk[j]);
      else                   pcut[j]=       1.       /10000./(xk[j+1]-xk[j]);
      if(DEBUG2) printf(" INI %d par %f sig %f histo %d\n",j,par[j],sigpar[j],histo[j]);
    }
    
    free(histo);

    iter_amo=Amoeba_d(n,x,y,k,par,sigpar,FTOL,MAXITER,Amoe_Funk_h_g_f_d_main); 
    if(DEBUG2) printf("\n INTERMEDIO par0 %.15g par1 %.15g \n\n",par[0],par[1]);
    iter_amo=MAXITER+1;   
    tries++;
    
    while(iter_amo>= MAXITER ) { 
      iter_amo=Amoeba_d(n,x,y,k,par,sigpar,FTOL,MAXITER,Amoe_Funk_h_g_f_d_main); 
      tries++;
      if(DEBUG) printf("\n BUSCANDOINTER par0 %.15g par1 %.15g \n\n",par[0],par[1]);
      for(j=0;j<kpar;j++) sigpar[j]=sqrt(fabs(par[j]));
    }
	


    if(DEBUG2) printf(" FINAL par0 %.15g par1 %.15g \n",par[0],par[1]);
    if(DEBUG2) printf(" FINAL par0 %.15g par1 %.15g par2 %.15g\n",par[0],par[1],par[2]);

/*     free(histo);  */
  }


  /* Aunque la hemos impuesto dentro de Anoeb_funk, volvemos a normalizar */
  norm=0;
  for(j=0;j<kpar;j++) norm+=par[j]*(xk[j+1]-xk[j]);
  for(j=0;j<kpar;j++) par[j]=par[j]/norm;

  MLmax=Amoe_Funk_h_g_f_d_main(n,x,y,par);


  /* Meto la solucion en Pk */
  for(j=0;j<kpar;j++) Pk[j]=par[j];


  /* Estimacion de los errores en los parametros */


  if(DEBUGPLOT) {
    cpgpage();  
    cpgswin(Pk[0]-0.6,Pk[0]+0.6,Pk[1]-0.6,Pk[1]+0.6);
    cpgbox("BCTNS",0,0,"BCTNS",0,0);
  }

  if(TRYEMPIRICAL) {
    conflim=exp(-.5/10.);
    EmpiricalCovars_h_g_f_d(n,x,errx,k,xk,Pk,sigpar,covar); 
  }
  if(CONTOURPLOT) {
    nemp_f++;
    conflim=exp(-.5*16.);    
    EmpiricalCovars_h_g_f_d(n,x,errx,k,xk,Pk,sigpar,covar); 
    nemp_f++;
    conflim=exp(-.5*9.);    
    EmpiricalCovars_h_g_f_d(n,x,errx,k,xk,Pk,sigpar,covar); 
    nemp_f++;
    conflim=exp(-.5*4.);    
    EmpiricalCovars_h_g_f_d(n,x,errx,k,xk,Pk,sigpar,covar); 
    nemp_f++;
    conflim=exp(-.5/1.);    
    EmpiricalCovars_h_g_f_d(n,x,errx,k,xk,Pk,sigpar,covar); 
    nemp_f++;
    conflim=exp(-.5/4.);    
    EmpiricalCovars_h_g_f_d(n,x,errx,k,xk,Pk,sigpar,covar); 
    nemp_f++;
    cpgsci(1);
    cpglab("P\\d1\\u","P\\d3\\u","Contornos de l�mites de confianza");
  }
  if(DEBUG) printf(" Calculo empirico\n");
  for(i=0;i<k;i++) {
    for(j=0;j<k;j++) {
      covPk[i][j]=covar[i][j];
      if(DEBUG2) printf(" covar %d %d = %g\n",i,j,covar[i][j]);
    }
  }
  
/*   conflim=readf(conflim); */

  
  free(y);
  if(DEBUG3) printf(" otro mas\n");
  free_matrix_d(covar,k,k); 
  if(DEBUG3) printf(" Este ya\n");
  free(par);
  if(DEBUG3) printf(" el par ta\n");
  free(sigpar);
  if(DEBUG3) printf(" tres mas\n");
/*   free(histo);  */
  if(DEBUG3) printf(" Queda uno \n");
  free(pcut); 

  if(DEBUG3) printf(" He Tu que salgo \n");
  if(DEBUG) printf(" MLF %g\n",MLmax); 

  if(iter_amo>=MAXITER-1) return(2);
  return(0);
}


double Amoe_Funk_h_g_f_d_main(int n, double *x, double *y, double *p) {

  int i,j;
  double logL=0.;
  double ltmp;
  double tmp1,priori;
  double norm=0;
  double scale,offset;
  double distnorm;
  double meanx,stdx;
  int    nstepnorm=20;
  int    nstep=20;
  (void)y;/* To avoid warning */
  (void)n;/* To avoid warning */


  logL=0.;

  /* Imponemos la normalizacion de la funcion de probabilidad */
  norm=0;
  for(j=0;j<kpar;j++) {
    if(DEBUG3) printf(" Antes p %d %f \n",j,p[j]);
    if(p[j]<0) p[j]=-p[j]; 
/*     if(p[j]==0) p[j]=1./(xklim[j+1]-xklim[j]); */
    norm+=p[j]*(xklim[j+1]-xklim[j]);
  }
  for(j=0;j<kpar;j++) p[j]=p[j]/norm;

  meanx=StWeightMedia_d(kpar,xklim,p,&stdx);

  pp=p;

  for(i=0;i<ndata;i++) {
    sigi=sig[i];
    xtmp=x[i];
    offset=meanx;
    scale=sqrt(sig[i]*sig[i]+stdx*stdx);
    distnorm=gaussinther_d(Funk1_norm_h_g_f_d,offset,scale,nstepnorm); 
/*     printf(" distnorm %f \n",distnorm);  */
      
    ltmp=0;
    priori=1;
    for(j=0;j<kpar;j++) {
      xtmp=x[i];
      priori*=exp(-pcut[j]/p[j]);  /* Distribucion a priori. Para evitar valores iguales a 0 */
      tmp1=p[j]*gaussintleg_d(Funk2_int_h_g_f_d,xklim[j],xklim[j+1],nstep);  
      ltmp+=tmp1;
      if(DEBUG3) printf("j %d ltmp %g tmp1 %g aqui %f p[j] %f pcut %f priori %f erf1 %f erf2 %f   x %f xk %f xk+1 %f arg1 %g arg2 %f\n",j,ltmp,tmp1,-p[j]/2.*(+erf((x[i]-xklim[j+1])/sig[i]/sqrt(2.))-erf((x[i]-xklim[j])/sig[i]/sqrt(2.))),p[j],pcut[j],priori,erf((x[i]-xklim[j+1])/sig[i]),erf((x[i]-xklim[j])/sig[i]),x[i],xklim[j],xklim[j+1],(x[i]-xklim[j+1])/sig[i],(x[i]-xklim[j])/sig[i]);
    }
    ltmp=ltmp/distnorm;
    if(DEBUG3) printf("i %d logL %20f ltmp %f \n",i,logL,ltmp);
    if(!USEPCUT)  priori=1;  
    if(ltmp!=0)  logL-=log(ltmp*priori);
    else         logL-=log(FTOL);
    if(DEBUG3) printf(" logL despues %20f\n",logL);
  }
  if(DEBUG) printf(" MLF iter %3d logL %.15g ",iter_m,logL);
  if(DEBUG) {
    printf(" Llamado con ");
    for(j=0;j<kpar;j++) printf(" %5.10f",p[j]);
    printf("\n");
  }
  
/*   Amoe_Funk_h_g_d_main(n,x,y,p); */
  
  /*   ML[iter_m]=logL;  */
  iter_m++;
  return(logL);
}


double Amoe_Funk_h_g_f_d_conf(int n, double *x, double *y, double *p) {
   
  int i,j;
  double logL=0.;
  double ltmp;
  double tmp1,priori;
  double norm=0;
  (void)y;/* To avoid warning */
  (void)n;/* To avoid warning */


  logL=0.;

  /* Imponemos la normalizacion de la funcion de probabilidad */
  norm=0;
  for(j=0;j<kpar;j++) {
    if(p[j]<0) p[j]=-p[j];
    norm+=p[j]*(xklim[j+1]-xklim[j]);
  }
  for(j=0;j<kpar;j++) p[j]=p[j]/norm;

  for(i=0;i<ndata;i++) {
    ltmp=0;
    priori=1;
    for(j=0;j<kpar;j++) {
      priori*=exp(-pcut[j]/p[j]);  /* Distribucion a priori. Para evitar valores iguales a 0 */
      tmp1=p[j]*(erf((x[i]-xklim[j+1])/sig[i]/sqrt(2.))-erf((x[i]-xklim[j])/sig[i]/sqrt(2.)));
      ltmp+=tmp1;
    }
    if(!USEPCUT)  priori=1;  
    if(ltmp!=0)  logL-=log(ltmp*priori);
    else         logL-=log(FTOL);
  }
  if(DEBUG3) printf(" MLFCONF iter %3d logL-logL0+log(conf) %.15g MLmin %.15g log(conf) %f\n ",iter_m,fabs(logL-(MLmax-log(conflim))),MLmax,log(conflim));
  if(DEBUG3) {
    printf(" Llamado con ");
    for(j=0;j<kpar;j++) printf(" %5.10f",p[j]);
    printf("\n");
  }

  if(DEBUG3) cpgpt1(p[0],p[1],1);

  return(fabs(logL-(MLmax-log(conflim))));
}



void   EmpiricalCovars_h_g_f_d(int n,double *x,double *errx,int k,double *xk,double *Pk,double *sigPk, double **covar) {

  int i,j,l;
  double *y;
  double **bb;
  double *parconf;
  double *sigparconf;
  double **invcovar;
  double **parelip;

  double **pareA;
  double **pareB;
  double **pareC;
  
  double **invcovA;
  double **invcovB;
  double **invcovC;
  double **covA;
  double **covB;
  double **covC;
  
  int nconflini;
  double first, median, third, *distmax;

/*   double a,b,c,d,f,e;  */
/*   double xcelip,ycelip,elipa,elipb,elipt;  */

/*   float xtmp,ytmp; */
  
  nconfl=NCONFL*kpar;
  nconflini=nconfl;

  bb=matrix_d(kpar-1,1);
  y=vector_d(n);
  parconf=vector_d(k);
  sigparconf=vector_d(k);

  parelip=matrix_d(k  ,nconflini);
  pareA  =matrix_d(k-1,nconflini);
  pareB  =matrix_d(k-1,nconflini);
  pareC  =matrix_d(k-1,nconflini);

  invcovar=matrix_d(kpar  ,kpar  );
  invcovA =matrix_d(kpar-1,kpar-1);
  invcovB =matrix_d(kpar-1,kpar-1);
  invcovC =matrix_d(kpar-1,kpar-1);
  covA    =matrix_d(kpar-1,kpar-1);
  covB    =matrix_d(kpar-1,kpar-1);
  covC    =matrix_d(kpar-1,kpar-1);

  distmax=vector_d(2);


/*   conflim=exp(-.5/100.);     */         //Los puntos corresponderan a 1 sigma entre 4 de desviacion para dist. normal en par
/*   conflim=exp(-.5);   */           //Los puntos corresponderan a 1 sigma entre 4 de desviacion para dist. normal en par
/*   conflim=exp(-.5/100.);    */


  if(DEBUG3) printf(" log(conflim) %f \n",log(conflim));
  
  for(i=0;i<nconfl;i++) {
    for(j=0;j<kpar;j++) {
      parconf[j]=Pk[j]+sigPk[j]*3.*Gasdev(); 
      sigparconf[j]=sigPk[j]; 
    }
    if(DEBUG3) printf(" INIPARTEST %d %d\n",i,(int)(nconfl/2.));
    if(i>(int)(nconfl/2.)) {
      for(j=0;j<kpar;j++) {
	parconf[j]=Pk[j]-((parelip[j])[(int)(i-nconfl/2.)+1]-Pk[j]);
	sigparconf[j]=((parelip[j])[(int)(i-nconfl/2.)+1]-Pk[j])/2.; 
      }
      if(DEBUG3) {
	cpgsci(3);
	cpgpt1((parconf[0]),(parconf[1]),5);
	cpgsci(1);
      }
    }

    if(DEBUG3) {
      cpgsci(2);
      cpgpt1((parconf[0]),(parconf[1]),4);
      cpgsci(1);
    }

    iter_c=0;
    iter_c=Amoeba_d(n,x,y,k,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_h_g_f_d_conf);

    if(DEBUG3) printf(" PARCONF %g %g iter_c %d AMF %g i %d\n",parconf[0],parconf[1],iter_c,Amoe_Funk_h_g_f_d_conf(n,x,y,parconf),i);

    for(j=0;j<kpar;j++)     (parelip[j])[i]=parconf[j];

    if(DEBUG3) {
      printf(" PARCONFEL ");
      for(j=0;j<i;j++)  {
	printf(" ## ");
	for(l=0;l<kpar;l++)    printf(" %g ",(parelip[l])[j]);
      }
      printf("\n");
    }
      
    if(iter_c==0 && Amoe_Funk_h_g_f_d_conf(n,x,y,parconf)>FTOL2 ) i--;
  }
  
  
  
  if(DEBUG) {
    printf(" PARCONFSALIDA ");
    for(j=0;j<nconfl;j++)  {
      printf(" %d ",j);
      for(l=0;l<kpar;l++)    printf(" %g ",(parelip[l])[j]-Pk[l]);
      printf(" \n");
      cpgsci(7);
      cpgpt1((parelip[0])[j],(parelip[1])[j],6);
    }
    printf("\n");
  }
/*   i=readi(i); */
  
  /* Supongo que el centro de la elipse es el valor que maximiza ML */

   for(i=0;i<nconfl;i++)    for(j=0;j<kpar;j++)       (parelip[j])[i]-=Pk[j];


  if(DEBUGPLOT) {
    for(i=0;i<nconfl;i++) {
      if(DEBUG3) printf(" PARCONFREAL %g %g\n",(parelip[0])[i],(parelip[1])[i]);
      cpgsci(nemp_f);
      cpgpt1((parelip[0])[i]+Pk[0],(parelip[1])[i]+Pk[1],1);
    }
  }

  /*   kpar=readi(kpar); */
  
  /* Detecto puntos que esten muy alejados de la supuesta elipse */
  for(j=0;j<kpar;j++) {
    Quartil_d(nconfl,parelip[j],&first,&median,&third);
    distmax[j]=maxf(fabs(first),fabs(third));
  }  
  for(i=0;i<nconfl;i++) {
    for(j=0;j<kpar;j++) {
      if(fabs((parelip[j])[i])>2*2*distmax[j]/1.35) {
	for(j=0;j<kpar;j++) {
	  memcpy(parelip[j]+i,parelip[j]+i+1,(nconfl-i-1)*sizeof(double));
	}
	i--;
	nconfl--;
	break;
      }
    }
  }
 
  /* Ajusto elipses en tres subespacios diferentes. De modo que al final obtengo 
     la matriz entera. Hay que tener en cuenta que esta matriz del hessiano tiene 
     determinante nulo y por lo tanto las superficies de sigma constante
     son formas cuadraticas de dimension kpar-1 en un espacio de dim kpar. Esto 
     es asi por la ligadura de que la normalizacion de las Pk*/

  for(i=0;i<nconfl;i++) {
    for(j=0;j<kpar-1;j++) {
      pareA[j][i]=parelip[j][i];
    }
  }
  MCElipN_d(nconfl,kpar-1,pareA,invcovA);
  for(i=0;i<nconfl;i++) {
    for(j=0;j<kpar-1;j++) {
      pareB[j][i]=parelip[j+1][i];
    }
  }
  MCElipN_d(nconfl,kpar-1,pareB,invcovB);
  for(i=0;i<nconfl;i++) {
    for(j=0;j<kpar-2;j++) {
      pareC[j][i]=parelip[j][i];
    }
    pareC[kpar-2][i]=parelip[kpar-1][i];
  }
  MCElipN_d(nconfl,kpar-1,pareC,invcovC);

  if(DEBUG) {
    printf(" ICOVA\n");
    for(i=0;i<k-1;i++) {
      for(j=0;j<k-1;j++) {
	printf("  %g ",invcovA[i][j]);  
      }
      printf("\n");
    }
    printf(" ICOVB\n");
    for(i=0;i<k-1;i++) {
      for(j=0;j<k-1;j++) {
	printf("  %g ",invcovB[i][j]);  
      }
      printf("\n");
    }
    printf(" ICOVC\n");
    for(i=0;i<k-1;i++) {
      for(j=0;j<k-1;j++) {
	printf("  %g ",invcovC[i][j]);  
      }
      printf("\n");
    }
  }
  for(i=0;i<kpar-1;i++) {
    for(j=0;j<kpar-1;j++) {
      covA[i][j]=invcovA[i][j];
      covB[i][j]=invcovB[i][j];
      covC[i][j]=invcovC[i][j];
    }
  } 
  gaussj_d(covA,kpar-1,bb,1);
  gaussj_d(covB,kpar-1,bb,1);
  gaussj_d(covC,kpar-1,bb,1);
  if(DEBUG2) {
    printf(" COVA\n");
    for(i=0;i<k-1;i++) {
      for(j=0;j<k-1;j++) {
	printf("   %g  ",covA[i][j]);  
      }
      printf("\n");
    }
    printf(" COVB\n");
    for(i=0;i<k-1;i++) {
      for(j=0;j<k-1;j++) {
	printf("   %g  ",covB[i][j]);  
      }
      printf("\n");
    }
    printf(" COVC\n");
    for(i=0;i<k-1;i++) {
      for(j=0;j<k-1;j++) {
	printf("   %g  ",covC[i][j]);  
      }
      printf("\n");
    }
  }

  /* Relleno la matriz de covarianza con las tres auxiliares */
  /* Uso covA para casi todo */
  for(i=0;i<kpar-1;i++) {
    for(j=0;j<kpar-1;j++) {
      covar[i][j]=covA[i][j];
    }
  } 

  if(DEBUG3) printf(" primer paso \n");
  /* Uso covB para la fila de abajo, la columna de la derecha y el extremo inferior derecha */

  for(j=1;j<kpar;j++)       covar[kpar-1][j]=covB[kpar-2][j-1];
  for(i=1;i<kpar;i++)       covar[i][kpar-1]=covB[i-1][kpar-2];
  if(DEBUG3) printf(" segund paso \n");
  /* Uso covC para el extremo superior derecha y el inferior izquierda  */
  covar[0][kpar-1]=covC[0][kpar-2];
  covar[kpar-1][0]=covC[kpar-2][0];
  if(DEBUG3) printf(" tercer paso \n");

  if(DEBUG) {
    for(i=0;i<k;i++) {
      for(j=0;j<k;j++) {
	printf(" covar  %01d %01d = %-15.12g",i,j,covar[i][j]);
      }
      printf("\n");
    }
  }
  



  /* Deshago el cambio para el limite de confidencia */
  for(i=0;i<kpar;i++) {
    for(j=0;j<kpar;j++) {
      covar[i][j]/=(-2*log(conflim));
    }
  } 


  



  free_matrix_d(bb,kpar-1,1);
  free(y);
  free(parconf);
  free(sigparconf);

  free_matrix_d(parelip,kpar  ,nconflini);
  free_matrix_d(pareA  ,kpar-1,nconflini);
  free_matrix_d(pareB  ,kpar-1,nconflini);
  free_matrix_d(pareC  ,kpar-1,nconflini);

  free_matrix_d(invcovar,kpar  ,kpar);
  free_matrix_d(invcovA ,kpar-1,kpar-1);
  free_matrix_d(invcovB ,kpar-1,kpar-1);
  free_matrix_d(invcovC ,kpar-1,kpar-1);

}





double Funk1_norm_h_g_f_d(double x) {
  double firstsum;
  int i;
  double tmp1;
  int nstep=20;

  xtmp=x;
  firstsum=0;
  for(i=0;i<kpar;i++) {
    tmp1=pp[i]*gaussintleg_d(Funk2_int_h_g_f_d,xklim[i],xklim[i+1],nstep);
    firstsum+=tmp1;
  }
  return(firstsum);
}


double Funk2_int_h_g_f_d(double x) {
  if(DEBUG3) {
    printf("Ret %f  x %f \n",Fermi(x,xf,Tf)*gaussian(xtmp,x,sigi),x);
    printf(" sigi %f xtmp %f\n",sigi,xtmp);
  }
  return(Fermi(x,xf,Tf)*gaussian(xtmp,x,sigi));
}
