#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_machine.h>
#include "alloc.h"
#include "mlhist.h"
#include "sthisto.h"

/*#define FTOL  1e-12  <-- Valor anterior*/
/*#define FTOL  1e-15  <-- Valos 2 */
/*#define FTOL  1e-17  <-- Valos 3 */
#define FTOL  1e-17
#define FTOL2 1e-6
#define FTOL3 1e-7
/* #define MAXITER  2000 >-- Valor anterior */
/* #define MAXITER  5000 >-- Valor2 */
/* #define MAXITER  7000 >-- Valor3 */
#define MAXITER  7000
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
#define TRYEMPIRICAL 0
#endif

/* #define USEPCUT 0 */
/* #define TOLERR 0.0001 */

/* Ahora mismo esta con el pcut definido, de modo que se aleja de valores 
   cercanos a 0. Puede afectar al calculo de las covarianzas!! */

double Amoe_Funk_hh_gg_d_main(int n, double *x, double *y, double *p);
double Amoe_Funk_hh_gg_d_conf(int n, double *x, double *y, double *p);
void   TeorCovars_hh_gg_d(int n,double *x,double *errx,double *y,double *erry,int kx,double *xk,int ky,double *yk,double *P,double **covar);
void   EmpiricalCovars_hh_gg_d(int n,double *x,double *errx,int k,double *xk,double *Pk,double *sigPk, double **covar);
/*void   EmpiricalCovars_hh_gg_d(int n,double *x,double *errx,double *y,double *erry,int kx,double *xk,int ky,double *yk,double *Pk,double *sigPk,double **covar);*/



double *_errx;
double *_erry;
double *_xk;
double *_yk;
int _ndata;
int _iter_m;
int _iter_c;
int _kx;
int _ky;
int _kpar;
int _nconfl;
double _conflim;
/* double ML[5*MAXTRIES*MAXITER]; */
double _MLmax;
/* double *_pcut; */

int _nemp=2;

int  MLA_hh_gg_d(int n,double *x,double *errx, double *y, double *erry, int kx, int ky, double *xk, double *yk, double *Pk, double **covPk) 
{

  double *par;
  double *sigpar;
  double norm;
  int **histo;
  int i,j,jx,jy;
  double **covar;
  int iter_amo;
  int tries=0;

  /* No una estas variables, asique las comento 
  double p0min,p0max;
  */

  /* Copias globales de estas variables para que se vean en toda la subrutina */
  _errx=errx;
  _erry=erry;
  _ndata=n;
  _kx=kx;
  _ky=ky;
  _xk=xk;
  _yk=yk;
  _kpar = _kx * _ky;

  par=vector_d(_kpar);
  sigpar=vector_d(_kpar);
  /* _pcut=vector_d(_kpar); */
  covar=matrix_d(_kpar,_kpar);

  _iter_m=0;

  if(DEBUG3) {
    for(i=0;i<n;i++) printf(" Entrada x %g errx %g y %g erry %g\n",x[i],errx[i],y[i],erry[i]);
  }
  
  if((histo=StHisto2D_d(n,x,y,kx,xk,xk+kx,ky,yk,yk+ky))==NULL) { 
    printf("I cannot dimension histo of %d elements \n",kx*ky);
    exit(1);
  }
  for(jx=0;jx<kx;jx++) 
  {
    for(jy=0;jy<ky;jy++) 
    {
      if(histo[jx][jy]==0)
        par[jx+kx*jy]= 1.  /(xk[jx+1]-xk[jx])/(yk[jy+1]-yk[jy]);
      else                  
        par[jx+kx*jy]=(double)histo[jx][jy]/(xk[jx+1]-xk[jx])/(yk[jy+1]-yk[jy]);
    }
  }
  norm=0;
  for(jx=0;jx!=kx;++jx) {
    for(jy=0;jy!=ky;++jy) {
      norm+=par[jx+kx*jy]*(xk[jx+1]-xk[jx])*(yk[jy+1]-yk[jy]);
    }
  }
  for(j=0;j<_kpar;j++)
  {
    par[j]=par[j]/norm;
  }

  /**** HASTA AQUI Kike, a partir de aqui empiezo yo */
  /* Paso a bidimensional, teniendo en cuenta que las cosas que dependen de las probabilidadades son un vector, no una matriz como el histograma */
    
  for(jx=0;jx!=kx;++jx)
  {
    for(jy=0;jy!=ky;++jy)
    {
        if(histo[jx][jy]==0) sigpar[jx+kx*jy]=sqrt(par[jx+kx*jy]*(xk[jx+1]-xk[jx])*(yk[jy+1]-yk[jy])*(1-par[jx+kx*jy]*(xk[jx+1]-xk[jx])*(yk[jy+1]-yk[jy]))/1.);
        else            sigpar[jx+kx*jy]=sqrt(par[jx+kx*jy]*(xk[jx+1]-xk[jx])*(yk[jy+1]-yk[jy])*(1-par[jx+kx*jy]*(xk[jx+1]-xk[jx])*(yk[jy+1]-yk[jy]))/histo[jx][jy]);
    }
  }
    
  free_matrix_i(histo,kx,ky);

  /* Cambio de variable tanto en par como en sigpar */
  for (j=0;j!=_kpar;++j)
  {
    sigpar[j] = sigpar[j]/par[j];
    par[j] = log(par[j]);
  }
 
  do{
    iter_amo=Amoeba_d(n,x,y,_kpar,par,sigpar,FTOL,MAXITER,Amoe_Funk_hh_gg_d_main);  /* Hay que introducir par2 = log(par) y sigpar2 = sigpar/par. Luego hay que deshacer el cambio...*/
    tries++;
    if(DEBUG) printf(" Exited from try %d with %d iterations\n",tries, iter_amo);
    for(j=0;j<_kpar;j++) sigpar[j]=sqrt(fabs(par[j]));
  }   while(iter_amo>= MAXITER && tries <=MAXTRIES );


  /* Aunque la hemos impuesto dentro de Anoeb_funk, volvemos a normalizar. La normalizaci�n se impone sobre las p originales, no sobre p' */
  
  norm=0;    
  for(jx=0;jx!=kx;++jx)  
      for(jy=0;jy!=ky;++jy) 
	norm+=exp(par[jx+kx*jy])*(xk[jx+1]-xk[jx])*(yk[jy+1]-yk[jy]);
  for(j=0;j<_kpar;j++) par[j]=log(exp(par[j])/norm);

  _MLmax=Amoe_Funk_hh_gg_d_main(n,x,y,par);
  
  /* Nueva soluci�n con el cambio de variable */
  for(j=0;j<_kpar;j++)
  {
    Pk[j] = exp(par[j]);
  }
  

  /* Estimacion de los errores en los parametros. Aqui no hace falta cambiar nada, uso los Pk recuperados con la formula. Ver si se ha escapado algo*/

  TeorCovars_hh_gg_d(n,x,errx,y,erry,kx,xk,ky,yk,Pk,covar);   
  if(DEBUG) printf(" Calculo teorico\n");
  /* Aqui realizo el cambio entre las covarianzas de p' y p multiplicando por exp(p'[i])=Pk[i] y exp(p'[j])=Pk[j] */
  for(i=0;i<_kpar;i++) {
    for(j=0;j<_kpar;j++) {
      /* covPk[i][j]=covar[i][j]*Pk[i]*Pk[j]; */
      covPk[i][j]=covar[i][j];
      if(DEBUG) printf(" covar %d %d = %g\n",i,j,covar[i][j]);  
    }
  }
  
  for(j=0;j<_kpar;j++)
  {
    Pk[j] = par[j];
  }

  if(TRYEMPIRICAL) 
  {
    _conflim=exp(-.5/10.);
    /*EmpiricalCovars_hh_gg_d(n,x,errx,y,erry,kx,xk,ky,yk,Pk,sigpar,covar); */
    if(DEBUG) printf(" Calculo empirico\n");
    for(i=0;i<_kpar;i++) {
      for(j=0;j<_kpar;j++) {
        covPk[i][j]=covar[i][j];
        if(DEBUG2) printf(" covar %d %d = %g\n",i,j,covar[i][j]);
      }
    }
  }

  free_matrix_d(covar,_kpar,_kpar); 
  free(par);
  free(sigpar);
  /* free(_pcut); */

  if(DEBUG) printf(" MLF %g\n",_MLmax); 

  if(iter_amo>=MAXITER-1) return(2);
  return(0);
}


double Amoe_Funk_hh_gg_d_main(int n, double *x, double *y, double *p) 
{

  int i,jx,jy,j;
  double logL=0.;
  double ltmp;
  double tmp1,priori;
  double norm=0;

  logL=0.;

  /* Imponemos la normalizacion de la funcion de probabilidad */
  norm=0;
  for(jx=0;jx<_kx;jx++) 
  {
    for(jy=0;jy<_ky;jy++) 
    {
      /* Esto fuera, siempre tengo exp(p) positivo 
      if(p[jx+_kx*jy]<0) p[jx+_kx*jy]=-p[jx+_kx*jy]; */
      norm+=exp(p[jx+_kx*jy])*(_xk[jx+1]-_xk[jx])*(_yk[jy+1]-_yk[jy]);
    }
  }
  for(j=0;j<_kpar;j++) p[j]=log(exp(p[j])/norm);

  for(i=0;i<_ndata;i++)
  {
    ltmp=0;
    priori=1;
    for(jx=0;jx<_kx;jx++) 
    {
      for(jy=0;jy<_ky;jy++) 
      {
        /* priori*=exp(-_pcut[jx+_kx*jy]/p[jx+_kx*jy]); Distribucion a priori.ME lo quito con el cambio de variable*/
	/* Nuevo tmp1. El cambio es exp(p[jx+_kx*jy]) donde pone p[jx+_kx*jy] */
	tmp1=exp(p[jx+_kx*jy])/4.*(erf((x[i]-_xk[jx+1])/_errx[i]/sqrt(2.))-erf((x[i]-_xk[jx])/_errx[i]/sqrt(2.)))*(erf((y[i]-_yk[jy+1])/_erry[i]/sqrt(2.))-erf((y[i]-_yk[jy])/_erry[i]/sqrt(2.)));
        ltmp+=tmp1;  /* Pasando de cambiar los DEBUG, al menos de momento. */
      }
    }
    /* if(!USEPCUT) priori=1;
    if(ltmp!=0)  logL-=log(ltmp*priori); Me quito el Pcut*/
    if(ltmp!=0)  logL-=log(ltmp);
    else         logL-=log(FTOL);
  }
  if(DEBUG) printf(" MLF iter %3d logL %.15g ",_iter_m,logL);
  if(DEBUG) {
    printf(" Llamado con ");
    for(j=0;j<_kpar;j++) printf(" %5.10f",p[j]);
    printf("\n");
  }

  _iter_m++;
  return(logL);
}


double Amoe_Funk_hh_gg_d_conf(int n, double *x, double *y, double *p)
{
  (void)y;/* To avoid warning */
  (void)n;/* To avoid warning */
  (void)x;
  (void)p;
   

/*   int i,j;
  double logL=0.;
  double ltmp;
  double tmp1,priori;
  double norm=0;
 
  logL=0.;
*/
  /* Imponemos la normalizacion de la funcion de probabilidad */

/*  norm=0;
  for(jx=0;jx<_kx;jx++) 
  {
    for(jy=0;jy<_ky;jy++) 
    {
      if(p[jx+_kx*jy]<0) p[jx+_kx*jy]=-p[jx+_kx*jy]; 
        norm+=p[jx+_kx*jy]*(_xk[jx+1]-_xk[jx])*(_yk[jy+1]-_yk[jy]);
    }
  }
  for(j=0;j<_kpar;j++) p[j]=p[j]/norm;

  for(i=0;i<_ndata;i++) {
    ltmp=0;
    priori=1;
    for(jx=0;jx<_kx;jx++) 
    {
      for(jy=0;jy<_ky;jy++) 
      {
        priori*=exp(-pcut[jx+_kx*jy]/p[jx+_kx*jy]);*/  /* Distribucion a priori. Para evitar valores iguales a 0 */
       /* tmp1=p[jx+_kx*jy]/4.*(erf((x[i]-_xk[jx+1])/_errx[i]/sqrt(2.))-erf((x[i]-_xk[jx])/_errx[i]/sqrt(2.)))*(erf((y[i]-_yk[jy+1])/_erry[i]/sqrt(2.))-erf((y[i]-_yk[jy])/_erry[i]/sqrt(2.)));
        ltmp+=tmp1; 
      }
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

  return(fabs(logL-(MLmax-log(conflim)))); */
  return 0;
}
 

void TeorCovars_hh_gg_d(int n,double *x,double *errx,double *y,double *erry,int kx,double *xk,int ky,double *yk,double *Pk,double **covar) 
{
  int i,jx,jy,l,lx,ly,m,mx,my;
  double tmp1,tmp2,tmp3,tmp4;
  double **hessian;
  double **b;

   
  hessian=matrix_d(_kpar+1,_kpar+1);
  b=matrix_d(_kpar+1,1);

  /* Definimos el hessiano de la funci�n minimizada y tenemos en cuenta la normalizaci�n */
  tmp2=0;
  for(ly=0;ly<ky;ly++) 
  {
    for(lx=0;lx<kx;lx++)
    {
      for(my=0;my<ky;my++)
      {   
        for(mx=0;mx<kx;mx++) 
        {
          tmp1=0;
          for(i=0;i<n;i++) 
          {
	     tmp2=0;
             for(jy=0;jy<ky;jy++)
	     {
	        for(jx=0;jx<kx;jx++)
	        {
	         tmp2+=Pk[jy*kx+jx]/4.*(erf((x[i]-xk[jx+1])/errx[i]/sqrt(2.))-erf((x[i]-xk[jx])/errx[i]/sqrt(2.)))*(erf((y[i]-yk[jy+1])/erry[i]/sqrt(2.))-erf((y[i]-yk[jy])/erry[i]/sqrt(2.)));
	        }
	     }
	     /* Introduzco un Pk[ly*kx+lx] y un Pk[my*kx+mx]*/
             tmp3=Pk[ly*kx+lx]/4.*(erf((x[i]-xk[lx+1])/errx[i]/sqrt(2.))-erf((x[i]-xk[lx])/errx[i]/sqrt(2.)))*(erf((y[i]-yk[ly+1])/erry[i]/sqrt(2.))-erf((y[i]-yk[ly])/erry[i]/sqrt(2.)));
	     tmp4=Pk[my*kx+mx]/4.*(erf((x[i]-xk[mx+1])/errx[i]/sqrt(2.))-erf((x[i]-xk[mx])/errx[i]/sqrt(2.)))*(erf((y[i]-yk[my+1])/erry[i]/sqrt(2.))-erf((y[i]-yk[my])/erry[i]/sqrt(2.)));
	     tmp1+=-tmp3*tmp4/(tmp2*tmp2);
          }
          hessian[ly*kx+lx][my*kx+mx]=tmp1;
        }
      }
    }  
  }
  for(ly=0;ly<ky;ly++) 
  {
    for(lx=0;lx<kx;lx++)
    {
      hessian[_kpar][ly*kx+lx]=-(xk[lx+1]-xk[lx])*(yk[ly+1]-yk[ly])*Pk[ly*kx+lx];
      hessian[ly*kx+lx][_kpar]=-(xk[lx+1]-xk[lx])*(yk[ly+1]-yk[ly])*Pk[ly*kx+lx];
    }
  }
  hessian[_kpar][_kpar]=0;
  
  
  printf("\n");

  /* Invierto el menos hessiano para meterlo en la matriz de covarianzas */
  gaussj_d(hessian,_kpar+1,b,1);   /* Me supongo que el hessiano invertido queda almacenado en hessiano */

  /* Pillo la matriz de covarianzas del hessiano invertido y le cambio el signo*/
  for(l=0;l<_kpar;l++) 
  {
    for(m=0;m<_kpar;m++) 
    {
      covar[l][m]=-hessian[l][m];        
    }
  }

  free_matrix_d(hessian,_kpar+1,_kpar+1);
  free_matrix_d(b,_kpar+1,1);
}



void   EmpiricalCovars_hh_gg_d(int n,double *x,double *errx,int k,double *xk,double *Pk,double *sigPk, double **covar) 
{

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
  
  int nconfl,nconflini;
  double first, median, third, *distmax;

/*   double a,b,c,d,f,e;  */
/*   double xcelip,ycelip,elipa,elipb,elipt;  */

/*   float xtmp,ytmp; */
  
  nconfl=NCONFL*_kpar;
  nconflini=nconfl;


  bb=matrix_d(_kpar-1,1);
  y=vector_d(n);
  parconf=vector_d(k);
  sigparconf=vector_d(k);

  parelip=matrix_d(_kpar  ,nconfl);
  pareA  =matrix_d(_kpar-1,nconfl);
  pareB  =matrix_d(_kpar-1,nconfl);
  pareC  =matrix_d(_kpar-1,nconfl);

  invcovar=matrix_d(_kpar  ,_kpar  );
  invcovA =matrix_d(_kpar-1,_kpar-1);
  invcovB =matrix_d(_kpar-1,_kpar-1);
  invcovC =matrix_d(_kpar-1,_kpar-1);
  covA    =matrix_d(_kpar-1,_kpar-1);
  covB    =matrix_d(_kpar-1,_kpar-1);
  covC    =matrix_d(_kpar-1,_kpar-1);
  distmax=vector_d(2);



/*   conflim=exp(-.5/100.);     */         //Los puntos corresponderan a 1 sigma entre 4 de desviacion para dist. normal en par
/*   conflim=exp(-.5);   */           //Los puntos corresponderan a 1 sigma entre 4 de desviacion para dist. normal en par
/*   conflim=exp(-.5/100.);    */


  if(DEBUG3) printf(" log(conflim) %f \n",log(_conflim));
  
  for(i=0;i<nconfl;i++) {
    for(j=0;j<_kpar;j++) {
      parconf[j]=Pk[j]+sigPk[j]*3.*Gasdev(); 
      sigparconf[j]=sigPk[j]; 
    }
    if(DEBUG3) printf(" INIPARTEST %d %d\n",i,(int)(nconfl/2.));
    if(i>(int)(nconfl/2.)) {
      for(j=0;j<_kpar;j++) {
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

    _iter_c=0;
    _iter_c=Amoeba_d(n,x,y,k,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_hh_gg_d_conf);

    for(j=0;j<_kpar;j++) {
      (parelip[j])[i]=parconf[j];
/*       phi=M_PI*Gasdev(); */
/*       teta=M_PI*Gasdev(); */
/*       (parelip[i])[0]=Pk[0]+1*cos(phi)*sin(teta); */
      /*       (parelip[i])[1]=Pk[1]+1*cos(phi)*cos(teta); */
      /*       (parelip[i])[2]=Pk[2]+1*sin(phi); */
    }
    if(DEBUG3) {
      printf(" PARCONFEL ");
      for(j=0;j<i;j++)  {
	printf(" ## ");
	for(l=0;l<_kpar;l++)    printf(" %g ",(parelip[l])[j]);
      }
      printf("\n");
    }
      
    if(_iter_c==0 && Amoe_Funk_hh_gg_d_conf(n,x,y,parconf)>FTOL2 ) i--;
  }
  
  
  /* Supongo que el centro de la elipse es el valor que maximiza ML */

  for(i=0;i<nconfl;i++)    for(j=0;j<_kpar;j++)       (parelip[j])[i]-=Pk[j];


  if(DEBUGPLOT) {
    for(i=0;i<nconfl;i++) {
      if(DEBUG3) printf(" PARCONFREAL %g %g\n",(parelip[0])[i],(parelip[1])[i]);
      cpgsci(_nemp);
      cpgpt1((parelip[0])[i]+Pk[0],(parelip[1])[i]+Pk[1],1);
    }
  }

/*   kpar=readi(kpar); */

  /* Detecto puntos que esten muy alejados de la supuesta elipse */
  for(j=0;j<_kpar;j++) {
    Quartil_d(nconfl,parelip[j],&first,&median,&third);
    distmax[j]=maxf(fabs(first),fabs(third));
  }  
  for(i=0;i<nconfl;i++) {
    for(j=0;j<_kpar;j++) {
      if(fabs((parelip[j])[i])>2*2*distmax[j]/1.35) {
	for(j=0;j<_kpar;j++) {
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
    for(j=0;j<_kpar-1;j++) {
      pareA[j][i]=parelip[j][i];
    }
  }
  MCElipN_d(nconfl,_kpar-1,pareA,invcovA);
  for(i=0;i<nconfl;i++) {
    for(j=0;j<_kpar-1;j++) {
      pareB[j][i]=parelip[j+1][i];
    }
  }
  MCElipN_d(nconfl,_kpar-1,pareB,invcovB);
  for(i=0;i<nconfl;i++) {
    for(j=0;j<_kpar-2;j++) {
      pareC[j][i]=parelip[j][i];
    }
    pareC[_kpar-2][i]=parelip[_kpar-1][i];
  }
  MCElipN_d(nconfl,_kpar-1,pareC,invcovC);

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
  for(i=0;i<_kpar-1;i++) {
    for(j=0;j<_kpar-1;j++) {
      covA[i][j]=invcovA[i][j];
      covB[i][j]=invcovB[i][j];
      covC[i][j]=invcovC[i][j];
    }
  } 
  gaussj_d(covA,_kpar-1,bb,1);
  gaussj_d(covB,_kpar-1,bb,1);
  gaussj_d(covC,_kpar-1,bb,1);
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
  for(i=0;i<_kpar-1;i++) {
    for(j=0;j<_kpar-1;j++) {
      covar[i][j]=covA[i][j];
    }
  } 

  if(DEBUG3) printf(" primer paso \n");
  /* Uso covB para la fila de abajo, la columna de la derecha y el extremo inferior derecha */

  for(j=1;j<_kpar;j++)       covar[_kpar-1][j]=covB[_kpar-2][j-1];
  for(i=1;i<_kpar;i++)       covar[i][_kpar-1]=covB[i-1][_kpar-2];
  if(DEBUG3) printf(" segund paso \n");
  /* Uso covC para el extremo superior derecha y el inferior izquierda  */
  covar[0][_kpar-1]=covC[0][_kpar-2];
  covar[_kpar-1][0]=covC[_kpar-2][0];
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
  for(i=0;i<_kpar;i++) {
    for(j=0;j<_kpar;j++) {
      covar[i][j]/=(-2*log(_conflim));
    }
  } 


 

  free_matrix_d(bb,_kpar-1,1);
  free(y);
  free(parconf);
  free(sigparconf);

  free_matrix_d(parelip,_kpar  ,nconflini);
  free_matrix_d(pareA  ,_kpar-1,nconflini);
  free_matrix_d(pareB  ,_kpar-1,nconflini);
  free_matrix_d(pareC  ,_kpar-1,nconflini);

  free_matrix_d(invcovar,_kpar  ,_kpar);
  free_matrix_d(invcovA ,_kpar-1,_kpar-1);
  free_matrix_d(invcovB ,_kpar-1,_kpar-1);
  free_matrix_d(invcovC ,_kpar-1,_kpar-1);

}
