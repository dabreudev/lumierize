#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_machine.h>
#include <cpgplot.h>
#include "alloc.h"
#include "mlhist.h"
#include "sthisto.h"
#include "amoeba.h"
#include "random.h"
#include "gaussj.h"
#include "minmax.h"
#include "elip.h"
#include "quartil.h"

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
#define TRYEMPIRICAL 0
#endif

/* #define USEPCUT 0 */
/* #define TOLERR 0.0001 */

/* Ahora mismo esta con el pcut definido, de modo que se aleja de valores 
   cercanos a 0. Puede afectar al calculo de las covarianzas!! */

double Amoe_Funk_ff_gg_d_main(int n, double *x, double *y, double *p);
double Amoe_Funk_ff_gg_d_conf(int n, double *x, double *y, double *p);
void   TeorCovars_ff_gg_d(int n,double *x,double *errx,double *y,double *erry,int kx,double *xk,double yff,double yinter,double *P,double **covar);


double *_errx_MLA_ff_gg_d;
double *_erry_MLA_ff_gg_d;
double *_xk_MLA_ff_gg_d;
double _yff_MLA_ff_gg_d;
double _yinter_MLA_ff_gg_d;
int _ndata_MLA_ff_gg_d;
int _iter_m_MLA_ff_gg_d;
int _iter_c_MLA_ff_gg_d;
int _kx_MLA_ff_gg_d;
int _ky_MLA_ff_gg_d;
int _kpar_MLA_ff_gg_d;
int _nconfl_MLA_ff_gg_d;
double _conflim_MLA_ff_gg_d;
/* double ML[5*MAXTRIES*MAXITER]; */
double _MLmax_MLA_ff_gg_d;
/* double *_pcut; */

int _nemp_MLA_ff_gg_d=2;

int  MLA_ff_gg_d(int n,double *x,double *errx, double *y, double *erry, int kx, double *xk, double yff, double yinter, double *Pk, double **covPk) 
{

  double *par;
  double *sigpar;
  double norm;
  double sparA;
  int **histo;
  int i,j,jx;
  double **covar;
  int iter_amo;
  int tries=0;

  /* No uso estas variables, asique las comento 
  double p0min,p0max;
  */
  
  /* Copias globales de estas variables para que se vean en toda la subrutina */
  _errx_MLA_ff_gg_d=errx;
  _erry_MLA_ff_gg_d=erry;
  _ndata_MLA_ff_gg_d=n;
  _kx_MLA_ff_gg_d=kx;
  _yff_MLA_ff_gg_d=yff;
  _yinter_MLA_ff_gg_d=yinter;
  _xk_MLA_ff_gg_d=xk;
  _kpar_MLA_ff_gg_d = _kx_MLA_ff_gg_d * 2;

  par=vector_d(_kpar_MLA_ff_gg_d);
  sigpar=vector_d(_kpar_MLA_ff_gg_d);
  /* _pcut=vector_d(_kpar); */
  covar=matrix_d(_kpar_MLA_ff_gg_d,_kpar_MLA_ff_gg_d);

  _iter_m_MLA_ff_gg_d=0;

  if(DEBUG3) {
    for(i=0;i<n;i++) printf(" Entrada x %g errx %g y %g erry %g\n",x[i],errx[i],y[i],erry[i]);
  }
  
  
  if((histo=StHisto2DFF_d(n,x,y,kx,xk,xk+kx,yff))==NULL) { 
    printf("I cannot dimension histo of %d elements \n",kx*2);
    exit(1);
  }
  
  
  /* Determinacion de las probabilidades iniciales en hk y ffk */
  for(jx=0;jx<kx;jx++)
  {
      if((histo[jx][0]+histo[jx][1])==0) 
           par[jx]= 1.  / (xk[jx+1]-xk[jx]);
      else                  
           par[jx]=((double)histo[jx][0]+(double)histo[jx][1])/(xk[jx+1]-xk[jx]);
  }
  for(jx=0;jx<kx;jx++)
  {
      if (histo[jx][0]==0) 
           par[kx + jx]= 1.;
      else
           if (histo[jx][1]==0)
	        par[kx + jx]= 0.001;                     
           else
	        par[kx + jx]= (double)histo[jx][1] / (xk[jx+1]-xk[jx]) / yinter / par[jx];
  }
  /* Valores muy cercanos a cero producen errores. Se adopta un valor inicial m�nimo del 0.1% */
  for(jx=0;jx<kx;jx++)
  {
      if (par[kx + jx] < 0.001)
           par[kx + jx] = 0.001;
  }
  
  /* Normalizo las pk */
  norm=0;
  for(jx=0;jx!=kx;++jx) {
       norm+=par[jx]*(xk[jx+1]-xk[jx]);
  }
  for(jx=0;jx!=kx;++jx) {
       par[jx]=par[jx]/norm;
  }

  /* Paso a bidimensional, teniendo en cuenta que las cosas que dependen de las probabilidadades son un vector, no una matriz como el histograma */   
  /* Ejemplo de las sigmas del m�todo bidimensional puro  
  for(jx=0;jx!=kx;++jx)
  {
    for(jy=0;jy!=ky;++jy)
    {
        if(histo[jx][jy]==0) sigpar[jx+kx*jy]=sqrt(par[jx+kx*jy]*(xk[jx+1]-xk[jx])*(yk[jy+1]-yk[jy])*(1-par[jx+kx*jy]*(xk[jx+1]-xk[jx])*(yk[jy+1]-yk[jy]))/1.);
        else            sigpar[jx+kx*jy]=sqrt(par[jx+kx*jy]*(xk[jx+1]-xk[jx])*(yk[jy+1]-yk[jy])*(1-par[jx+kx*jy]*(xk[jx+1]-xk[jx])*(yk[jy+1]-yk[jy]))/histo[jx][jy]);
    }
  }
  */
  
  /* Definici�n de las sigmas iniciales. Esto le da una idea de por donde empezar al m�todo */
  for(jx=0;jx!=kx;++jx)
  {
       if((histo[jx][0]+histo[jx][1])==0) 
           sigpar[jx]=sqrt(par[jx]*(xk[jx+1]-xk[jx])*(1-par[jx]*(xk[jx+1]-xk[jx])));
       else                  
           sigpar[jx]=sqrt(par[jx]*(xk[jx+1]-xk[jx])*(1-par[jx]*(xk[jx+1]-xk[jx]))/((double)histo[jx][0]+(double)histo[jx][1]));
  }
  
  for(jx=0;jx!=kx;++jx)
  {
       if((histo[jx][1])==0) 
       {
           sigpar[kx + jx]=sqrt(par[kx + jx]*par[jx]*(xk[jx+1]-xk[jx]*yinter)*(1 - par[kx + jx]*par[jx]*(xk[jx+1]-xk[jx])*yinter)); 
       }
       else 
       {
           sparA = sqrt(par[kx + jx]*par[jx]*(xk[jx+1]-xk[jx])*yinter*(1 - par[kx + jx]*par[jx]*(xk[jx+1]-xk[jx])*yinter)) / (double)histo[jx][1];
           sigpar[kx + jx]=par[kx + jx]*sqrt( (sparA * par[jx] / par[kx + jx])*(sparA * par[jx] / par[kx + jx]) + (sigpar[jx] / par[jx])*(sigpar[jx] / par[jx]));
       }
  }
    
  free_matrix_i(histo,kx,2);
  
  /* Cambio de variable tanto en par como en sigpar */
  for (j=0;j!=_kpar_MLA_ff_gg_d;++j)
  {
    sigpar[j] = sigpar[j]/par[j];
    par[j] = log(par[j]);
  }
 
  do{
    iter_amo=Amoeba_d(n,x,y,_kpar_MLA_ff_gg_d,par,sigpar,FTOL,MAXITER,Amoe_Funk_ff_gg_d_main); 
    tries++;
    if(DEBUG) printf(" Exited from try %d with %d iterations\n",tries, iter_amo);
    for(j=0;j<_kpar_MLA_ff_gg_d;j++) sigpar[j]=sqrt(fabs(par[j]));
  }   while(iter_amo>= MAXITER && tries <=MAXTRIES );


  /* Aunque la hemos impuesto dentro de Anoeb_funk, volvemos a normalizar. La normalizaci�n se impone sobre las p originales, no sobre p' */
  
  norm=0;    
  for(jx=0;jx!=kx;++jx)
  { 
       norm+=exp(par[jx])*(xk[jx+1]-xk[jx]);
  }
  for(jx=0;jx!=kx;++jx)
  {
       par[jx]=log(exp(par[jx])/norm);
  }

  _MLmax_MLA_ff_gg_d=Amoe_Funk_ff_gg_d_main(n,x,y,par);
  
  /* Guardo la soluci�n "par" en "Pk"*/
  for(j=0;j<_kpar_MLA_ff_gg_d;j++)
  {
    Pk[j] = par[j];	
  }
  

  /* Estimacion de los errores en los parametros. Aqui no hace falta cambiar nada, uso los Pk recuperados con la formula. Ver si se ha escapado algo. De momento lo comento*/

    
  TeorCovars_ff_gg_d(n,x,errx,y,erry,kx,xk,yff,yinter,Pk,covar);   
  if(DEBUG) printf(" Calculo teorico\n");
  /* Aqui pinto las covarianzas que me molan */
   
  for(i=0;i<_kpar_MLA_ff_gg_d;i++) {
    for(j=0;j<_kpar_MLA_ff_gg_d;j++) {
      /* covPk[i][j] = 1.; */
      covPk[i][j]=covar[i][j];
      if(DEBUG) printf(" covar %d %d = %g\n",i,j,covar[i][j]);  
    }
  }

  /* Paso del c�lculo emp�rico 
  if(TRYEMPIRICAL) 
  {
    _conflim=exp(-.5/10.);
    EmpiricalCovars_hh_gg_d(n,x,errx,y,erry,kx,xk,ky,yk,Pk,sigpar,covar);
    if(DEBUG) printf(" Calculo empirico\n");
    for(i=0;i<_kpar;i++) {D.yff=(double)readf(HD.yff);
 
      for(j=0;j<_kpar;j++) {
        covPk[i][j]=covar[i][j];
        if(DEBUG2) printf(" covar %d %d = %g\n",i,j,covar[i][j]);
      }
    }
  }

  free_matrix_d(covar,_kpar,_kpar); 
  free(par);
  free(sigpar);

  if(DEBUG) printf(" MLF %g\n",_MLmax); 
  */

  if(iter_amo>=MAXITER-1) return(2);
  return(0);
}


double Amoe_Funk_ff_gg_d_main(int n, double *x, double *y, double *p) 
{

  int i,jx,j;
  double logL=0.;
  double ltmp;
  double tmp1,priori;
  double norm=0;
  (void)n;


  logL=0.;

  /* Imponemos la normalizacion de la funcion de probabilidad */
  norm=0;
  for(jx=0;jx<_kx_MLA_ff_gg_d;jx++)
  {
      /* Esto fuera, siempre tengo exp(p) positivo 
      if(p[jx+_kx*jy]<0) p[jx+_kx*jy]=-p[jx+_kx*jy]; */
      norm += exp(p[jx])*(_xk_MLA_ff_gg_d[jx+1]-_xk_MLA_ff_gg_d[jx]);
  }
 
  for(jx=0;jx<_kx_MLA_ff_gg_d;jx++) p[jx]=log(exp(p[jx])/norm);

  for(i=0;i<_ndata_MLA_ff_gg_d;i++)
  {
 
    ltmp=0;
    priori=1;
    for(jx=0;jx<_kx_MLA_ff_gg_d;jx++) 
    {
        /* M�todo 1 para la funci�n ML. Redefino el histograma y tomo como +/- inf los l�mites en A */
	tmp1 =  -exp(p[jx])/4. * (erf((x[i]-_xk_MLA_ff_gg_d[jx+1]) / _errx_MLA_ff_gg_d[i]/sqrt(2.)) - erf((x[i]-_xk_MLA_ff_gg_d[jx]) / _errx_MLA_ff_gg_d[i]/sqrt(2.))) * (1 + (2*exp(p[_kx_MLA_ff_gg_d+jx])-1)*erf((y[i]-_yff_MLA_ff_gg_d)/_erry_MLA_ff_gg_d[i]/sqrt(2.))); 
	
	/* M�todo 2. Solo redefino el histograma */
	/* tmp1 = exp(p[jx])/4. * (erf((x[i]-_xk_MLA_ff_gg_d[jx+1]) / _errx_MLA_ff_gg_d[i] / sqrt(2.)) - erf((x[i]-_xk_MLA_ff_gg_d[jx]) / _errx_MLA_ff_gg_d[i] / sqrt(2.))) * (exp(p[_kx_MLA_ff_gg_d+jx]) * (erf((y[i]-_yff_MLA_ff_gg_d-_yinter_MLA_ff_gg_d) / _erry_MLA_ff_gg_d[i] / sqrt(2.)) - erf((y[i]-_yff_MLA_ff_gg_d) / _erry_MLA_ff_gg_d[i] / sqrt(2.))) + (1 - exp(p[_kx_MLA_ff_gg_d+jx])) * (erf((y[i]-_yff_MLA_ff_gg_d) / _erry_MLA_ff_gg_d[i] / sqrt(2.)) - erf((y[i]-_yff_MLA_ff_gg_d+_yinter_MLA_ff_gg_d) / _erry_MLA_ff_gg_d[i] / sqrt(2.)))); */
        ltmp+=tmp1;  /* Pasando de cambiar los DEBUG, al menos de momento. */
    }
    /* if(!USEPCUT) priori=1;
    if(ltmp!=0)  logL-=log(ltmp*priori); Me quito el Pcut */
    if(ltmp!=0)  logL-=log(ltmp);
    else         logL-=log(FTOL);
  }
  if(DEBUG) printf(" MLF iter %3d logL %.15g ",_iter_m_MLA_ff_gg_d,logL);
  if(DEBUG) {
    printf(" Llamado con ");
    for(j=0;j<_kpar_MLA_ff_gg_d;j++) printf(" %5.10f",p[j]);
    printf("\n");
  }

  _iter_m_MLA_ff_gg_d++;
  return(logL);
}


 double Amoe_Funk_ff_gg_d_conf(int n, double *x, double *y, double *p)
{
   (void)n;
   (void)x;
   (void)y;
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
 

void TeorCovars_ff_gg_d(int n,double *x,double *errx,double *y,double *erry,int kx,double *xk,double yff, double yinter, double *Pk,double **covar) 
{
  
  int i,l,lx,m,mx,jx;
  double tmp1,tmp2,tmp3,tmp4;
  double **hessian;
  double **b;
  (void)yinter;
   
  hessian=matrix_d(_kpar_MLA_ff_gg_d+1,_kpar_MLA_ff_gg_d+1);
  b=matrix_d(_kpar_MLA_ff_gg_d+1,1);

  /* Definimos el hessiano de la funci�n minimizada y tenemos en cuenta la normalizaci�n */
  /* Definici�n de los elementos del Hessiano correspondientes a las pk - pk */ 
  for(lx=0;lx<kx;lx++) 
  {
     for(mx=0;mx<kx;mx++)
     {
        tmp1 = 0.;
        for(i=0;i<n;i++)
	{
           tmp2=0;
	   for(jx=0;jx<kx;jx++)
	   {
	         tmp2 += -exp(Pk[jx])/4. * (erf((x[i] - xk[jx+1])/ errx[i]/sqrt(2.)) - erf((x[i] - xk[jx]) / errx[i]/sqrt(2.))) * (1 + (2*exp(Pk[kx+jx])-1) * erf((y[i] - yff) / erry[i] / sqrt(2.)));
           }
	   tmp3 = exp(Pk[lx])/4. * (erf((x[i] - xk[lx+1])/ errx[i]/sqrt(2.)) - erf((x[i] - xk[lx]) / errx[i]/sqrt(2.))) * (1 + (2*exp(Pk[kx+lx])-1) * erf((y[i]- yff) / erry[i] / sqrt(2.)));
	   tmp4 = exp(Pk[mx])/4. * (erf((x[i] - xk[mx+1])/ errx[i]/sqrt(2.)) - erf((x[i] - xk[mx]) / errx[i]/sqrt(2.))) * (1 + (2*exp(Pk[kx+mx])-1) * erf((y[i]- yff) / erry[i] / sqrt(2.)));
	   if (lx == mx)
	        /* � Tengo un t�rmino libre que depende de lamnda en est� caso ! Lo escojo de forma que anule el primer t�rmino (prueba)
	        tmp1 += (-tmp3*tmp2 - tmp3*tmp4)/(tmp2*tmp2);*/
		tmp1 += (-tmp3*tmp4)/(tmp2*tmp2);
	   else
	        tmp1 += (-tmp3*tmp4)/(tmp2*tmp2);
       }
     hessian[lx][mx]=-tmp1; 
     }
  }
   
   /* Definici�n de los elementos del Hessiano correspondientes a las pk - ffk*/
  for(lx=0;lx<kx;lx++) 
  {
     for(mx=0;mx<kx;mx++)
     {
        tmp1 = 0.;
        for(i=0;i<n;i++)
	{
           tmp2=0;
	   for(jx=0;jx<kx;jx++)
	   {
	         tmp2 += -exp(Pk[jx])/4. * (erf((x[i] - xk[jx+1])/ errx[i]/sqrt(2.)) - erf((x[i] - xk[jx]) / errx[i]/sqrt(2.))) * (1 + (2*exp(Pk[kx+jx])-1) * erf((y[i]- yff) / erry[i] / sqrt(2.)));
           }
	   tmp3 = exp(Pk[lx])/4. * (erf((x[i] - xk[lx+1])/ errx[i]/sqrt(2.)) - erf((x[i] - xk[lx]) / errx[i]/sqrt(2.))) * (1 + (2*exp(Pk[kx+lx])-1) * erf((y[i]- yff) / erry[i] / sqrt(2.)));
	   tmp4 = exp(Pk[mx])*exp(Pk[kx + mx])/2. * (erf((x[i] - xk[mx+1])/ errx[i]/sqrt(2.)) - erf((x[i] - xk[mx]) / errx[i]/sqrt(2.))) * erf((y[i]- yff) / erry[i] / sqrt(2.));
           if (lx == mx)
	        tmp1 += (-tmp4*tmp2 - tmp3*tmp4)/(tmp2*tmp2);
	   else
	        tmp1 += (-tmp3*tmp4)/(tmp2*tmp2);
        }
      hessian[lx][kx + mx]=-tmp1;
      hessian[kx + lx][mx]=-tmp1;
      }
   }
   
  /* Definici�n de los elementos del Hessiano correspondientes a las ffk - ffk*/
  for(lx=0;lx<kx;lx++) 
  {
     for(mx=0;mx<kx;mx++)
     {
        tmp1 = 0.;
        for(i=0;i<n;i++)
	{
           tmp2=0;
	   for(jx=0;jx<kx;jx++)
	   {
	         tmp2 += -exp(Pk[jx])/4. * (erf((x[i] - xk[jx+1])/ errx[i]/sqrt(2.)) - erf((x[i] - xk[jx]) / errx[i]/sqrt(2.))) * (1 + (2*exp(Pk[kx+jx])-1) * erf((y[i]- yff) / erry[i] / sqrt(2.)));
           }
	   tmp3 = exp(Pk[lx])*exp(Pk[kx + mx])/2. * (erf((x[i] - xk[lx+1])/ errx[i]/sqrt(2.)) - erf((x[i] - xk[lx]) / errx[i]/sqrt(2.))) * erf((y[i]- yff) / erry[i] / sqrt(2.));
	   tmp4 = exp(Pk[mx])*exp(Pk[kx + mx])/2. * (erf((x[i] - xk[mx+1])/ errx[i]/sqrt(2.)) - erf((x[i] - xk[mx]) / errx[i]/sqrt(2.))) * erf((y[i]- yff) / erry[i] / sqrt(2.));
	   if (lx == mx)
	        tmp1 += (-tmp3*tmp2 - tmp3*tmp3)/(tmp2*tmp2);
	   else
	        tmp1 += (-tmp3*tmp4)/(tmp2*tmp2);

        }
      hessian[kx + lx][kx + mx]=-tmp1;
      }
   }
  
  /* Definici�n de los elementos del Hessiano correspondientes a las pk - lamda; ffk - lambda y lamda - lambda*/ 
   for(lx=0;lx<kx;lx++) 
  {
      hessian[_kpar_MLA_ff_gg_d][lx]=(xk[lx+1] - xk[lx])*exp(Pk[lx]);
      hessian[lx][_kpar_MLA_ff_gg_d]=(xk[lx+1] - xk[lx])*exp(Pk[lx]);
      hessian[_kpar_MLA_ff_gg_d][kx+lx]=0;
      hessian[kx+lx][_kpar_MLA_ff_gg_d]=0;
  }
  hessian[_kpar_MLA_ff_gg_d][_kpar_MLA_ff_gg_d]=0;
    
  
  /* Invierto el menos hessiano para meterlo en la matriz de covarianzas */
  gaussj_d(hessian,_kpar_MLA_ff_gg_d,b,1);  /* Me supongo que el hessiano invertido queda almacenado en hessiano */

  /* Pillo la matriz de covarianzas del hessiano invertido y le cambio el signo*/
  for(l=0;l<_kpar_MLA_ff_gg_d;l++) 
  {
    for(m=0;m<_kpar_MLA_ff_gg_d;m++) 
    {
      covar[l][m]=hessian[l][m];     
    }
  }
  free_matrix_d(hessian,_kpar_MLA_ff_gg_d+1,_kpar_MLA_ff_gg_d+1);
  free_matrix_d(b,_kpar_MLA_ff_gg_d+1,1);
}

