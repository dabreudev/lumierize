#include "modulos.h"
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
void   EmpiricalCovars_ff_gg_d(int n,double *x,double *errx,int k,double *xk,double *Pk,double *sigPk, double **covar);
/*void   EmpiricalCovars_hh_gg_d(int n,double *x,double *errx,double *y,double *erry,int kx,double *xk,int ky,double *yk,double *Pk,double *sigPk,double **covar);*/
int  **StHisto2DFF_d(int n, double *x, double *y, int nbinx, double *xmin, double *xmax, double yff);



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
  /* Valores muy cercanos a cero producen errores. Se adopta un valor inicial mínimo del 0.1% */
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
  /* Ejemplo de las sigmas del método bidimensional puro  
  for(jx=0;jx!=kx;++jx)
  {
    for(jy=0;jy!=ky;++jy)
    {
        if(histo[jx][jy]==0) sigpar[jx+kx*jy]=sqrt(par[jx+kx*jy]*(xk[jx+1]-xk[jx])*(yk[jy+1]-yk[jy])*(1-par[jx+kx*jy]*(xk[jx+1]-xk[jx])*(yk[jy+1]-yk[jy]))/1.);
        else            sigpar[jx+kx*jy]=sqrt(par[jx+kx*jy]*(xk[jx+1]-xk[jx])*(yk[jy+1]-yk[jy])*(1-par[jx+kx*jy]*(xk[jx+1]-xk[jx])*(yk[jy+1]-yk[jy]))/histo[jx][jy]);
    }
  }
  */
  
  /* Definición de las sigmas iniciales. Esto le da una idea de por donde empezar al método */
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


  /* Aunque la hemos impuesto dentro de Anoeb_funk, volvemos a normalizar. La normalización se impone sobre las p originales, no sobre p' */
  
  norm=0;    
  for(jx=0;jx!=kx;++jx)
  { 
       norm+=exp(par[jx])*(xk[jx+1]-xk[jx]);
  }
  for(jx=0;jx!=kx;++jx)
  {
       par[jx]=log(exp(par[jx])/norm);
  }

  _MLmax_MLA_ff_gg_d=Amoe_Funk_hh_gg_d_main(n,x,y,par);
  
  /* Guardo la solución "par" en "Pk"*/
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

  /* Paso del cálculo empírico 
  if(TRYEMPIRICAL) 
  {
    _conflim=exp(-.5/10.);
    /*EmpiricalCovars_hh_gg_d(n,x,errx,y,erry,kx,xk,ky,yk,Pk,sigpar,covar);
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
  /* free(_pcut); 

  if(DEBUG) printf(" MLF %g\n",_MLmax); 
  */

  if(iter_amo>=MAXITER-1) return(2);
  return(0);
}


double Amoe_Funk_ff_gg_d_main(int n, double *x, double *y, double *p) 
{

  int i,jx,jy,j;
  double logL=0.;
  double ltmp;
  double tmp1,priori;
  double norm=0;

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
        /* Método 1 para la función ML. Redefino el histograma y tomo como +/- inf los límites en A */
	tmp1 =  -exp(p[jx])/4. * (erf((x[i]-_xk_MLA_ff_gg_d[jx+1]) / _errx_MLA_ff_gg_d[i]/sqrt(2.)) - erf((x[i]-_xk_MLA_ff_gg_d[jx]) / _errx_MLA_ff_gg_d[i]/sqrt(2.))) * (1 + (2*exp(p[_kx_MLA_ff_gg_d+jx])-1)*erf((y[i]-_yff_MLA_ff_gg_d)/_erry_MLA_ff_gg_d[i]/sqrt(2.))); 
	
	/* Método 2. Solo redefino el histograma */
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
   
  hessian=matrix_d(_kpar_MLA_ff_gg_d+1,_kpar_MLA_ff_gg_d+1);
  b=matrix_d(_kpar_MLA_ff_gg_d+1,1);

  /* Definimos el hessiano de la función minimizada y tenemos en cuenta la normalización */
  /* Definición de los elementos del Hessiano correspondientes a las pk - pk */ 
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
	        /* ¡ Tengo un término libre que depende de lamnda en esté caso ! Lo escojo de forma que anule el primer término (prueba)
	        tmp1 += (-tmp3*tmp2 - tmp3*tmp4)/(tmp2*tmp2);*/
		tmp1 += (-tmp3*tmp4)/(tmp2*tmp2);
	   else
	        tmp1 += (-tmp3*tmp4)/(tmp2*tmp2);
       }
     hessian[lx][mx]=-tmp1; 
     }
  }
   
   /* Definición de los elementos del Hessiano correspondientes a las pk - ffk*/
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
   
  /* Definición de los elementos del Hessiano correspondientes a las ffk - ffk*/
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
  
  /* Definición de los elementos del Hessiano correspondientes a las pk - lamda; ffk - lambda y lamda - lambda*/ 
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



void   EmpiricalCovars_ff_gg_d(int n,double *x,double *errx,int k,double *xk,double *Pk,double *sigPk, double **covar) 
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
  
  nconfl=NCONFL*_kpar_MLA_ff_gg_d;
  nconflini=nconfl;


  bb=matrix_d(_kpar_MLA_ff_gg_d-1,1);
  y=vector_d(n);
  parconf=vector_d(k);
  sigparconf=vector_d(k);

  parelip=matrix_d(_kpar_MLA_ff_gg_d  ,nconfl);
  pareA  =matrix_d(_kpar_MLA_ff_gg_d-1,nconfl);
  pareB  =matrix_d(_kpar_MLA_ff_gg_d-1,nconfl);
  pareC  =matrix_d(_kpar_MLA_ff_gg_d-1,nconfl);

  invcovar=matrix_d(_kpar_MLA_ff_gg_d  ,_kpar_MLA_ff_gg_d  );
  invcovA =matrix_d(_kpar_MLA_ff_gg_d-1,_kpar_MLA_ff_gg_d-1);
  invcovB =matrix_d(_kpar_MLA_ff_gg_d-1,_kpar_MLA_ff_gg_d-1);
  invcovC =matrix_d(_kpar_MLA_ff_gg_d-1,_kpar_MLA_ff_gg_d-1);
  covA    =matrix_d(_kpar_MLA_ff_gg_d-1,_kpar_MLA_ff_gg_d-1);
  covB    =matrix_d(_kpar_MLA_ff_gg_d-1,_kpar_MLA_ff_gg_d-1);
  covC    =matrix_d(_kpar_MLA_ff_gg_d-1,_kpar_MLA_ff_gg_d-1);
  distmax=vector_d(2);



/*   conflim=exp(-.5/100.);     */         //Los puntos corresponderan a 1 sigma entre 4 de desviacion para dist. normal en par
/*   conflim=exp(-.5);   */           //Los puntos corresponderan a 1 sigma entre 4 de desviacion para dist. normal en par
/*   conflim=exp(-.5/100.);    */


  if(DEBUG3) printf(" log(conflim) %f \n",log(_conflim_MLA_ff_gg_d));
  
  for(i=0;i<nconfl;i++) {
    for(j=0;j<_kpar_MLA_ff_gg_d;j++) {
      parconf[j]=Pk[j]+sigPk[j]*3.*Gasdev(); 
      sigparconf[j]=sigPk[j]; 
    }
    if(DEBUG3) printf(" INIPARTEST %d %d\n",i,(int)(nconfl/2.));
    if(i>(int)(nconfl/2.)) {
      for(j=0;j<_kpar_MLA_ff_gg_d;j++) {
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

    _iter_c_MLA_ff_gg_d=0;
    _iter_c_MLA_ff_gg_d=Amoeba_d(n,x,y,k,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_ff_gg_d_conf);

    for(j=0;j<_kpar_MLA_ff_gg_d;j++) {
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
	for(l=0;l<_kpar_MLA_ff_gg_d;l++)    printf(" %g ",(parelip[l])[j]);
      }
      printf("\n");
    }
      
    if(_iter_c_MLA_ff_gg_d==0 && Amoe_Funk_hh_gg_d_conf(n,x,y,parconf)>FTOL2 ) i--;
  }
  
  
  /* Supongo que el centro de la elipse es el valor que maximiza ML */

  for(i=0;i<nconfl;i++)    for(j=0;j<_kpar_MLA_ff_gg_d;j++)       (parelip[j])[i]-=Pk[j];


  if(DEBUGPLOT) {
    for(i=0;i<nconfl;i++) {
      if(DEBUG3) printf(" PARCONFREAL %g %g\n",(parelip[0])[i],(parelip[1])[i]);
      cpgsci(_nemp_MLA_ff_gg_d);
      cpgpt1((parelip[0])[i]+Pk[0],(parelip[1])[i]+Pk[1],1);
    }
  }

/*   kpar=readi(kpar); */

  /* Detecto puntos que esten muy alejados de la supuesta elipse */
  for(j=0;j<_kpar_MLA_ff_gg_d;j++) {
    Quartil_d(nconfl,parelip[j],&first,&median,&third);
    distmax[j]=maxf(fabs(first),fabs(third));
  }  
  for(i=0;i<nconfl;i++) {
    for(j=0;j<_kpar_MLA_ff_gg_d;j++) {
      if(fabs((parelip[j])[i])>2*2*distmax[j]/1.35) {
	for(j=0;j<_kpar_MLA_ff_gg_d;j++) {
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
    for(j=0;j<_kpar_MLA_ff_gg_d-1;j++) {
      pareA[j][i]=parelip[j][i];
    }
  }
  MCElipN_d(nconfl,_kpar_MLA_ff_gg_d-1,pareA,invcovA);
  for(i=0;i<nconfl;i++) {
    for(j=0;j<_kpar_MLA_ff_gg_d-1;j++) {
      pareB[j][i]=parelip[j+1][i];
    }
  }
  MCElipN_d(nconfl,_kpar_MLA_ff_gg_d-1,pareB,invcovB);
  for(i=0;i<nconfl;i++) {
    for(j=0;j<_kpar_MLA_ff_gg_d-2;j++) {
      pareC[j][i]=parelip[j][i];
    }
    pareC[_kpar_MLA_ff_gg_d-2][i]=parelip[_kpar_MLA_ff_gg_d-1][i];
  }
  MCElipN_d(nconfl,_kpar_MLA_ff_gg_d-1,pareC,invcovC);

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
  for(i=0;i<_kpar_MLA_ff_gg_d-1;i++) {
    for(j=0;j<_kpar_MLA_ff_gg_d-1;j++) {
      covA[i][j]=invcovA[i][j];
      covB[i][j]=invcovB[i][j];
      covC[i][j]=invcovC[i][j];
    }
  } 
  gaussj_d(covA,_kpar_MLA_ff_gg_d-1,bb,1);
  gaussj_d(covB,_kpar_MLA_ff_gg_d-1,bb,1);
  gaussj_d(covC,_kpar_MLA_ff_gg_d-1,bb,1);
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
  for(i=0;i<_kpar_MLA_ff_gg_d-1;i++) {
    for(j=0;j<_kpar_MLA_ff_gg_d-1;j++) {
      covar[i][j]=covA[i][j];
    }
  } 

  if(DEBUG3) printf(" primer paso \n");
  /* Uso covB para la fila de abajo, la columna de la derecha y el extremo inferior derecha */

  for(j=1;j<_kpar_MLA_ff_gg_d;j++)       covar[_kpar_MLA_ff_gg_d-1][j]=covB[_kpar_MLA_ff_gg_d-2][j-1];
  for(i=1;i<_kpar_MLA_ff_gg_d;i++)       covar[i][_kpar_MLA_ff_gg_d-1]=covB[i-1][_kpar_MLA_ff_gg_d-2];
  if(DEBUG3) printf(" segund paso \n");
  /* Uso covC para el extremo superior derecha y el inferior izquierda  */
  covar[0][_kpar_MLA_ff_gg_d-1]=covC[0][_kpar_MLA_ff_gg_d-2];
  covar[_kpar_MLA_ff_gg_d-1][0]=covC[_kpar_MLA_ff_gg_d-2][0];
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
  for(i=0;i<_kpar_MLA_ff_gg_d;i++) {
    for(j=0;j<_kpar_MLA_ff_gg_d;j++) {
      covar[i][j]/=(-2*log(_conflim_MLA_ff_gg_d));
    }
  } 


 

  free_matrix_d(bb,_kpar_MLA_ff_gg_d-1,1);
  free(y);
  free(parconf);
  free(sigparconf);

  free_matrix_d(parelip,_kpar_MLA_ff_gg_d  ,nconflini);
  free_matrix_d(pareA  ,_kpar_MLA_ff_gg_d-1,nconflini);
  free_matrix_d(pareB  ,_kpar_MLA_ff_gg_d-1,nconflini);
  free_matrix_d(pareC  ,_kpar_MLA_ff_gg_d-1,nconflini);

  free_matrix_d(invcovar,_kpar_MLA_ff_gg_d  ,_kpar_MLA_ff_gg_d);
  free_matrix_d(invcovA ,_kpar_MLA_ff_gg_d-1,_kpar_MLA_ff_gg_d-1);
  free_matrix_d(invcovB ,_kpar_MLA_ff_gg_d-1,_kpar_MLA_ff_gg_d-1);
  free_matrix_d(invcovC ,_kpar_MLA_ff_gg_d-1,_kpar_MLA_ff_gg_d-1);

}
