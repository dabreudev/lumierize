#include "modulos.h"
#define FTOL  1e-12
//#define FTOL2 1e-6
#define FTOL2 1e-4
#define FTOL3 1e-7
#define MAXITER  2000
//#define MAXITER2 300
#define MAXITER2 100
#define MAXTRIES   5
#define DEBUG  1
#define DEBUG2 0
#define DEBUG3 0
#define DEBUGPLOT 0
#define CONTOURPLOT 0
#define NCONFL 5
#define TOLERR 0.07 

#ifdef TRYEMPIRICAL
#else 
#define TRYEMPIRICAL 1
#endif

/* #define TOLERR 0.0001 */

double Amoe_Funk_SWML_M_main(int n, double *x, double *y, double *p);
double Amoe_Funk_SWML_M_norm(int n, double *x, double *y, double *p);
double Amoe_Funk_SWML_M_prior(int n, double *x, double *y, double *p);
double Amoe_Funk_SWML_M_conf(int n, double *x, double *y, double *p);
double Funk2_int_SWML_M(double x);
double Funk1_norm_SWML_M(double x);
void   EmpiricalCovars_SWML_M(int n,double *magn,double *z,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Steplf_M *lf);

void   ComputeNorma_SWML_M(int n, double mlim, double strrad, double zlow, double zup, struct Steplf_M *lf);
struct cosmo_param *co;
double magl;

int ndata;
int iter_m;
int iter_c;
int nconfl;
double conflim;
double *pp;
double MLmax;
double xf,Tf;
double sigi;
double xtmp;
double *magbin;
int nbin;


int  MLA_SWML_M(int n,double *magn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_M *lf) {

  double *y;
  double *par;
  double *sigpar;
  double **covar;
  int i,j;
  int iter_amo;
  double *magabs;
  int *histo;
  double norm;
  double maghistomin,maghistomax;

  if(DEBUG) printf(" Estoy auiq\n");

  /* Copias globales de estas variables para que se vean en toda la subrutina */
  ndata=n;
  co=&cosmo;
  magl=mlim;
  nbin=lf->nbin;
  magbin=lf->magni;

  magabs=vector_d(n);

  y=vector_d(n);
  par=vector_d(lf->nbin);
  sigpar=vector_d(lf->nbin);
  covar=matrix_d(lf->nbin,lf->nbin);
  iter_m=0;

    for(i=0;i<n;i++) {
      magabs[i]=Mag(z[i],magn[i],*co);
      if(DEBUG3) printf(" Entrada x %g\n",magn[i]);
    }


  iter_amo=MAXITER+1;
  if(DEBUG3) printf(" iter_amo %d\n",iter_amo);

  maghistomin=magbin[0];
  maghistomax=magbin[lf->nbin];

  if((histo=StHisto2_d(n,magabs,lf->nbin,&(maghistomin),&(maghistomax)))==NULL) { printf("I cannot dimension histo of %d elements \n",lf->nbin);exit(1);}
  
  for(j=0;j<nbin;j++) {
    if(histo[j]==0)  par[j]=   1.           /(magbin[j+1]-magbin[j]);
    else             par[j]=(double)histo[j]/(magbin[j+1]-magbin[j]);
  }
  norm=0;
  for(j=0;j<nbin;j++) norm+=par[j]*(magbin[j+1]-magbin[j]);
  for(j=0;j<nbin;j++) { 
    par[j]=par[j]/norm;
    if(histo[j]==0) sigpar[j]=3*sqrt(par[j]*(magbin[j+1]-magbin[j])*(1-magbin[j]*(magbin[j+1]-magbin[j]))/1.);
    else            sigpar[j]=3*sqrt(par[j]*(magbin[j+1]-magbin[j])*(1-par[j]*(magbin[j+1]-magbin[j]))/histo[j]);
    if(DEBUG) printf(" %d par %f sig %f hisot %d magni %f - %f\n",j,par[j],sigpar[j],histo[j],magbin[j],magbin[j+1]);
    if(DEBUG) cpgsch(3.);
    if(DEBUG) cpgpt1((float)((magbin[j+1]+magbin[j])/2.),(float)(log10(par[j])),6);
    if(DEBUG) cpgsch(1.);
  }
  free(histo);

  for(j=0;j<nbin;j++) sigpar[j]=0.2;   /* Utilizo logaritmos mejor */
  for(j=0;j<nbin;j++) {
    if(histo[j]==0) par[j]=-30;
    else par[j]=log(par[j]);   /* Utilizo logaritmos mejor */
  }
/*   par[0]=-4.8;sigpar[0]=2.2; */
/*   par[1]=-3.5;sigpar[1]=2.2; */
/*   par[2]=-0.6;sigpar[2]=2.2; */
/*   par[3]=0.3;sigpar[3]=0.2; */

  if(DEBUG) for(j=0;j<nbin;j++) printf(" BEFORE NORM LOG %d par %f sig %f\n",j,par[j],sigpar[j]);

  norm=0;
  for(j=0;j<nbin;j++) norm+=exp(par[j])*(magbin[j+1]-magbin[j]);
  for(j=0;j<nbin;j++)   par[j]=par[j]-log(norm);
  if(DEBUG) for(j=0;j<nbin;j++) printf(" AFTER  NORM LOG %d par %f sig %f\n",j,par[j],sigpar[j]);


  


  printf(" Computing LF...\n");

  iter_amo=Amoeba_d(n,magn,z,lf->nbin,par,sigpar,FTOL,MAXITER,Amoe_Funk_SWML_M_main);


  MLmax=Amoe_Funk_SWML_M_main(n,magn,z,par);


  /* Meto la solucion en la salida */


  for(j=0;j<nbin;j++) lf->lf[j]=par[j];  


  /* Estimacion de los errores en los parametros */

  if(DEBUG) {
    printf(" FINAL : ");
    for(j=0;j<nbin;j++) printf(" %g ",par[j]); 
    printf(" \n");
  }



  if(TRYEMPIRICAL) {
    if(DEBUG) printf(" llamantry \n");
    conflim=exp(-.5/10.);
    EmpiricalCovars_SWML_M(n,magn,z,par,sigpar,mlim,cosmo,lf);  
    if(DEBUG) printf(" sale \n");
  }
  if(CONTOURPLOT) {
/*     nemp_f++; */
/*     conflim=exp(-.5*16.);     */
/*     EmpiricalCovars_SWML_M(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*9.);     */
/*     EmpiricalCovars_SWML_M(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*4.);     */
/*     EmpiricalCovars_SWML_M(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/1.);     */
/*     EmpiricalCovars_SWML_M(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/4.);     */
/*     EmpiricalCovars_SWML_M(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     cpgsci(1); */
/*     cpglab("P\\d1\\u","P\\d3\\u","Contornos de límites de confianza"); */
  }
  if(DEBUG) printf(" Calculo empirico\n");

  if(DEBUG) {
    printf(" Solution before norm:\n");
    for(j=0;j<nbin;j++) printf(" Inter_mag %11g - %11g    LF %11g (log=%7g) Err_LF %11g (log=%7g)\n",lf->magni[j],lf->magni[j+1],exp(lf->lf[j]),log10(exp(lf->lf[j])),exp(lf->lf[j])*lf->errlf[j],lf->errlf[j]*log(10.));
    for(j=0;j<nbin;j++)  printf(" %d par %f sig %f hisot %d magni %f - %f\n",j,par[j],sigpar[j],histo[j],magbin[j],magbin[j+1]);
  }
  ComputeNorma_SWML_M(n, mlim,  strrad, zlow,  zup, lf);
  
  if(DEBUG) {
    printf(" Solution after norm:\n");
    for(j=0;j<nbin;j++) printf(" Inter_mag %11g - %11g    LF %11g (log=%7g) Err_LF %11g (log=%7g)\n",lf->magni[j],lf->magni[j+1],exp(lf->lf[j]),log10(exp(lf->lf[j])),exp(lf->lf[j])*lf->errlf[j],lf->errlf[j]*log(10.));
    for(j=0;j<nbin;j++)  printf(" %d par %f sig %f hisot %d magni %f - %f\n",j,par[j],sigpar[j],histo[j],magbin[j],magbin[j+1]);
  }
  
/*   conflim=readf(conflim); */

  free(y);
  free(par);
  free(sigpar);
  free_matrix_d(covar,lf->nbin,lf->nbin);

  if(iter_amo>=MAXITER-1) return(2);
  return(0);
}


double Amoe_Funk_SWML_M_main(int n, double *x, double *y, double *p) {

  int i,j; 
  int jend;
  double logL=0.;
  double Mabs;
  double intsch;
  double Mlow;
  double norm;
  double funl;

  int jbin;

/*   for(j=0;j<nbin;j++) printf(" EN %f \n ",p[j]); */

  norm=0;
  for(j=0;j<nbin;j++) norm+=exp(p[j])*(magbin[j+1]-magbin[j]);
  for(j=0;j<nbin;j++) p[j]=p[j]-log(norm);

  if(DEBUG) {
    cpgsci(2);
    cpgsch(2.);
    for(j=0;j<nbin;j++) cpgpt1((float)((magbin[j+1]+magbin[j])/2.),(float)(log10(exp(p[j]))),5);
    if(DEBUG3) {
      printf("FL:\n");
      for(j=0;j<nbin;j++) printf("  %f %f %f ",(magbin[j+1]+magbin[j])/2.,exp(p[j]),log10(exp(p[j])));
      printf("\n");
    }
  }
  logL=0.;
  for(i=0;i<ndata;i++) {
    Mabs=Mag(y[i],x[i],*co);
    Mlow=Mag(y[i],magl,*co);
    intsch=0;
    if(Mlow<magbin[0])  jend=-1;
    if(Mlow>magbin[nbin])   jend=nbin;
    for(j=0;j<nbin;j++) if(Mlow>magbin[j] && Mlow<magbin[j+1]) jend=j;
    if(jend!=-1   && jend!=nbin)   intsch+=exp(p[jend])*(Mlow-magbin[jend]);
    for(j=0;j<jend;j++) 	      intsch+=exp(p[j])*(magbin[j+1]-magbin[j]);
    if(Mabs<magbin[0] || Mabs>magbin[nbin]) {
      jbin=-1;
    }
    else {
      for(j=0;j<nbin;j++) {
	if(DEBUG3) printf(" j %d magbin[j] %f magbin[j+1] %f \n",j,magbin[j],magbin[j+1]);
	if(magbin[j]<Mabs && Mabs<magbin[j+1]) {
	  funl=exp(p[j]);
	  jbin=j;
	  break;
	}
      }
      logL-= log(funl/intsch);
    }
    if(DEBUG3) printf(" ndata %d  logL %g x %f  y %f Mabs %f  jbin %d funl %g int %g\n",i,logL,x[i],y[i],Mabs,jbin,funl,intsch);   
/*     printf(" LF: lfamo.Mstar %g lfamo.phistar %f lf.alfa  %f\n",lfamo.Mstar,lfamo.phistar,lfamo.alfa); */
  }


  if(DEBUG) {
    cpgsci(0);
    cpgsch(2.);
    for(j=0;j<nbin;j++) cpgpt1((float)((magbin[j+1]+magbin[j])/2.),(float)(log10(exp(p[j]))),5);
    if(DEBUG3) {
      printf(" iter %d logL %f       ",iter_m,logL);
      for(j=0;j<nbin;j++) printf(" %g ",p[j]);
      printf("\n");
    }
  }
  iter_m++;
  return(logL);
}


double Amoe_Funk_SWML_M_conf(int n, double *x, double *y, double *p) {

/*   printf(" MLM %g conf %f   %f %f\n",MLmax,conflim,fabs(Amoe_Funk_SWML_M_main(n,x,y,p)-(MLmax-log(conflim))),Amoe_Funk_SWML_M_main(n,x,y,p)); */
   
  return(fabs(Amoe_Funk_SWML_M_main(n,x,y,p)-(MLmax-log(conflim)))); 
}



void   EmpiricalCovars_SWML_M(int n,double *magn,double *z,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Steplf_M *lf) {

  int i,j;  
  double *parconf; 
  double *sigparconf; 
  double **invcovar;
  double **covar;
  double **parelip;
  double *y;

  double **bb;
  double **pareA;
  double **pareB;
  double **pareC;
  
  double **invcovA;
  double **invcovB;
  double **invcovC;
  double **covA;
  double **covB;
  double **covC;

  double first, median, third, *distmax;

  int nconfl,nconflini;

  if(DEBUG3) printf(" n vale %d \n",n);
  nconfl=NCONFL*lf->nbin;
  nconfl=2*lf->nbin*(lf->nbin+1)/2;
  nconflini=nconfl;

  bb=matrix_d(lf->nbin-1,1);
  for(i=0;i<lf->nbin-1;i++) (bb[i])[0]=0;
  y=vector_d(n);
  parconf=vector_d(lf->nbin);
  sigparconf=vector_d(lf->nbin);
  parelip=matrix_d(lf->nbin  ,nconflini);
  pareA  =matrix_d(lf->nbin-1,nconflini);
  pareB  =matrix_d(lf->nbin-1,nconflini);
  pareC  =matrix_d(lf->nbin-1,nconflini);

  invcovar=matrix_d(lf->nbin  ,lf->nbin  );
  covar   =matrix_d(lf->nbin  ,lf->nbin  );
  invcovA =matrix_d(lf->nbin-1,lf->nbin-1);
  invcovB =matrix_d(lf->nbin-1,lf->nbin-1);
  invcovC =matrix_d(lf->nbin-1,lf->nbin-1);
  covA    =matrix_d(lf->nbin-1,lf->nbin-1);
  covB    =matrix_d(lf->nbin-1,lf->nbin-1);
  covC    =matrix_d(lf->nbin-1,lf->nbin-1);
  distmax =vector_d(lf->nbin);
 

  printf(" Error step  ");
  for(i=0;i<nconfl;i++) printf(".");
  for(i=0;i<nconfl;i++) printf("\b");
  fflush(NULL);

  for(i=0;i<nconfl;i++) {
    printf("#");
/*     printf(" Por %d %d \n",i,nconfl); */
    fflush(NULL);
    for(j=0;j<lf->nbin;j++) {
      parconf[j]=par[j]+3*sigpar[0]*Gasdev();
      sigparconf[j]=sigpar[j];
    }
    if(i>(int)(nconfl/2.)) {
      for(j=0;j<lf->nbin;j++) {
	parconf[j]=par[j]-((parelip[j])[(int)(i-nconfl/2.)+1]-par[j]);
	sigparconf[j]=((parelip[j])[(int)(i-nconfl/2.)+1]-par[j])/2.; 
      }
    }
    iter_c=0;
    if(DEBUG3) printf(" antes a nmo\n"); 
    iter_c=Amoeba_d(n,magn,z,lf->nbin,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_SWML_M_conf);
    if(DEBUG) printf(" %d CONF iter %d par0 %10f par1 %10f par2 %10f\n", i,  iter_c ,parconf[0],parconf[1],parconf[2]);
    for(j=0;j<lf->nbin;j++)  (parelip[j])[i]=parconf[j];
    if(iter_c==0 && Amoe_Funk_SWML_M_conf(n,magn,z,parconf)>FTOL2 ) {
      printf("\b.\b");
      i--;
      fflush(NULL);
    }
  }

  if(DEBUG3) printf(" Sali del bucle nconfl\n");
  
  if(DEBUG3) {
    for(i=0;i<nconfl;i++) {
      printf(" %d ",i);
      for(j=0;j<lf->nbin;j++)    printf("  %f ",(parelip[j])[i]);
      printf(" \n");
    }
    
  }

  /* Supongo que el centro de la elipse es el valor que maximiza ML */
  for(i=0;i<nconfl;i++) {
    for(j=0;j<lf->nbin;j++)       (parelip[j])[i]-=par[j];
  }

  if(DEBUG3) printf(" Estamos a punto de Elip\n");
  
  /* Detecto puntos que esten muy alejados de la supuesta elipse */
  for(j=0;j<lf->nbin;j++) {
    Quartil_d(nconfl,parelip[j],&first,&median,&third);
    distmax[j]=maxf(fabs(first),fabs(third));
  }  
  for(i=0;i<nconfl;i++) {
    for(j=0;j<lf->nbin;j++) {
      if(fabs((parelip[j])[i])>2*2*distmax[j]/1.35) {
	for(j=0;j<lf->nbin;j++) {
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
    for(j=0;j<lf->nbin-1;j++) {
      if(DEBUG3) printf(" %f ",parelip[j][i]);
      pareA[j][i]=parelip[j][i];
    }
    if(DEBUG3) printf("\n");
  }
  MCElipN_d(nconfl,lf->nbin-1,pareA,invcovA);
  for(i=0;i<nconfl;i++) {
    for(j=0;j<lf->nbin-1;j++) {
      pareB[j][i]=parelip[j+1][i];
    }
  }
  MCElipN_d(nconfl,lf->nbin-1,pareB,invcovB);
  for(i=0;i<nconfl;i++) {
    for(j=0;j<lf->nbin-2;j++) {
      pareC[j][i]=parelip[j][i];
    }
    pareC[lf->nbin-2][i]=parelip[lf->nbin-1][i];
  }
  MCElipN_d(nconfl,lf->nbin-1,pareC,invcovC);

  if(DEBUG3) printf(" Ya he calculado las elip\n");
  


  for(i=0;i<lf->nbin-1;i++) {
    for(j=0;j<lf->nbin-1;j++) {
      if(DEBUG3) printf(" A %f B %f C %f ",invcovA[i][j],invcovB[i][j],invcovC[i][j]);
      covA[i][j]=invcovA[i][j];
      covB[i][j]=invcovB[i][j];
      covC[i][j]=invcovC[i][j];
    }
    if(DEBUG3) printf("\n");
  } 
  if(DEBUG3) printf(" COVA antes gauuss\n");
  for(i=0;i<lf->nbin-1;i++) {
    for(j=0;j<lf->nbin-1;j++) {
      if(DEBUG3) printf("   %g  ",covA[i][j]);  
    }
    if(DEBUG3) printf("\n");
  }
  gaussj_d(covA,lf->nbin-1,bb,1);
  gaussj_d(covB,lf->nbin-1,bb,1);
  gaussj_d(covC,lf->nbin-1,bb,1);
  if(DEBUG3) printf(" COVA despues gauss\n");
  for(i=0;i<lf->nbin-1;i++) {
    for(j=0;j<lf->nbin-1;j++) {
      if(DEBUG3) printf("   %g  ",covA[i][j]);  
    }
    if(DEBUG3) printf("\n");
  }


  if(DEBUG3) printf(" ya he hecho los gauus\n");
  
  /* Relleno la matriz de covarianza con las tres auxiliares */
  /* Uso covA para casi todo */
  for(i=0;i<lf->nbin-1;i++) {
    for(j=0;j<lf->nbin-1;j++) {
      covar[i][j]=covA[i][j];
    }
  } 
  /* Uso covB para la fila de abajo, la columna de la derecha y el extremo inferior derecha */
  for(j=1;j<lf->nbin;j++)       covar[lf->nbin-1][j]=covB[lf->nbin-2][j-1];
  for(i=1;i<lf->nbin;i++)       covar[i][lf->nbin-1]=covB[i-1][lf->nbin-2];
  /* Uso covC para el extremo superior derecha y el inferior izquierda  */
  covar[0][lf->nbin-1]=covC[0][lf->nbin-2];
  covar[lf->nbin-1][0]=covC[lf->nbin-2][0];

  if(DEBUG3) printf(" Termine de rellenar covar\n");


  
  /* Deshago el cambio para el limite de confidencia */
  for(i=0;i<lf->nbin;i++) {
    for(j=0;j<lf->nbin;j++) {
      covar[i][j]/=(-2*log(conflim));
      lf->covarlf[i][j]=covar[i][j];
   }
  }

  for(i=0;i<lf->nbin;i++) lf->errlf[i]=sqrt(covar[i][i]);
  for(i=0;i<lf->nbin;i++) parconf[i]=par[i];

  free_matrix_d(bb,lf->nbin-1,1);
  free(y);
  free(sigparconf);
  free(parconf);
  free(distmax);

  free_matrix_d(parelip,lf->nbin  ,nconflini);
  free_matrix_d(pareA  ,lf->nbin-1,nconflini);
  free_matrix_d(pareB  ,lf->nbin-1,nconflini);
  free_matrix_d(pareC  ,lf->nbin-1,nconflini);

  free_matrix_d(invcovar,lf->nbin  ,lf->nbin);
  free_matrix_d(   covar,lf->nbin  ,lf->nbin);
  free_matrix_d(invcovA ,lf->nbin-1,lf->nbin-1);
  free_matrix_d(invcovB ,lf->nbin-1,lf->nbin-1);
  free_matrix_d(invcovC ,lf->nbin-1,lf->nbin-1);
  free_matrix_d(   covA ,lf->nbin-1,lf->nbin-1);
  free_matrix_d(   covB ,lf->nbin-1,lf->nbin-1);
  free_matrix_d(   covC ,lf->nbin-1,lf->nbin-1);
  
  printf("\n");
}




 
void   ComputeNorma_SWML_M(int n, double mlim, double strrad, double zlow, double zup, struct Steplf_M *lf) {

  double Ntot;
  double pi=3.1415926535897932384;
  int j;

  Ntot=Int_step_M(*lf,zlow,zup,mlim,*co)*strrad/4./pi; 
  
  if(DEBUG) printf(" Normalacing Ntot %f\n",Ntot*4*pi/strrad);

  for(j=0;j<lf->nbin;j++) {
    if(DEBUG) printf(" From %f ",lf->lf[j]);
    lf->lf[j]=lf->lf[j]+log((float)(n)/Ntot);
    if(DEBUG) printf(" to %f \n",lf->lf[j]);
    /* El error, como es logaritmico, sigue intacto */
    lf->errlf[j]=lf->errlf[j];
  }

}
