 #include "modulos.h"
#define FTOL  1e-12
#define FTOL2 1e-6
#define FTOL3 1e-7
#define MAXITER  2000
#define MAXITER2 300
#define MAXTRIES   5
#define DEBUG  0
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

double Amoe_Funk_SWML_L_main(int n, double *x, double *y, double *p);
double Amoe_Funk_SWML_L_norm(int n, double *x, double *y, double *p);
double Amoe_Funk_SWML_L_prior(int n, double *x, double *y, double *p);
double Amoe_Funk_SWML_L_conf(int n, double *x, double *y, double *p);
double Funk2_int_SWML_L(double x);
double Funk1_norm_SWML_L(double x);
void   EmpiricalCovars_SWML_L(int n,double *flux,double *z,double *par, double *sigpar,double flim, struct cosmo_param cosmo,struct Steplf_L *lf);

void   ComputeNorma_SWML_L(int n, double fluxlim, double strrad, double zlow, double zup, struct Steplf_L *lf);
struct cosmo_param *co;
double fluxl;

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
double *lumbin;
int nbin;


int  MLA_SWML_L(int n,double *flux,double *z,double flim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Steplf_L *lf) {

  double *y;
  double *par;
  double *sigpar;
  int j;
  int iter_amo;
  double norm;
  int arezeros; 
  if(DEBUG) printf(" Estoy auiq SWML_L\n");

  /* Copias globales de estas variables para que se vean en toda la subrutina */
  ndata=n;
  co=&cosmo;
  fluxl=flim;
  nbin=lf->nbin;
  lumbin=lf->lumi;


  y=vector_d(n);
  par=vector_d(lf->nbin);
  sigpar=vector_d(lf->nbin);
  iter_m=0;
  iter_amo=0;


  VVmax_L(n,flux,z,flim,strrad,zlow,zup,cosmo,lf);
  if(DEBUG) for(j=0;j<nbin;j++) printf(" VVMAX Lum %g - %g LF %g\n",lf->lumi[j]/log(10),lf->lumi[j+1]/log(10),lf->lf[j]/log(10));
  arezeros=1;
  while(arezeros) {
    arezeros=0;
    for(j=0;j<nbin;j++) {
      if(lf->lf[j]==0 && lf->errlf[j]==0) {
	arezeros=1;
	if(j==0) {
	  lf->lf[0]=lf->lf[1];
	  lf->errlf[0]=lf->errlf[1];
	}
	else if(j==nbin-1) {
	  lf->lf[nbin-1]=lf->lf[nbin-2];
	  lf->errlf[nbin-1]=lf->errlf[nbin-2];
	}
	else {
	  lf->lf[j]=(lf->lf[j+1]+lf->lf[j-1])/2.;
	  lf->errlf[j]=(lf->errlf[j+1]+lf->errlf[j-1])/2.;
	}
      }
    }
  }
  if(DEBUG) for(j=0;j<nbin;j++) printf(" VVMAX Lum %g - %g LF %g\n",lf->lumi[j]/log(10),lf->lumi[j+1]/log(10),lf->lf[j]/log(10));
  printf(" Computing LF...\n");
  while(iter_amo==0) { 
    for(j=0;j<nbin;j++) {
      par[j]=lf->lf[j];
      sigpar[j]=4*lf->errlf[j];
    }
    norm=0;
    for(j=0;j<nbin;j++) norm+=exp(par[j])*(lumbin[j+1]-lumbin[j]);
    for(j=0;j<nbin;j++)   par[j]=par[j]-log(norm);
    iter_amo=Amoeba_d(n,flux,z,lf->nbin,par,sigpar,FTOL,MAXITER,Amoe_Funk_SWML_L_main);
  }

  MLmax=Amoe_Funk_SWML_L_main(n,flux,z,par);


  /* Meto la solucion en la salida */


  for(j=0;j<nbin;j++) lf->lf[j]=par[j];  


  /* Estimacion de los errores en los parametros */



  if(TRYEMPIRICAL) {
    if(DEBUG) printf(" llamantry \n");
    conflim=exp(-.5/10.);
    EmpiricalCovars_SWML_L(n,flux,z,par,sigpar,flim,cosmo,lf);  
    if(DEBUG) printf(" sale \n");


  }
  if(CONTOURPLOT) {
/*     nemp_f++; */
/*     conflim=exp(-.5*16.);     */
/*     EmpiricalCovars_SWML_L(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*9.);     */
/*     EmpiricalCovars_SWML_L(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*4.);     */
/*     EmpiricalCovars_SWML_L(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/1.);     */
/*     EmpiricalCovars_SWML_L(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/4.);     */
/*     EmpiricalCovars_SWML_L(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     cpgsci(1); */
/*     cpglab("P\\d1\\u","P\\d3\\u","Contornos de límites de confianza"); */
  }
  if(DEBUG) printf(" Calculo empirico\n");

  ComputeNorma_SWML_L(n, flim,  strrad, zlow,  zup, lf);

  
/*   conflim=readf(conflim); */

  
  free(y);
  free(par);
  free(sigpar);

  return(iter_amo);
}


double Amoe_Funk_SWML_L_main(int n, double *x, double *y, double *p) {

  int i,j; 
  int jini=0;
  double logL=0.;
  double Lumabs;
  double intstep;
  double Llow;
  double norm;
  double funl;

/*   for(j=0;j<nbin;j++) printf(" EN %f \n ",p[j]); */

  norm=0;
  for(j=0;j<nbin;j++) norm+=exp(p[j])*(exp(lumbin[j+1])-exp(lumbin[j]));
  for(j=0;j<nbin;j++) p[j]=p[j]-log(norm);

  if(DEBUG) {
    cpgsci(2);
    cpgsch(2.);
    for(j=0;j<nbin;j++) cpgpt1((float)((lumbin[j+1]+lumbin[j])/2./log(10)),(float)(log10(exp(p[j]))),5);
    if(DEBUG3) {
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
    if(jini!=-1   && jini!=nbin)   intstep+=exp(p[jini]+lumbin[jini])-exp(p[jini]+log(Llow));
    for(j=jini+1;j<nbin;j++) 	      intstep+=exp(p[j]+lumbin[j+1])-exp(p[j]+lumbin[j]);
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
/*     printf(" LF: lfamo.Mstar %g lfamo.phistar %f lf.alfa  %f\n",lfamo.Mstar,lfamo.phistar,lfamo.alfa); */
  }


  if(DEBUG) {
    cpgsci(0);
    cpgsch(2.);
    for(j=0;j<nbin;j++) cpgpt1((float)((lumbin[j+1]+lumbin[j])/2./log(10)),(float)(log10(exp(p[j]))),5);
    if(DEBUG3)  printf(" iter %d logL %f\n",iter_m,logL);
  }
  if(DEBUG)  printf(" iter %d logL %f\n",iter_m,logL);
  if(DEBUG) {
    printf("FL:\n");
    for(j=0;j<nbin;j++) printf("  L %f  LF %g %f ",(lumbin[j+1]+lumbin[j])/2./log(10),exp(p[j]),log10(exp(p[j])));
    printf("\n");
  }


  iter_m++;
  return(logL);
}


double Amoe_Funk_SWML_L_conf(int n, double *x, double *y, double *p) {

/*   printf(" MLM %g conf %f   %f %f\n",MLmax,conflim,fabs(Amoe_Funk_SWML_L_main(n,x,y,p)-(MLmax-log(conflim))),Amoe_Funk_SWML_L_main(n,x,y,p)); */
   
  return(fabs(Amoe_Funk_SWML_L_main(n,x,y,p)-(MLmax-log(conflim)))); 
}



void   EmpiricalCovars_SWML_L(int n,double *flux,double *z,double *par, double *sigpar,double flim, struct cosmo_param cosmo,struct Steplf_L *lf) {


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

  if(DEBUG) printf(" n vale %d \n",n);
  nconfl=NCONFL*lf->nbin;
  nconfl=2*lf->nbin*(lf->nbin+1)/2;
  nconflini=nconfl;

  bb=matrix_d(lf->nbin-1,1);
  for(i=0;i<lf->nbin-1;i++) (bb[i])[0]=0;
  y=vector_d(n);
  parconf=vector_d(lf->nbin);
  sigparconf=vector_d(lf->nbin);
  parelip=matrix_d(lf->nbin,  nconflini);
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
    iter_c=Amoeba_d(n,flux,z,lf->nbin,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_SWML_L_conf);
    if(DEBUG) printf(" %d SOL iter %d par0 %f    par1 %f\n", i,  iter_c ,parconf[0],   parconf[1]);
    for(j=0;j<lf->nbin;j++)  (parelip[j])[i]=parconf[j];
    if(iter_c==0 && Amoe_Funk_SWML_L_conf(n,flux,z,parconf)>FTOL2 ) {
      printf("\b.\b");
      i--;
      fflush(NULL);
    }
  }

  /* Supongo que el centro de la elipse es el valor que maximiza ML */
  for(i=0;i<nconfl;i++) {
    for(j=0;j<lf->nbin;j++)       (parelip[j])[i]-=par[j];
  }
  
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





void   ComputeNorma_SWML_L(int n, double flim, double strrad, double zlow, double zup, struct Steplf_L *lf) {

  double Ntot;
  double pi=3.1415926535897932384;
  int j;

  Ntot=Int_step_L(*lf,zlow,zup,flim,*co)*strrad/4./pi; 
  

  for(j=0;j<lf->nbin;j++) {
    if(DEBUG) printf(" From %f ",lf->lf[j]);
    lf->lf[j]=lf->lf[j]+log((float)(n)/Ntot);
    if(DEBUG) printf(" to %f \n",lf->lf[j]);
    /* El error, como es logaritmico, sigue intacto */
    lf->errlf[j]=lf->errlf[j];
  }

}
