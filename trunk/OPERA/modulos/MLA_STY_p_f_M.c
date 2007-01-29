#include "modulos.h"
#define FTOL  1e-9
#define FTOL2 1e-6
#define FTOL3 1e-7
#define MAXITER  400
#define MAXITER2 120
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


double Amoe_Funk_STY_p_f_M_main(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_p_f_M_norm(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_p_f_M_prior(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_p_f_M_conf(int n, double *x, double *y, double *p);
void   EmpiricalCovars_STY_p_f_M(int n,double *magn, double *z, double *par, double *sigpar,struct fermifsel_M fsel, struct cosmo_param cosmo,struct Schlf_M *lf);

void   ComputeNorma_STY_p_f_M(int n, struct fermifsel_M fsel, double strrad, double zlow, double zup, struct Schlf_M *lf);

struct Schlf_M *lf_STY_p_f_M;
struct cosmo_param *co_STY_p_f_M;
double z_STY_p_f_M;
double magn_STY_p_f_M;
struct fermifsel_M fsel_STY_p_f_M;

int ndata;
int iter_m;
int iter_c;
int nconfl;
double conflim;
/* double ML[5*MAXTRIES*MAXITER]; */
double *pp;
double MLmax;
double zlow_STY_p_f_M,zup_STY_p_f_M,strrad_STY_p_f_M;


int  MLA_STY_p_f_M(int n,double *magn,double *z,struct fermifsel_M fsel, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf) {

  double *y;
  double par[3];
  double sigpar[3];
  int i;
  int iter_amo;

  struct Steplf_M  lfvvmax;
  double minMabs,maxMabs;
  double *Mabs;
  double chisq;
  struct Schlf_M  lffit;

  Mabs=vector_d(n);

  /* Copias globales de estas variables para que se vean en toda la subrutina */
  
  ndata=n;
  co_STY_p_f_M=&cosmo;
  fsel_STY_p_f_M=fsel;
  strrad_STY_p_f_M=strrad;
  zlow_STY_p_f_M=zlow;
  zup_STY_p_f_M=zup;

  y=vector_d(n);
  iter_m=0;

  if(DEBUG3) {
    for(i=0;i<n;i++) printf(" Entrada x %g\n",magn[i]);
  }

  lfvvmax.nbin=(int)(n/3);
  if(lfvvmax.nbin>30) lfvvmax.nbin=30;

  lfvvmax.magni   =vector_d(lfvvmax.nbin+1);
  lfvvmax.errmagni=vector_d(lfvvmax.nbin+1);
  lfvvmax.lf      =vector_d(lfvvmax.nbin);
  lfvvmax.errlf   =vector_d(lfvvmax.nbin);
  lfvvmax.covarlf =matrix_d(lfvvmax.nbin,lfvvmax.nbin);

  for(i=0;i<n;i++)   Mabs[i]=Mag(z[i],magn[i],cosmo);
  MinMax_d(n,Mabs,&minMabs,&maxMabs);
  for(i=0;i<=lfvvmax.nbin;i++) lfvvmax.magni[i]=minMabs+i*(maxMabs-minMabs)/lfvvmax.nbin;
  VVmax_M(n,magn,z,fsel_STY_p_f_M.magcut,strrad,zlow,zup,cosmo,&lfvvmax);
/*   cpgopen("?"); */
/*   PlotStepLF_M(lfvvmax); */
  if(DEBUG) printf(" Salida VVmax\n");
  if(DEBUG) for(i=0;i<lfvvmax.nbin;i++) printf(" Mabs %g - %g LF %g\n",lfvvmax.magni[i],lfvvmax.magni[i+1],lfvvmax.lf[i]/log(10));
  FitSch2StepLF_M(lfvvmax,&lffit, &chisq);
  if(DEBUG) {
    printf(" Desuess ajuste MRQ\n");
    printf(" Schechter FIT\n");
    printf(" Mstar (:  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",lffit.Mstar,lffit.alfa,lffit.phistar,log10(lffit.phistar));
    printf(" E_Mstar:    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lffit.errMstar,lffit.erralfa,lffit.errphistar,lffit.errphistar/lffit.phistar/log(10.));
    printf(" Covar(Lstar,alpha): %g Covar(Mstar,alpha): %g\n",lffit.covaralfaMstar,lffit.covaralfaMstar);
    printf(" Covar(Mstar,Phistar): %g Covar(Mstar,log(Phistar)): %g\n",lffit.covarphistarMstar,lffit.covarphistarMstar/lffit.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lffit.covaralfaphistar,lffit.covaralfaphistar/lffit.phistar/log(10.));
  }

/*   PlotStepSchLF_M(lfvvmax,lffit); */

/*   cpgclos(); */

  printf(" Computing LF...\n");

  iter_amo=0;
  while(iter_amo==0) { 
    par[0]=lffit.alfa;  /* Alpha */
/*     par[0]=-0.9; */
    par[1]=lffit.Mstar;
    par[2]=log(lffit.phistar);
/*     par[0]=-0.7; */
/*     par[1]=-21.5; */
/*     par[2]=log(0.0033); */
    

    sigpar[0]=6.*lffit.erralfa;
    sigpar[1]=6.*lffit.errMstar;
    sigpar[2]=6.*lffit.errphistar/lffit.phistar;
    if(sigpar[0]>0.8) sigpar[0]=0.8;
    if(sigpar[1]>3) sigpar[1]=3;
    if(sigpar[1]<0.5) sigpar[1]=0.5;
    if(sigpar[2]>1) sigpar[2]=1;
    iter_amo=Amoeba_d(n,magn,z,3,par,sigpar,FTOL,MAXITER,Amoe_Funk_STY_p_f_M_main);  
    if(DEBUG) printf(" iteramo %d\n",iter_amo);
    lf->alfa=par[0];
    lf->Mstar=par[1];
    lf->phistar=exp(par[2]);
    if(DEBUG) {
      printf(" Solucion MALA\n");
      printf(" Mstar :  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",lf->Mstar,lf->alfa,lf->phistar,log10(lf->phistar));
      printf(" E_log(Lstar):    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lf->errMstar,lf->erralfa,lf->errphistar,lf->errphistar/lf->phistar/log(10.));
    }
  }
  MLmax=Amoe_Funk_STY_p_f_M_main(n,magn,z,par);

  if(DEBUG) printf(" DESP MLMAX\n");
  /* Meto la solucion en la salida */
 
  lf->alfa=par[0];
  lf->Mstar=par[1];
  lf->phistar=exp(par[2]);

  /* Estimacion de los errores en los parametros */

  if(DEBUG) printf(" ANTES TRY\n");


  if(TRYEMPIRICAL) {
    if(DEBUG) printf(" llamantry \n");
    conflim=exp(-.5/10.);
    EmpiricalCovars_STY_p_f_M(n,magn,z,par,sigpar,fsel,cosmo,lf); 
    if(DEBUG) printf(" sale \n");

    if(DEBUG) printf(" Solucion final: Mstar %.15g +/- %.15g alpha %.4g +/- %.4g\n",lf->Mstar,lf->errMstar,lf->alfa,lf->erralfa);

  }

/*   ComputeNorma_STY_p_f_M(n, flim,  strrad, zlow,  zup, lf); */

  if(DEBUG) printf(" ANTES FREE\n");  

  free(y);
  free(Mabs);
  free(lfvvmax.magni);
  free(lfvvmax.errmagni);
  free(lfvvmax.lf);
  free(lfvvmax.errlf);
  free_matrix_d(lfvvmax.covarlf,lfvvmax.nbin,lfvvmax.nbin);
  return(iter_amo);
}


double Amoe_Funk_STY_p_f_M_main(int n, double *x, double *y, double *p) {

  int i; 
  double logL=0.;
  double logLold=0.;
  struct Schlf_M lfamo;
  double Mabs;
  double Lstar;
  double M;
  double magtmp,magmin;
  double Mleft,Mright;
  double xleft,xright;
/*   double Mup=-30; */
/*   double Lup; */
  double Ntot;
  double intabajo,inttemp;
  double logprobarriba;
  double probabajo;

  double intsup;
  int j;
  int nM_fs=200;

  lfamo.alfa=p[0];
  lfamo.Mstar=p[1];
  lfamo.phistar=exp(p[2]);

  lf_STY_p_f_M=&lfamo;

  logL=0.;
  logLold=0.;

  if(DEBUG3) printf(" Entra con  par0 %g par1 %g par2 %g\n",p[0],p[1],p[2]);

  Lstar=pow(10.,-0.4*lfamo.Mstar);
  intsup=incom(1+lfamo.alfa,200.);
  for(i=0;i<ndata;i++) {
    z_STY_p_f_M=y[i];
    Mabs=Mag(y[i],x[i],*co_STY_p_f_M);
    magmin=fsel_STY_p_f_M.magcut-5*fsel_STY_p_f_M.deltamag;
    Mleft=Mag(z_STY_p_f_M,magmin,*co_STY_p_f_M);
    Mright=Mag(z_STY_p_f_M,fsel_STY_p_f_M.magcut+5*fsel_STY_p_f_M.deltamag,*co_STY_p_f_M);
    xright=pow(10.,-0.4*Mright)/Lstar;
    xleft =pow(10.,-0.4*Mleft)/Lstar;
    if(DEBUG2) printf(" Mleft %f Mright %f Mabs %f  xleft %f xright %f\n",Mleft,Mright,Mabs,xleft,xright);
    if(xleft>=200) intsup=incom(1+lfamo.alfa,xleft+50.);
    else intsup=incom(1+lfamo.alfa,200.);
    if(xright>=200) {
      if(DEBUG2) printf(" KAKITA xright > 200\n");
      logL+=60; 
    }
    else {
      intabajo=0;
      if(Mleft<lfamo.Mstar && Mright>lfamo.Mstar) {
	for(j=0;j<nM_fs;j++) {
	  M=Mleft+j*(lfamo.Mstar-Mleft)/(nM_fs-1.);
	  magtmp=mag(z_STY_p_f_M,M,*co_STY_p_f_M);
	  intabajo=intabajo+Schechter_M(M,lfamo)*Fermi(magtmp,fsel_STY_p_f_M.magcut,fsel_STY_p_f_M.deltamag);
	}
	intabajo=intabajo/nM_fs*(lfamo.Mstar-Mleft);
	if(DEBUG2) printf(" Caso 1. Primera primera eta intabajo %g\n",intabajo);
	
	inttemp=0;
	for(j=0;j<nM_fs;j++) {
	  M=lfamo.Mstar+j*(Mright-lfamo.Mstar)/(nM_fs-1.);
	  magtmp=mag(z_STY_p_f_M,M,*co_STY_p_f_M);
	  inttemp=inttemp+Schechter_M(M,lfamo)*Fermi(magtmp,fsel_STY_p_f_M.magcut,fsel_STY_p_f_M.deltamag);
	}
	inttemp=inttemp/nM_fs*(Mright-lfamo.Mstar);  

	intabajo+=inttemp;
	if(DEBUG2) printf(" Caso 1. Primera segunda eta intabajo %g\n",intabajo);
      }
      else {
	/* Primero la integral a pelo de Lleft a Lright */
	for(j=0;j<nM_fs;j++) {
	  M=Mleft+j*(Mright-Mleft)/(nM_fs-1.);
	  magtmp=mag(z_STY_p_f_M,M,*co_STY_p_f_M);
	  intabajo=intabajo+Schechter_M(M,lfamo)*Fermi(magtmp,fsel_STY_p_f_M.magcut,fsel_STY_p_f_M.deltamag);
	}
	intabajo=intabajo/nM_fs*(Mright-Mleft);  
	if(DEBUG2) printf(" Primera etapa intabajo %g \n",intabajo);
      }
      /* Ahora la integral desde Mleft a -inf */
      if(xleft<200)  intabajo+=lfamo.phistar*(intsup-incom(1+lfamo.alfa,xleft));
      if(DEBUG2) printf(" Segunda etapa intabajo %g xleft %f  intsup %f\n",intabajo,xleft,intsup);
      probabajo=intabajo;
      
      logprobarriba=log(Schechter_M(Mabs,lfamo)) + log(Fermi(magn_STY_p_f_M,fsel_STY_p_f_M.magcut,-fsel_STY_p_f_M.deltamag));
      

      if(Fermi(magn_STY_p_f_M,fsel_STY_p_f_M.magcut,fsel_STY_p_f_M.deltamag)==0 || probabajo==0) logL+=50;
      else logL-= logprobarriba  - log(probabajo);   /* Perfectamente testado */
    }
    if(DEBUG2) printf(" iobj %d logL %f  probaba %f logprobarri %f\n",i,logL,probabajo,logprobarriba);
  }
  /* Aqui viene la parte de la poissoniana de npob */
  Ntot=Int_sch_f_M(lfamo,zlow_STY_p_f_M,zup_STY_p_f_M,fsel_STY_p_f_M,*co_STY_p_f_M)*strrad_STY_p_f_M/4./M_PI; 
  logL-=    (ndata*log(Ntot) - Ntot - gammln((double)ndata+1.)); 
/*   logLold-= (ndata*log(Ntot) - Ntot - gammln((double)ndata)+1.);  */
  
  if(DEBUG2) printf(" NTOT %f ndata*log(Ntot) %f gamm %f\n",Ntot,ndata*log(Ntot),gammln((double)ndata)+1.);

  iter_m++;
  if(DEBUG) printf(" Iter %d  logL %f logLold %f par0 %g par1 %g par2 %g\n",iter_m,logL,logLold,p[0],p[1],p[2]);
  return(logL); 
}


double Amoe_Funk_STY_p_f_M_conf(int n, double *x, double *y, double *p) {

/*   printf(" MLM %g conf %f   %f %f\n",MLmax,conflim,fabs(Amoe_Funk_STY_p_f_M_main(n,x,y,p)-(MLmax-log(conflim))),Amoe_Funk_STY_p_f_M_main(n,x,y,p)); */
   
  return(fabs(Amoe_Funk_STY_p_f_M_main(n,x,y,p)-(MLmax-log(conflim)))); 
}



void   EmpiricalCovars_STY_p_f_M(int n,double *magn,double *z,double *par, double *sigpar,struct fermifsel_M fsel, struct cosmo_param cosmo,struct Schlf_M *lf) {

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
    parconf[0]=par[0]+2*sigpar[0]*Gasdev();
    parconf[1]=par[1]+2*sigpar[1]*Gasdev();
    parconf[2]=par[2]+2*sigpar[2]*Gasdev();
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
    iter_c=Amoeba_d(n,magn,z,3,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_STY_p_f_M_conf);
    if(DEBUG) printf(" %d SOL iter %d par0 %f    par1 %f\n", i,  iter_c ,parconf[0],   parconf[1]);
    (parelip[0])[i]=parconf[0];
    (parelip[1])[i]=parconf[1];
    (parelip[2])[i]=parconf[2];
    if(iter_c==0 && Amoe_Funk_STY_p_f_M_conf(n,magn,z,parconf)>FTOL2 ) {
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
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      covar[i][j]=invcovar[i][j];
    }
  }
  gaussj_d(covar,3,bb,1);
  /* Deshago el cambio para el limite de confidencia */
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      covar[i][j]/=(-2*log(conflim));
    }
  }

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


