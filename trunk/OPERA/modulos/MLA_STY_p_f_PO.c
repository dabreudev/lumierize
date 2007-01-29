#include "modulos.h"
#define FTOL  1e-10
#define FTOL2 1e-6
#define FTOL3 1e-7
#define MAXITER  1000
#define MAXITER2 120
#define NBSNORMA  50
#define MAXTRIES   5
#define VERBOSE 0
#define DEBUG  1
#define DEBUG2 0
#define DEBUG3 0
#define DEBUGPLOT 0
#define CONTOURPLOT 0
#define NCONFL 18
#define TOLERR 0.07 

#ifdef TRYEMPIRICAL
#else 
#define TRYEMPIRICAL 0
#endif


/* Ahora mismo esta con el pcut definido, de modo que se aleja de valores 
   cercanos a 0. Puede afectar al calculo de las covarianzas!! */

double Amoe_Funk_STY_p_f_PO_main(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_p_f_PO_norm(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_p_f_PO_prior(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_p_f_PO_conf(int n, double *x, double *y, double *p);
double Funk2_int_STY_p_f_PO(double x);
double Funk1_norm_STY_p_f_PO(double x);
void   EmpiricalCovars_STY_p_f_PO(int n,double *magn,double *ew,double *z, double *par, double *sigpar, char photband[51], float gamma, float delta, float Kcoc, struct poselfunc fsel, struct SurveyItem si, double ewlim, struct Histdist ewd, struct cosmo_param cosmo,struct Schlf_L *lf);

double *ew_STY_p_f_PO;
struct cosmo_param *co_STY_p_f_PO;
double ewl_STY_p_f_PO;
struct Histdist *ewd_STY_p_f_PO;
struct poselfunc *fsel_STY_p_f_PO;
struct SurveyItem *si_STY_p_f_PO;

int ndata;
int nconfl;
double conflim;
double *pp;
int iter_m,iter_c;
double MLmax;
double sigi;
double xtmp;
double strrad_STY_p_f_PO;
char photband_STY_p_f_PO[51];
double gamma_STY_p_f_PO,delta_STY_p_f_PO,Kcoc_STY_p_f_PO;


int  MLA_STY_p_f_PO  (int n,double *magn,double *ew,double *z,char photband[51], float gamma, float delta, float Kcoc, struct poselfunc fsel, struct SurveyItem si, double ewlim, struct Histdist ewd, double strrad, struct cosmo_param cosmo, struct Schlf_L *lf) {


  double *y;
  double par[3];
  double sigpar[3];
  int i;
  int iter_amo;
  
  struct Steplf_L  lfvvmax;
  double minlum,maxlum;
  double *lumi;
  double *flux;
  double chisq;
  struct Schlf_L  lffit;
  double flim;
  double fluxmin,fluxmax;
  double ewmin,ewmax;
  double magnmin,magnmax;
  double ewmoda, magnmoda;
  double magnlim;
  

  lumi=vector_d(n);
  flux=vector_d(n);

  /* Copias globales de estas variables para que se vean en toda la subrutina */
  ndata=n;
  co_STY_p_f_PO=&cosmo;
  ewd_STY_p_f_PO=&ewd;
  fsel_STY_p_f_PO=&fsel;
  si_STY_p_f_PO=&si;
  ewl_STY_p_f_PO=ewlim;
  strrad_STY_p_f_PO=strrad;
  ew_STY_p_f_PO=ew;
  strcpy(photband_STY_p_f_PO,photband);
  gamma_STY_p_f_PO=gamma;
  delta_STY_p_f_PO=delta;
  Kcoc_STY_p_f_PO=Kcoc;

  y=vector_d(n);
  iter_m=0;

  if(DEBUG3) {
    for(i=0;i<n;i++) printf(" Entrada x %g\n",flux[i]);
  }

  lfvvmax.nbin=(int)(n/3);
  if(lfvvmax.nbin>30) lfvvmax.nbin=30;

  lfvvmax.lumi   =vector_d(lfvvmax.nbin+1);
  lfvvmax.errlumi=vector_d(lfvvmax.nbin+1);
  lfvvmax.lf     =vector_d(lfvvmax.nbin);
  lfvvmax.errlf  =vector_d(lfvvmax.nbin);
  lfvvmax.covarlf=matrix_d(lfvvmax.nbin,lfvvmax.nbin);

  for(i=0;i<n;i++)   {
    flux[i]=Flux_ew_mag(ew[i],magn[i],photband,gamma,delta,Kcoc);
    lumi[i]=log(Lum(z[i],flux[i],cosmo));
  }

  MinMax_d(n,flux,&fluxmin,&fluxmax);
  magnmin=0;
  magnmax=0;
  ewmin=0;
  ewmax=0;
  if(n>60) {
    magnmoda=StModa_d(n,magn,30,&magnmin,&magnmax);
    ewmoda= StModa_d(n,ew,30,&ewmin,&ewmax);
  }
  else {
    magnmoda=StModa_d(n,magn,(int)(n/2.),&magnmin,&magnmax);
    ewmoda= StModa_d(n,ew,(int)(n/2),&ewmin,&ewmax);
  }
  MinMax_d(n,magn,&magnmin,&magnmax);
  magnlim=(magnmax+magnmoda)/2;
  if(DEBUG) printf(" magnlim %f magnmin %f magnmoda %f\n",magnlim,magnmin,magnmoda);

  MinMax_d(n,lumi,&minlum,&maxlum);
  for(i=0;i<=lfvvmax.nbin;i++) lfvvmax.lumi[i]=minlum+i*(maxlum-minlum)/lfvvmax.nbin;
  flim=Flux_ew_mag(ewmoda,magnlim,photband,gamma,delta,Kcoc);
  VVmax_L(n,flux,z,flim,strrad,fsel.zbin[0],fsel.zbin[fsel.nz-1],cosmo,&lfvvmax);
/*   cpgopen("?"); */
/*   PlotStepLF_L(lfvvmax); */
  if(DEBUG) printf(" Salida VVmax\n");
  if(DEBUG) for(i=0;i<lfvvmax.nbin;i++) printf(" Lum %g - %g LF %g\n",lfvvmax.lumi[i]/log(10),lfvvmax.lumi[i+1]/log(10),lfvvmax.lf[i]/log(10));
  FitSch2StepLF_L(lfvvmax,&lffit, &chisq);
  if(DEBUG) {
    printf(" Desuess ajuste MRQ\n");
    printf(" Schechter FIT\n");
    printf(" log(Lstar (W)):  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",log10(lffit.Lstar),lffit.alfa,lffit.phistar,log10(lffit.phistar));
    printf(" E_log(Lstar):    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lffit.errLstar/lffit.Lstar/log(10.),lffit.erralfa,lffit.errphistar,lffit.errphistar/lffit.phistar/log(10.));
    printf(" Covar(Lstar,alpha): %g Covar(log(Lstar),alpha): %g\n",lffit.covaralfaLstar,lffit.covaralfaLstar/lffit.Lstar/log(10.));
    printf(" Covar(Lstar,Phistar): %g Covar(log(Lstar),log(Phistar)): %g\n",lffit.covarphistarLstar,lffit.covarphistarLstar/lffit.Lstar/log(10.)/lffit.phistar/log(10.));
    printf(" Covar(alpha,Phistar): %g Covar(alpha,log(Phistar)): %g\n",lffit.covaralfaphistar,lffit.covaralfaphistar/lffit.phistar/log(10.));
  }

/*   PlotStepSchLF_L(lfvvmax,lffit); */

/*   cpgclos(); */
  
  
  printf(" Computing LF...\n");

  iter_amo=0;
  while(iter_amo==0) { 
    par[0]=lffit.alfa;  /* Alpha */
    par[1]=log10(lffit.Lstar);
    par[2]=log(lffit.phistar);
/*     par[0]=-0.7; */
/*     par[1]=34; */
/*     par[2]=log(0.0033); */
    
    sigpar[0]=5.*lffit.erralfa;
    sigpar[1]=3.*lffit.errLstar/lffit.Lstar/log(10.);
    sigpar[2]=3.*lffit.errphistar/lffit.phistar;
    if(sigpar[0]>0.8) sigpar[0]=0.8;
    if(sigpar[1]>2) sigpar[1]=2;
    if(sigpar[1]<0.5) sigpar[1]=0.5;
    if(sigpar[2]>1) sigpar[2]=1;
    iter_amo=Amoeba_d(n,magn,z,3,par,sigpar,FTOL,MAXITER,Amoe_Funk_STY_p_f_PO_main);  
    if(DEBUG) printf(" iteramo %d\n",iter_amo);
    lf->alfa=par[0];
    lf->Lstar=pow(10.,par[1]);
    lf->phistar=exp(par[2]);
    if(DEBUG) {
      printf(" Solucion MALA\n");
      printf(" log(Lstar (W)):  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",log10(lf->Lstar),lf->alfa,lf->phistar,log10(lf->phistar));
      printf(" E_log(Lstar):    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lf->errLstar/lf->Lstar/log(10.),lf->erralfa,lf->errphistar,lf->errphistar/lf->phistar/log(10.));
    }
  } 
  MLmax=Amoe_Funk_STY_p_f_PO_main(n,magn,z,par);

  if(DEBUG) printf(" MLMAX %f\n",MLmax);

  /* Meto la solucion en la salida */
 
  lf->alfa=par[0];
  lf->Lstar=pow(10.,par[1]);
  lf->phistar=exp(par[2]);

  if(DEBUG) {
    printf(" Despoues ML  FIT\n");
    printf(" log(Lstar (W)):  %g   alpha:    %g  Phistar  :  %g (log=%g)\n",log10(lf->Lstar),lf->alfa,lf->phistar,log10(lf->phistar));
    printf(" E_log(Lstar):    %g   E_alpha:  %g  E_Phistar:  %g (log=%g) \n",lf->errLstar/lf->Lstar/log(10.),lf->erralfa,lf->errphistar,lf->errphistar/lf->phistar/log(10.));
  }

  /* Estimacion de los errores en los parametros */



  if(TRYEMPIRICAL) {
    if(DEBUG) printf(" llamantry \n");
    conflim=exp(-.5/10.);
    EmpiricalCovars_STY_p_f_PO(n,magn,ew,z,par,sigpar,photband,gamma,delta,Kcoc,fsel,si,ewlim,ewd,cosmo,lf); 
    if(DEBUG) printf(" sale \n");

    if(DEBUG) printf(" Solucion final: Lstar %.15g +/- %.15g alpha %.4g +/- %.4g\n",lf->Lstar,lf->errLstar,lf->alfa,lf->erralfa);

  }
  if(CONTOURPLOT) {
/*     nemp_f++; */
/*     conflim=exp(-.5*16.);     */
/*     EmpiricalCovars_STY_p_f_PO(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*9.);     */
/*     EmpiricalCovars_STY_p_f_PO(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5*4.);     */
/*     EmpiricalCovars_STY_p_f_PO(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/1.);     */
/*     EmpiricalCovars_STY_p_f_PO(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     conflim=exp(-.5/4.);     */
/*     EmpiricalCovars_STY_p_f_PO(n,x,errx,k,xk,Pk,sigpar,covar);  */
/*     nemp_f++; */
/*     cpgsci(1); */
/*     cpglab("P\\d1\\u","P\\d3\\u","Contornos de límites de confianza"); */
  }
  if(DEBUG) printf(" Calculo empirico\n");

  
/*   conflim=readf(conflim); */

  
  free(y);
  free(lumi);
  free(flux);
  free(lfvvmax.lumi);
  free(lfvvmax.errlumi);
  free(lfvvmax.lf);
  free(lfvvmax.errlf);
  free_matrix_d(lfvvmax.covarlf,lfvvmax.nbin,lfvvmax.nbin);
  return(iter_amo);
}


double Amoe_Funk_STY_p_f_PO_main(int n, double *x, double *y, double *p) {

  int i,j,k; 
  double logL=0.;
  struct Schlf_L lfamo;
  struct Schlf_M lfamo_M;
  double Lumi;
  double Mabs;
  double Llow;
  double dLdm;
  double xmin;
/*   double Mup=-30; */
/*   double Lup; */
  double intsup;
  double Ntot;
  double flim;
  double z,ew,magn,flux;
  double logprobarriba;
  double probabajo;
  double ewint,zint,Lint,magnint;
  int nEW=150;
  int nz=300;
  double normaewd;
  double pd;


  lfamo.alfa=p[0];
  lfamo.Lstar=pow(10.,p[1]);
  lfamo.phistar=exp(p[2]);

  lfamo_M.alfa=p[0];
  lfamo_M.Mstar=-2.5*p[1];
  lfamo_M.phistar=exp(p[2]);
  normaewd=0;
  for(i=0;i<ewd_STY_p_f_PO->k;i++) {
    if(ewl_STY_p_f_PO>ewd_STY_p_f_PO->xk[i+1]);
    else if(ewl_STY_p_f_PO>ewd_STY_p_f_PO->xk[i])  normaewd+=ewd_STY_p_f_PO->Pk[i]*(ewd_STY_p_f_PO->xk[i+1]-ewl_STY_p_f_PO);
    else normaewd+=ewd_STY_p_f_PO->Pk[i]*(ewd_STY_p_f_PO->xk[i+1]-ewd_STY_p_f_PO->xk[i]);
  }


  logL=0.;
  if(nEW<2*(ewd_STY_p_f_PO->k)) nEW=2*(ewd_STY_p_f_PO->k);

  intsup=incom(1+lfamo.alfa,200.);
  probabajo=Int_sch_f_PO(lfamo, photband_STY_p_f_PO,gamma_STY_p_f_PO, delta_STY_p_f_PO,Kcoc_STY_p_f_PO,*fsel_STY_p_f_PO,*si_STY_p_f_PO,ewl_STY_p_f_PO,*ewd_STY_p_f_PO,*co_STY_p_f_PO);
  if(probabajo==0) {
    printf(" MLA_STY_p_f_PO: probaabajo=0. Exiting\n");
    exit(1);
  }
  for(i=0;i<ndata;i++) {
    z=y[i];
    magn=x[i];
    ew=ew_STY_p_f_PO[i];
    flux=Flux_ew_mag(ew,magn,photband_STY_p_f_PO,gamma_STY_p_f_PO,delta_STY_p_f_PO,Kcoc_STY_p_f_PO);
    Lumi=Lum(z,flux,*co_STY_p_f_PO);
    Mabs=-2.5*log10(Lumi);
    flim=Flux_ew_mag(ewl_STY_p_f_PO,fsel_STY_p_f_PO->magbin[fsel_STY_p_f_PO->nmag-1],photband_STY_p_f_PO,gamma_STY_p_f_PO,delta_STY_p_f_PO,Kcoc_STY_p_f_PO);
    Llow=Lum(z,flim,*co_STY_p_f_PO);
    xmin=Llow/lfamo.Lstar;
      
    if(DEBUG3) printf(" ew %f Hist %f\n",ew,Histfunc(ew,*ewd_STY_p_f_PO));
    dLdm=Lumi*0.4*log(10); /* Esto es exacto */
    pd=prob_poselfunc_scale(*fsel_STY_p_f_PO,*si_STY_p_f_PO,ew,z,magn);
    if(DEBUG2) printf(" TESTFL _L %g _M %g\n",Schechter_L(Lumi,lfamo)+log(dLdm),log(Schechter_M(Mabs,lfamo_M)));
    logprobarriba=Schechter_L(Lumi,lfamo)+log(Histfunc(ew,*ewd_STY_p_f_PO)/normaewd)+log(dLdm)+log(dVdz(z,*co_STY_p_f_PO)/1e18)+log(pd);
    if(probabajo==0 || pd==0) logL+=100;
    else  logL-= logprobarriba - log(probabajo);
    if(DEBUG2) printf(" FACTORES Sch %f His %f dLdm %f dVdz %f \n",Schechter_L(Lumi,lfamo),log(Histfunc(ew,*ewd_STY_p_f_PO)),log(dLdm),log(dVdz(z,*co_STY_p_f_PO)));
    if(DEBUG2) printf(" obj %d logL %f logprobarriba %f probabajo %f oldint %f\n",i,logL,logprobarriba,probabajo,lfamo.phistar*(incom(1+lfamo.alfa,200)-incom(1+lfamo.alfa,xmin)));
  }
  /* Aqui viene la parte de la poissoniana de npob */
/*   Ntot=Int_sch_f_PO(lfamo, photband_STY_p_f_PO,gamma_STY_p_f_PO, delta_STY_p_f_PO,Kcoc_STY_p_f_PO,*fsel_STY_p_f_PO,*si_STY_p_f_PO,ewl_STY_p_f_PO,*ewd_STY_p_f_PO,*co_STY_p_f_PO); */
  Ntot=probabajo*strrad_STY_p_f_PO/4./M_PI; ;  /* Es lo mismo */
  logL-= (ndata*log(Ntot) - Ntot - gammln((double)ndata+1.));
  
  /*   printf(" Ntot %f\n",Ntot); */
  
  if(DEBUG2) printf(" NTOT %f ndata*log(Ntot) %f gamm %f\n",Ntot,ndata*log(Ntot),gammln((double)ndata)+1.);
  
  iter_m++;
  if(DEBUG) printf(" Iter %d  logL %f par0 %g par1 %g par2 %g\n",iter_m,logL,p[0],p[1],p[2]); 
/*   if(DEBUG) printf(" Iter %d  logL %f par0 %g par1 %g\n",iter_m,logL,p[0],p[1]); */
  return(logL);
}


double Amoe_Funk_STY_p_f_PO_conf(int n, double *x, double *y, double *p) {

/*   printf(" MLM %g conf %f   %f %f\n",MLmax,conflim,fabs(Amoe_Funk_STY_p_f_PO_main(n,x,y,p)-(MLmax-log(conflim))),Amoe_Funk_STY_p_f_PO_main(n,x,y,p)); */
   
  return(fabs(Amoe_Funk_STY_p_f_PO_main(n,x,y,p)-(MLmax-log(conflim)))); 
}



void   EmpiricalCovars_STY_p_f_PO(int n,double *magn,double *ew,double *z,double *par, double *sigpar,char photband[51], float gamma, float delta, float Kcoc, struct poselfunc fsel, struct SurveyItem si, double ewlim, struct Histdist ewd, struct cosmo_param cosmo,struct Schlf_L *lf) {

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
    iter_c=Amoeba_d(n,magn,z,3,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_STY_p_f_PO_conf);
    if(DEBUG) printf(" %d SOL iter %d par0 %f    par1 %f  par2 %f\n", i,  iter_c ,parconf[0],   parconf[1], parconf[2]);
    (parelip[0])[i]=parconf[0];
    (parelip[1])[i]=parconf[1];
    (parelip[2])[i]=parconf[2];
    if(iter_c==0 && Amoe_Funk_STY_p_f_PO_conf(n,magn,z,parconf)>FTOL2 ) {
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
  lf->errLstar=sqrt(covar[1][1])*lf->Lstar*log(10.);
  lf->errphistar=sqrt(covar[2][2])*lf->phistar; /* Ya que p[2] esta en logaritmos */
  lf->covaralfaLstar=covar[0][1]*lf->Lstar*log(10.);
  lf->covaralfaphistar=covar[0][2]*lf->phistar;  /* Por la misma razon */
  lf->covarphistarLstar=covar[1][2]*lf->phistar*lf->Lstar*log(10.);

  free_matrix_d(bb,3,1);
  free(y);
  free(parconf);
  free(sigparconf);

  free_matrix_d(parelip,3,nconflini);
  free_matrix_d(invcovar,3  ,3);
  free_matrix_d(   covar,3  ,3);

  printf("\n");
}


