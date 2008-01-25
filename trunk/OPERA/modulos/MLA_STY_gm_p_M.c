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
#define DEBUG2 1
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


double Amoe_Funk_STY_gm_p_M_main(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_gm_p_M_norm(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_gm_p_M_prior(int n, double *x, double *y, double *p);
double Amoe_Funk_STY_gm_p_M_conf(int n, double *x, double *y, double *p);
double Funk2_intmag_STY_gm_p_M(double fluxreal);
double Funk1_intMag_STY_gm_p_M(double x);
void   EmpiricalCovars_STY_gm_p_M(int n,double *magn,double *errmagn, double *z, double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf);


struct Schlf_M *_lf_STY_gm_p_M;
double *_errmagn_STY_gm_p_M;
struct cosmo_param *_cosmo_STY_gm_p_M;
double _mlim_STY_gm_p_M;
double _z_i_STY_gm_p_M;
double _magn_i_STY_gm_p_M;
double _errmagn_i_STY_gm_p_M;

int _ndata_STY_gm_p_M;
int _iter_m_STY_gm_p_M;
int _iter_c_STY_gm_p_M;
int _nconfl_STY_gm_p_M;
double _conflim_STY_gm_p_M;
/* double ML[5*MAXTRIES*MAXITER]; */
double *_pp_STY_gm_p_M;
double _MLmax_STY_gm_p_M;
double _zlow_STY_gm_p_M;
double _zup_STY_gm_p_M;
double _strrad_STY_gm_p_M;


int  MLA_STY_gm_p_M(int n,double *magn,double *errmagn,double *z,double mlim, double strrad, double zlow, double zup, struct cosmo_param cosmo,struct Schlf_M *lf, struct MLProcessInfo *mlinfo) {

  double par[3];
  double sigpar[3];
  int i;
  int iter_amo;

  /* Variables to use for vvmax fitting */
  struct Steplf_M  lfvvmax;
  double minMabs,maxMabs;
  double *Mabs;
  double chisq;
  struct Schlf_M  lffit;

  Mabs=vector_d(n);

  /* Copias globales de estas variables para que se vean en toda la subrutina */
  
  _errmagn_STY_gm_p_M=errmagn;
  _ndata_STY_gm_p_M=n;
  _cosmo_STY_gm_p_M=&cosmo;
  _mlim_STY_gm_p_M=mlim;
  _strrad_STY_gm_p_M=strrad;
  _zlow_STY_gm_p_M=zlow;
  _zup_STY_gm_p_M=zup;

  _iter_m_STY_gm_p_M=0;

  if(DEBUG3) {
    for(i=0;i<n;i++) printf(" Entrada x %g\n",magn[i]);
  }

  lfvvmax.nbin=(int)(n/3);
  if(lfvvmax.nbin>30) lfvvmax.nbin=30;

  lfvvmax.magni     =vector_d(lfvvmax.nbin+1);
  lfvvmax.errmagni  =vector_d(lfvvmax.nbin+1);
  lfvvmax.lnlf      =vector_d(lfvvmax.nbin);
  lfvvmax.errlnlf   =vector_d(lfvvmax.nbin);
  lfvvmax.covarlnlf =matrix_d(lfvvmax.nbin,lfvvmax.nbin);

  for(i=0;i<n;i++)   Mabs[i]=Mag(z[i],magn[i],cosmo);
  MinMax_d(n,Mabs,&minMabs,&maxMabs);
  for(i=0;i<=lfvvmax.nbin;i++) lfvvmax.magni[i]=minMabs+i*(maxMabs-minMabs)/lfvvmax.nbin;
  VVmax_M(n,magn,magn,z,mlim,strrad,zlow,zup,cosmo,&lfvvmax);
/*   cpgopen("?"); */
/*   PlotStepLF_M(lfvvmax); */
  if(DEBUG) printf(" Salida VVmax\n");
  if(DEBUG) for(i=0;i<lfvvmax.nbin;i++) printf(" Mabs %g - %g LF %g\n",lfvvmax.magni[i],lfvvmax.magni[i+1],lfvvmax.lnlf[i]/log(10));
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
  free(Mabs);
  free(lfvvmax.magni);
  free(lfvvmax.errmagni);
  free(lfvvmax.lnlf);
  free(lfvvmax.errlnlf);
  free_matrix_d(lfvvmax.covarlnlf,lfvvmax.nbin,lfvvmax.nbin);

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


    sigpar[0]=10.*lffit.erralfa;
    sigpar[1]=10.*lffit.errMstar;
    sigpar[2]=10.*lffit.errphistar/lffit.phistar;
    iter_amo=Amoeba_d(n,magn,z,3,par,sigpar,FTOL,MAXITER,Amoe_Funk_STY_gm_p_M_main);
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
  _MLmax_STY_gm_p_M=Amoe_Funk_STY_gm_p_M_main(n,magn,z,par);


  /* Meto la solucion en la salida */
 
  lf->alfa=par[0];
  lf->Mstar=par[1];
  lf->phistar=exp(par[2]);

  /* Estimacion de los errores en los parametros */



  if(TRYEMPIRICAL) {
    if(DEBUG) printf(" llamantry \n");
    _conflim_STY_gm_p_M=exp(-.5/10.);
    EmpiricalCovars_STY_gm_p_M(n,magn,errmagn,z,par,sigpar,mlim,cosmo,lf); 
    if(DEBUG) printf(" sale \n");

    if(DEBUG) printf(" Solucion final: Mstar %.15g +/- %.15g alpha %.4g +/- %.4g\n",lf->Mstar,lf->errMstar,lf->alfa,lf->erralfa);

  }

  return(iter_amo);
}


double Amoe_Funk_STY_gm_p_M_main(int n, double *x, double *y, double *p) {

  int i; 
  double logL=0.;
  double logLold=0.;
  struct Schlf_M lfamo;
  double Lstar;
  double Mabs;
  double Mlow;
  double xmin;
/*   double Mup=-30; */
/*   double Lup; */
  double Ntot;
  double offset;
  double scale;
  double x1,x2;
  int npa=21;
  int npb=21;
  int npc=31;
  double probarriba;
  double probabajo;

  double intsup;

  lfamo.alfa=p[0];
  lfamo.Mstar=p[1];
  lfamo.phistar=exp(p[2]);

  _lf_STY_gm_p_M=&lfamo;

  logL=0.;
  logLold=0.;

  Lstar=pow(10.,-0.4*lfamo.Mstar);
  
  if(DEBUG3) printf(" Entra con  par0 %g par1 %g par2 %g\n",p[0],p[1],p[2]);

  intsup=incom(1+lfamo.alfa,100.);
  for(i=0;i<_ndata_STY_gm_p_M;i++) 
  {
    _z_i_STY_gm_p_M=y[i];
    _magn_i_STY_gm_p_M=x[i];
    _errmagn_i_STY_gm_p_M=_errmagn_STY_gm_p_M[i];
    Mabs=Mag(y[i],x[i],*_cosmo_STY_gm_p_M);
    Mlow=Mag(y[i],_mlim_STY_gm_p_M,*_cosmo_STY_gm_p_M);
    xmin =pow(10.,-0.4*Mlow)/Lstar;
    if(xmin>=100)
    {
      logL+=10; 
    }
    else 
    {
      x1=lfamo.Mstar-5;
      x2=Mag(y[i],_mlim_STY_gm_p_M-6*_errmagn_i_STY_gm_p_M,*_cosmo_STY_gm_p_M);
      probabajo=gaussintleg_d(Funk1_intMag_STY_gm_p_M,x1,x2,npb);
      if(DEBUG2) printf(" Primer abajo %g con 1000 %g x1 %f x2 %f\n",probabajo,gaussintleg_d(Funk1_intMag_STY_gm_p_M,x1-10,x2,1000),x1,x2);
      x2=Mag(y[i],_mlim_STY_gm_p_M+6*_errmagn_i_STY_gm_p_M,*_cosmo_STY_gm_p_M);
      x1=Mag(y[i],_mlim_STY_gm_p_M-6*_errmagn_i_STY_gm_p_M,*_cosmo_STY_gm_p_M);
      probabajo+=gaussintleg_d(Funk1_intMag_STY_gm_p_M,x1,x2,npc);
      if(DEBUG2) printf(" Segundo abajo %g con 1000 %g x1 %f x2 %f\n",gaussintleg_d(Funk1_intMag_STY_gm_p_M,x1,x2,npc),gaussintleg_d(Funk1_intMag_STY_gm_p_M,x1,Mag(y[i],_mlim_STY_gm_p_M+12*_errmagn_i_STY_gm_p_M,*_cosmo_STY_gm_p_M),1000),x1,x2);

      if(DEBUG2) printf(" Calculo abajo %g old %g \n",probabajo,lfamo.phistar*(incom(1+lfamo.alfa,200.)-incom(1+lfamo.alfa,pow(10.,-0.4*Mlow)/Lstar)));

      /* We have to compute them again because Funk? may have overriden them */
      _z_i_STY_gm_p_M=y[i];
      _magn_i_STY_gm_p_M=x[i];
      _errmagn_i_STY_gm_p_M=_errmagn_STY_gm_p_M[i];

      offset=_magn_i_STY_gm_p_M;
      scale=_errmagn_i_STY_gm_p_M*sqrt(2.);
      x2=_mlim_STY_gm_p_M; 
      if(_mlim_STY_gm_p_M>(_magn_i_STY_gm_p_M+6*_errmagn_i_STY_gm_p_M)) x2=_magn_i_STY_gm_p_M+6*_errmagn_i_STY_gm_p_M; 
      x1=_magn_i_STY_gm_p_M-6*_errmagn_i_STY_gm_p_M; 
      if(_mlim_STY_gm_p_M<(_magn_i_STY_gm_p_M-5*_errmagn_i_STY_gm_p_M)) x1=x2-6*_errmagn_i_STY_gm_p_M;
      probarriba=gaussintleg_d(Funk2_intmag_STY_gm_p_M,x1,x2,npa);
      if(DEBUG2) printf(" parriba %g  con 1000 %g con her %g\n",probarriba,gaussintleg_d(Funk2_intmag_STY_gm_p_M,x1,x2,1000),gaussinther_d(Funk2_intmag_STY_gm_p_M,offset,scale,100));
      if(_errmagn_i_STY_gm_p_M==0) probarriba=Schechter_M(Mabs,*_lf_STY_gm_p_M);

      if(DEBUG2) printf(" Calculo arriba %g old %g x1 %g  x2 % g magn %g err %g magnl %g\n",probarriba,Schechter_M(Mabs,*_lf_STY_gm_p_M),x1,x2,_magn_i_STY_gm_p_M,_errmagn_i_STY_gm_p_M,_mlim_STY_gm_p_M);


/*       intsch=lfamo.phistar*(intsup-incom(1+lfamo.alfa,xmin)); */
/*       printf(" El de sche %g y el nuevo %g    con flux %g err %g  tanto %f\n",Schechter_M(Lumi,lfamo)+ log(dLumdflux(y[i],*co_STY_gm_p_M)),log(probarriba),flux_STY_gm_p_M,errflux_STY_gm_p_M,errflux_STY_gm_p_M/flux_STY_gm_p_M);  */
/*       printf(" La incom %g y el probabajo %g\n",log(intsch),log(probabajo)); */
/*       printf("   El de incom %g  y el nuevo %g\n",lfamo.phistar*(intsup-incom(1+lfamo.alfa,xmin)),probabajo); */

      /* Una vez comprobado que Schechter_M funciona bien, lo hago con esa */
      /*       logL-= log(probarriba)  -log(probabajo);  */    /* Esta hay que decomentarla */
/*       logLold-= Schechter_M(Lumi,lfamo)-log(intsch); */
      if(probarriba==0 || probabajo==0) logL+=10;
      logL-= log(probarriba)  -log(probabajo);   /* Perfectamente testado */
      if(DEBUG3) printf(" iobj %d logL %f loglold %f      sch %g    pa %g pb %g (%g)  xmin %g x1 %g x2 %g magn %g err %g\n",i,logL,logLold,Schechter_M(Mabs,lfamo),log(probarriba),log(probabajo),probabajo,xmin,x1,x2,_magn_i_STY_gm_p_M,_errmagn_i_STY_gm_p_M);
/*       printf(" iobj %d logL %f loglold %f      sch %g int %g (%g)   pa %g pb %g (%g)  xmin %g x1 %f x2 %f\n",i,logL,logLold,Schechter_M(Lumi,lfamo),log(intsch),intsch,log(probarriba),log(probabajo),probabajo,xmin,x1,x2); */
    }
  }
  /* Aqui viene la parte de la poissoniana de npob */
  Ntot=Int_sch_M(lfamo,_zlow_STY_gm_p_M,_zup_STY_gm_p_M,_mlim_STY_gm_p_M,*_cosmo_STY_gm_p_M)*_strrad_STY_gm_p_M/4./M_PI; 
  logL-=    (_ndata_STY_gm_p_M*log(Ntot) - Ntot - gammln((double)_ndata_STY_gm_p_M+1.)); 
/*   logLold-= (ndata*log(Ntot) - Ntot - gammln((double)ndata)+1.);  */
  
  if(DEBUG2) printf(" NTOT %f ndata*log(Ntot) %f gamm %f\n",Ntot,_ndata_STY_gm_p_M*log(Ntot),gammln((double)_ndata_STY_gm_p_M)+1.);

  _iter_m_STY_gm_p_M++;
  if(DEBUG) printf(" Iter %d  logL %f logLold %f par0 %g par1 %g par2 %g\n",_iter_m_STY_gm_p_M,logL,logLold,p[0],p[1],p[2]);
  return(logL); 
}


double Amoe_Funk_STY_gm_p_M_conf(int n, double *x, double *y, double *p) {

/*   printf(" MLM %g conf %f   %f %f\n",MLmax,conflim,fabs(Amoe_Funk_STY_gm_p_M_main(n,x,y,p)-(MLmax-log(conflim))),Amoe_Funk_STY_gm_p_M_main(n,x,y,p)); */
   
  return(fabs(Amoe_Funk_STY_gm_p_M_main(n,x,y,p)-(_MLmax_STY_gm_p_M-log(_conflim_STY_gm_p_M)))); 
}



void   EmpiricalCovars_STY_gm_p_M(int n,double *magn,double *errmagn,double *z,double *par, double *sigpar,double mlim, struct cosmo_param cosmo,struct Schlf_M *lf) {

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
    _iter_c_STY_gm_p_M=0;
    if(DEBUG) printf(" antes a nmo\n"); 
    _iter_c_STY_gm_p_M=Amoeba_d(n,magn,z,3,parconf,sigparconf,FTOL2,MAXITER2,Amoe_Funk_STY_gm_p_M_conf);
    if(DEBUG) printf(" %d SOL iter %d par0 %f    par1 %f\n", i,  _iter_c_STY_gm_p_M ,parconf[0],   parconf[1]);
    (parelip[0])[i]=parconf[0];
    (parelip[1])[i]=parconf[1];
    (parelip[2])[i]=parconf[2];
    if(_iter_c_STY_gm_p_M==0 && Amoe_Funk_STY_gm_p_M_conf(n,magn,z,parconf)>FTOL2 ) 
    {
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
      covar[i][j]/=(-2*log(_conflim_STY_gm_p_M));
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


double Funk1_intMag_STY_gm_p_M(double Mabs) 
{
  /* Lo suyo es que sea en magnitudes absolutas, ya que 
     la FL está expresada en intervalos de Mabs. 
     Si no, tendría que multiplicar por dMabs/dmag */

  int npa=21;
/*   double scale,offset; */
  double x1,x2;
  double firstint;

  _magn_i_STY_gm_p_M=mag(_z_i_STY_gm_p_M,Mabs,*_cosmo_STY_gm_p_M);
  x2=_mlim_STY_gm_p_M; 
  if(_mlim_STY_gm_p_M>_magn_i_STY_gm_p_M+6*_errmagn_i_STY_gm_p_M) x2=_magn_i_STY_gm_p_M+6*_errmagn_i_STY_gm_p_M; 
  x1=_magn_i_STY_gm_p_M-6*_errmagn_i_STY_gm_p_M; 
  if(_mlim_STY_gm_p_M<_magn_i_STY_gm_p_M-5*_errmagn_i_STY_gm_p_M) x1=x2-6*_errmagn_i_STY_gm_p_M;
  firstint=gaussintleg_d(Funk2_intmag_STY_gm_p_M,x1,x2,npa);
  if(_errmagn_i_STY_gm_p_M==0) firstint=Schechter_M(Mabs,*_lf_STY_gm_p_M);
/*   printf(" Segunda integral %g\n",firstint); */
/*   printf(" El bueno %g y el sch %g lum %g\n",log(firstint),Schechter_M(Lumi,*lf_STY_gm_p_M),lumlog);  */
/*   printf(" Sale con %g\n",firstint*Lumi); */
  return(firstint); /* Este producto es para tener en cuenta que hacemos la integral en log(flux) */
}

double Funk2_intmag_STY_gm_p_M(double magnreal) 
{
  double Mabs;
  double logfacLF,logfacerr;
  if(magnreal>_mlim_STY_gm_p_M) return(0);
  else 
  {
    Mabs=Mag(_z_i_STY_gm_p_M,magnreal,*_cosmo_STY_gm_p_M);
    logfacLF = log(Schechter_M(Mabs,*_lf_STY_gm_p_M));
    logfacerr= log(gaussian(magnreal,_magn_i_STY_gm_p_M,_errmagn_i_STY_gm_p_M));
    return(exp(logfacLF+logfacerr));
  }
}
