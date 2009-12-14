#include "modulos.h"

#define DEBUG 0
#define DEBUG2 0 

static double (*amofunc_mlerr_amo_d_main_ptr)(int, double *, double *, double *);
static double amofunc_mlerr_amo_d_conf(int n, double *x, double *y, double *p);
static double  MLmax_mlerr_amo_d; 
static double conflim_mlerr_amo_d;

int mlerr_amo_d(int npt, double *x, double *y, int npar, double *par, double *sigpar, double ftol, int niter, double (*amofunc_main)(int, double *, double *, double *), int nconfl, double **covarpar) {
  
  
  int i,j;
  double *parconf;
  double *sigparconf; 
  double **invcovar;
  double **covar;
  double **parelip;
  double **bb;
  int iter_c;
  double first, median, third, *distmax;
  double *scale;
  int nconfmax;
  int iconf=0;


  nconfmax=4*nconfl;
  amofunc_mlerr_amo_d_main_ptr=amofunc_main;
  MLmax_mlerr_amo_d=amofunc_mlerr_amo_d_main_ptr(npt,x,y,par);
  conflim_mlerr_amo_d=exp(-.5/10.);
  
  bb=matrix_d(npar,1);
  parconf=vector_d(npar);
  sigparconf=vector_d(npar);
  parelip=matrix_d(npar,nconfl);
  invcovar=matrix_d(npar ,npar );
  covar=matrix_d(npar ,npar );
  distmax =vector_d(npar);
  scale =vector_d(npar);
  
  for(i=0;i<nconfl && iconf<=nconfmax;i++,iconf++) {
    for(j=0;j<npar;j++) {
/*       parconf[j]=par[j]+2*sigpar[j]*Gasdev(); */
      parconf[j]=par[j];
      sigparconf[j]=sigpar[j];
    }
    if(i>(int)(nconfl/2.)) {
      for(j=0;j<npar;j++) {
	parconf[j]=par[j]-((parelip[j])[(int)(i-nconfl/2.)+1]-par[j]);
        sigparconf[j]=((parelip[j])[(int)(i-nconfl/2.)+1]-par[j])/2.; 
      }
    }
    iter_c=0;
    iter_c=Amoeba_d(npt,x,y,npar,parconf,sigparconf,ftol,niter,amofunc_mlerr_amo_d_conf);
    if(DEBUG) {
      printf(" SOLCONF \n");
      for(j=0;j<npar;j++) printf(" par%d  %f",j,parconf[j]);
      printf(" EN iter %d ml %f\n",iter_c,amofunc_mlerr_amo_d_conf(npt,x,y,parconf));
    }
    for(j=0;j<npar;j++) (parelip[j])[i]=parconf[j];
    if(iter_c==0 || amofunc_mlerr_amo_d_conf(npt,x,y,parconf)>ftol )  {
      i--;
      if(DEBUG) printf(" DELETED\n");
    }
    if(DEBUG) printf(" i vale %d iconf %d/%d\n",i,iconf,nconfmax);
  }

  if(DEBUG) printf("  A este punto %d %d \n",i,iconf);

  nconfl=i;
 
  if(nconfl < (int)(npar*3))   return(0);


  if(DEBUG) {
    printf("\n Directamente\n");
    for(i=0;i<nconfl;i++) {
      for(j=0;j<npar;j++)   printf(" Par%d  %f",j,(parelip[j])[i]);
      printf("\n");
    }
  }
  for(i=0;i<nconfl;i++) {
    for(j=0;j<npar;j++)    (parelip[j])[i]-=par[j];
  }
  if(DEBUG) {
    printf("\n Ya he restado\n");
    for(i=0;i<nconfl;i++) {
      for(j=0;j<npar;j++)   printf(" Par%d  %f",j,(parelip[j])[i]);
      printf("\n");
    }
  }
  for(j=0;j<npar;j++) {
    Quartil_d(nconfl,parelip[j],&first,&median,&third);
    distmax[j]=maxf(fabs(first),fabs(third));
    if(DEBUG) printf(" Par %d first %f med %f thirf %f  distm %f\n",j,first,median,third,distmax[j]);
  }  
  for(i=0;i<nconfl;i++) {
    for(j=0;j<npar;j++) {
      if(fabs((parelip[j])[i])>3.5*2*distmax[j]/1.35) {
	for(j=0;j<npar;j++) {
	  memcpy(parelip[j]+i,parelip[j]+i+1,(nconfl-i-1)*sizeof(double));
	}
	i--;
	nconfl--;
	break;
      }
    }
  }
  for(j=0;j<npar;j++) {
    Quartil_d(nconfl,parelip[j],&first,&median,&third);
    scale[j]=(third-first);
    for(i=0;i<nconfl;i++) (parelip[j])[i]/=scale[j];
  }  
  if(DEBUG) {
    printf(" Despues scale\n");
    for(i=0;i<nconfl;i++) {
      for(j=0;j<npar;j++)   printf(" Par%d  %f",j,(parelip[j])[i]);
      printf("\n");
    }
  }
  MCElipN_d(nconfl,npar,parelip,invcovar);
  for(i=0;i<npar;i++) {
    for(j=0;j<npar;j++) {
      covar[i][j]=invcovar[i][j];
    }
  }
  gaussj_d(covar,npar,bb,1);
  /* Deshago el cambio para el limite de confidencia */
  for(i=0;i<npar;i++) {
    for(j=0;j<npar;j++) {
      covar[i][j]/=(-2*log(conflim_mlerr_amo_d))/(scale[i]*scale[j]);
      covarpar[i][j]=covar[i][j];
    }
  }
  
  
  free_matrix_d(bb,npar,1);
  free(parconf);
  free(sigparconf);
  free(distmax);
  free(scale);
  
  free_matrix_d(parelip,npar,nconfl);
  free_matrix_d(invcovar,npar ,npar);
  free_matrix_d(   covar,npar ,npar);
  
  return(1);
}


static double amofunc_mlerr_amo_d_conf(int n, double *x, double *y, double *p) {
  return(fabs(amofunc_mlerr_amo_d_main_ptr(n,x,y,p)-(MLmax_mlerr_amo_d-log(conflim_mlerr_amo_d)))); 
}


