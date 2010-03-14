#include <math.h>
#include "amoeba.h"
#include "mlerramoeba.h"
#include "alloc.h"

#define DEBUG 0
#define DEBUG2 0 

static float (*amofunc_mlerr_amo_main_ptr)(int, float *, float *, float *);
static float amofunc_mlerr_amo_conf(int n, float *x, float *y, float *p);
static float  MLmax_mlerr_amo; 
static float conflim_mlerr_amo;

int mlerr_amo(int npt, float *x, float *y, int npar, float *par, float *sigpar, float ftol, int niter, float (*amofunc_main)(int, float *, float *, float *), int nconfl, float **covarpar) {
  
  
  int i,j;
  float *parconf;
  float *sigparconf; 
  float **invcovar;
  float **covar;
  float **parelip;
  float **bb;
  int iter_c;
  float first, median, third, *distmax;
  float *scale;
  int nconfmax;
  int iconf=0;

  nconfmax=4*nconfl;
  amofunc_mlerr_amo_main_ptr=amofunc_main;
  MLmax_mlerr_amo=amofunc_mlerr_amo_main_ptr(npt,x,y,par);
  conflim_mlerr_amo=exp(-.5/10.);
  
  bb=matrix_f(npar,1);
  parconf=vector_f(npar);
  sigparconf=vector_f(npar);
  parelip=matrix_f(npar,nconfl);
  invcovar=matrix_f(npar ,npar );
  covar=matrix_f(npar ,npar );
  distmax =vector_f(npar);
  scale =vector_f(npar);
 
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
    iter_c=Amoeba(npt,x,y,npar,parconf,sigparconf,ftol,niter,amofunc_mlerr_amo_conf);
    if(DEBUG) {
      printf(" SOLCONF \n");
      for(j=0;j<npar;j++) printf(" par%d  %f",j,parconf[j]);
      printf(" EN iter %d ml %f\n",iter_c,amofunc_mlerr_amo_conf(npt,x,y,parconf));
    }
    for(j=0;j<npar;j++) (parelip[j])[i]=parconf[j];
    if(iter_c==0 || amofunc_mlerr_amo_conf(npt,x,y,parconf)>ftol )  {
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
    Quartil(nconfl,parelip[j],&first,&median,&third);
    distmax[j]=maxf(fabs(first),fabs(third));
    if(DEBUG) printf(" Par %d first %f med %f thirf %f  distm %f\n",j,first,median,third,distmax[j]);
  }  
  for(i=0;i<nconfl;i++) {
    for(j=0;j<npar;j++) {
      if(fabs((parelip[j])[i])>3.5*2*distmax[j]/1.35) {
	for(j=0;j<npar;j++) {
	  memcpy(parelip[j]+i,parelip[j]+i+1,(nconfl-i-1)*sizeof(float));
	}
	i--;
	nconfl--;
	break;
      }
    }
  }
  if(DEBUG) {
    printf(" Despues\n");
    for(i=0;i<nconfl;i++) {
      for(j=0;j<npar;j++)   printf(" Par%d  %f",j,(parelip[j])[i]);
      printf("\n");
    }
  }
  for(j=0;j<npar;j++) {
    Quartil(nconfl,parelip[j],&first,&median,&third);
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
  MCElipN(nconfl,npar,parelip,invcovar);
  for(i=0;i<npar;i++) {
    for(j=0;j<npar;j++) {
      covar[i][j]=invcovar[i][j];
      if(DEBUG) printf(" covar %f   ",covar[i][j]);
    }
    if(DEBUG)  printf(" = %f \n",bb[i][1]);
  }
  if(DEBUG) printf(" Despues MCElip\n");
  gaussj(covar,npar,bb,1);
  /* Deshago el cambio para el limite de confidencia */
  for(i=0;i<npar;i++) {
    for(j=0;j<npar;j++) {
      covar[i][j]/=(-2*log(conflim_mlerr_amo))/(scale[i]*scale[j]);
      covarpar[i][j]=covar[i][j];
      if(DEBUG) printf(" COVAR %f   ",covar[i][j]);
    }
    if(DEBUG)  printf("  \n");
  }
  
  
  
  free_matrix_f(bb,npar,1);
  free(parconf);
  free(sigparconf);
  free(distmax);
  free(scale);

  free_matrix_f(parelip,npar,nconfl);
  free_matrix_f(invcovar,npar ,npar);
  free_matrix_f(   covar,npar ,npar);

  return(1);
}


static float amofunc_mlerr_amo_conf(int n, float *x, float *y, float *p) {
  if(DEBUG2)   printf(" p0 %f p1 %f p2 %f MLmax %f   ML %f\n",p[0],p[1],p[2],MLmax_mlerr_amo,fabs(amofunc_mlerr_amo_main_ptr(n,x,y,p)-(MLmax_mlerr_amo-log(conflim_mlerr_amo)))); 
  return(fabs(amofunc_mlerr_amo_main_ptr(n,x,y,p)-(MLmax_mlerr_amo-log(conflim_mlerr_amo)))); 
}

