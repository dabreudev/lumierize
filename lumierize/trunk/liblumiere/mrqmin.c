#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mrqmin.h"
#include "alloc.h"
#include "gaussj.h"

#define DEBUG 0
#define MAXALAMDA 1e6
#define SWAP(a,b) {float swap_temp=(a);(a)=(b);(b)=swap_temp;}

void mrqmin(float  x[],float y[],float sig[],int ndata,float a[],int ia[],int ma, float **covar,float **alpha, float *chisq, void (*funcs)(float, float [], float *, float [],int), float *alamda);
void mrqcof(float x[],float y[],float sig[],int ndata,float a[],int ia[],int ma,float **alpha,float beta[],float *chisq,void (*funcs)(float ,float [],float *,float [],int));
void covsrt(float **covar,int ma,int ia[],int mfit);

int Mrq(float  x[],float y[],float sig[],int ndata,float a[],int ia[],int ma, float **covar, float *chisq, void (*funcs)(float, float [], float *, float [],int)) {
  
  float **alpha;
  float alamda;
  float oldalamda;
/*    int *ia; */
  float ochisq;
  int iter=0;

  if(DEBUG) printf("MRQ: Enter %d\n",ma); 

  alpha=matrix_f(ma,ma);

  if(DEBUG)   printf(" iter INI %f %f %f \n",a[0],a[1],a[2]);  

  
  alamda=-1;
  if(DEBUG) printf("entra\n"); 
  mrqmin(x,y,sig,ndata, a,ia,ma,covar,alpha, chisq, funcs,&alamda);
  iter++;
  ochisq=*chisq*2;
  oldalamda=alamda;
  if(DEBUG) printf("iter 0 par %11.8g %11.8g %11.8g %g\n",a[0],a[1],a[2],*chisq);
  while((fabs(*chisq-ochisq)/ochisq>0.001  || alamda>oldalamda) ) {
    if(DEBUG) printf("iter  COND 1st %d   2nd %d   3rd %d   lam %g old %g\n",fabs(*chisq-ochisq)/ochisq>0.001,alamda>=0.001,alamda>=oldalamda,alamda,oldalamda);
    ochisq=*chisq;
    oldalamda=alamda;
    mrqmin(x,y,sig,ndata, a,ia,ma,covar,alpha, chisq, funcs,&alamda);
    if(DEBUG) printf("iter %d par %11.8g %11.8g %11.8g chi %g lam %g ",iter,a[0],a[1],a[2],*chisq,alamda);
    iter++;
    if(DEBUG) printf(" ochisq %g chisq %g cond %g\n",ochisq,*chisq,fabs(*chisq-ochisq)/ochisq);
  }
  if(DEBUG)    printf("sal ite\n"); 
  alamda=0;
  mrqmin(x,y,sig,ndata, a,ia,ma,covar,alpha, chisq, funcs,&alamda);
  iter++;
  if(DEBUG) printf("iter %d END par %11.8g %11.8g %11.8g lam %g chi %g ochi %g c %g \n",iter,a[0],a[1],a[2],alamda,*chisq,ochisq,fabs(*chisq-ochisq)/ochisq);   
  return(iter);

  free_matrix_f(alpha,ma,ma);
}


void mrqmin(float  x[],float y[],float sig[],int ndata,float a[],int ia[],int ma, float **covar,float **alpha, float *chisq, void (*funcs)(float, float [], float *, float [],int), float *alamda) {
  int j,k,l,m;
  static int mfit;
  static float ochisq, *atry, *beta, *da, **oneda;
  if(*alamda < 0.0 ) {
    atry=vector_f(ma);
    beta=vector_f(ma);
    da  =vector_f(ma);
    
    for(mfit=0,j=1;j<=ma;j++) 
      if(ia[j-1]) mfit++;
    
    oneda=matrix_f(mfit,1);
    
    *alamda=0.001;
    /*      for(j=0;j<4;j++) printf(" fsolpar[%d] %f",j,a[j]);   */
    mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);
    /*      printf(" chisq %g\n",*chisq);  */
    ochisq=(*chisq);
    for(j=1;j<=ma;j++) atry[j-1]=a[j-1];
  }
  for(j=0,l=1;l<=ma;l++) {
    if (ia[l-1]) {
      for(j++,k=0,m=1;m<=ma;m++) {
	if(ia[m-1]) {
	  k++;
	  covar[j-1][k-1]=alpha[j-1][k-1];
	}
      }
      if(DEBUG) printf("covar %d %f\n",j,covar[j-1][j-1]); 
      covar[j-1][j-1]=alpha[j-1][j-1]*(1.0+(*alamda));
      oneda[j-1+(1-1)*mfit][0]=beta[j-1];
      if(DEBUG) printf("covar %d %f\n",j,covar[j-1][j-1]); 
    }
  }
  if(DEBUG) {
    for(j=1;j<=mfit;j++) {   
      for(k=1;k<=mfit;k++) {   
	printf(" %g ",covar[j-1][k-1]);   
      }   
      printf(" = %g\n",oneda[j-1][0]);   
    }   
    printf("antes gsuus\n");   
  }
  if(DEBUG) printf("El mfit es %d\n",mfit);
  if(*alamda>=10000) 
    for(j=1;j<=mfit;j++) oneda[j-1][0]=oneda[j-1][0]/covar[j-1][j-1];
  else   gaussj(covar,mfit,oneda,1);
  if(DEBUG) printf("sale gauus\n"); 
  for(j=1;j<=mfit;j++) da[j-1]=oneda[j-1][0];
  if(DEBUG) printf(" da calculado\n");
  if (*alamda == 0.0 ) {
    if(DEBUG) printf("Calculo covarianzas\n");  
    covsrt(covar,ma,ia,mfit);
    covsrt(alpha,ma,ia,mfit);
    if(DEBUG) printf("Cova \n");  
    free_matrix_f(oneda,mfit,1);
    if(DEBUG) printf("y ahora da\n");
    free(da);
    if(DEBUG) printf("fr2 \n");  
    free(beta);
    if(DEBUG) printf("ya\n");  
    free(atry); 
    if(DEBUG) printf("sale\n");  
    return;
  }
  for(j=0,l=1;l<=ma;l++) 
    if(ia[l-1]){
      j++;
      atry[l-1]=a[l-1]+da[j -1];
      if(DEBUG) printf("atry[%d] %f da %g\n",l-1,atry[l-1],da[j-1]);  
    }
  if(DEBUG) {
    for(j=0;j<4;j++) printf(" fsolpar[%d] %f",j,atry[j]);   
    printf("\n");
  }
  mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
  if(DEBUG) printf(" chisq %f old %f\n",*chisq,ochisq); 
  if(DEBUG) printf("second\n"); 
  if(*chisq < ochisq) {
    if(DEBUG) printf("Decrem\n"); 
    if(*alamda<1.e-9);
    else *alamda *= 0.1;
    ochisq=(*chisq);
    for (j=0,l=1;l<=ma;l++) {
      if(ia[l-1]) {
	for(j++,k=0,m=1;m<=ma;m++) {
	  if(ia[m-1]) {
	    k++;
	    alpha[j-1][k-1]=covar[j-1][k-1];
	  }
	}
	beta[j-1]=da[j-1];
	a[l-1]=atry[l-1];
      }
    }
  } else {
    if(DEBUG)  printf("Increme\n"); 
    if(*alamda>MAXALAMDA);
    else *alamda *= 10.;
    *chisq=ochisq;
  }
  
  
}




void mrqcof(float x[],float y[],float sig[],int ndata,float a[],int ia[],int ma,float **alpha,float beta[],float *chisq,void (*funcs)(float ,float [],float *,float [],int)) {
  
  int i,j,k,l,m,mfit=0;
  float ymod,wt,sig2i,dy,*dyda;

  if(DEBUG) {
    printf("mrqcof called with param: ");
    for(j=1;j<=ma;j++) printf(" %f ",a[j-1]);
    printf("\n");
  }

 
  if((dyda=malloc(ma*sizeof(float))) == NULL) {
    printf("I cannot dimension dyda of %d elements \n",ma);
    exit(1);
  }

  for(j=1;j<=ma;j++)
    if(ia[j-1]) mfit++;
  if(DEBUG)  printf(" mfit %d\n",mfit);
  for(j=1;j<=mfit;j++) {
    beta[j-1]=0.0;
    for(k=1;k<=j;k++) {
      alpha[j-1][k-1]=0.0;
      if(DEBUG)       printf(" %d-%d al %g b %g\n",j,k,alpha[j-1][k-1],beta[j-1]); 
    }
  }
  *chisq=0.0;
  for(i=1;i<=ndata;i++) {
    if(DEBUG) printf(" i %d x %f y %f sig %f\n",i,x[i-1],y[i-1],sig[i-1]);
    (*funcs)(x[i-1],a,&ymod,dyda,ma);
    if(DEBUG) {
      printf(" ymod %f ",ymod);  
      for(j=1;j<=ma;j++) printf(" dyfun %g",dyda[j-1]);   
      printf("\n");  
    }
    sig2i=1.0/(sig[i-1]*sig[i-1]);
    dy=y[i-1]-ymod;
    if(DEBUG) printf(" i %d y %f sigy %f dy %f sig2i %f\n",i,ymod,sig[i-1],dy,sig2i);  
    for(j=0,l=1;l<=ma;l++) {
      if(ia[l-1]) {
	wt=dyda[l-1]*sig2i;
	for(j++,k=0,m=1;m<=l;m++) {
	  if(ia[m-1]) {
	    k++;
	    alpha[j-1][k-1] += wt * dyda[m-1];
	    if(DEBUG) printf(" al %g wt %g\n",alpha[j-1][k-1],wt); 
	  }
	}
	beta[j-1] += dy*wt;
	if(DEBUG) printf(" beta %g\n",beta[j-1]);
      }
    }
    *chisq += dy*dy*sig2i;
    if(DEBUG) printf(" x %f y %f ymod %f sig %f chi2 %g\n",x[i-1],y[i-1],ymod,sig[i-1],*chisq); 
  }
  for(j=2;j<=mfit;j++)
    for(k=1;k<j;k++) alpha[k-1][j-1]=alpha[j-1][k-1];
  if(DEBUG) printf("and final chisq %f\n",*chisq);
  free(dyda);
}
	       



void covsrt(float **covar,int ma,int ia[],int mfit) {
  int i,j,k;
/*   float swap; */
  for(i=mfit+1;i<=ma;i++)
    for(j=1;j<=i;j++) covar[i-1][j-1]=covar[j-1][i-1]=0.0;
  k=mfit;
  for(j=ma;j>=1;j--) {
    if(ia[j-1]) {
      for( i=1;i<=ma;i++) SWAP(covar[i-1][k-1],covar[i-1][j-1]);
      for( i=1;i<=ma;i++) SWAP(covar[k-1][i-1],covar[j-1][i-1]);
      k--;
    }
  }
}
