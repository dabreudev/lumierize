#include "mrqmin.h"
#define DEBUG 0
#define DEBUG2 0
#define MAXALAMDA 1e6


int Mrq_d(double  x[],double y[],double sig[],int ndata,double a[],int ia[],int ma, double **covar, double *chisq, void (*funcs)(double, double [], double *, double [],int)) {
  
  double **alpha;
  double alamda;
  double oldalamda;
/*    int *ia; */
  double ochisq;
  double chisqTol = 0.0001;
  int iter=0;

  if(DEBUG) printf("MRQ: Enter %d\n",ma); 

  alpha=matrix_d(ma,ma);

  if(DEBUG)   printf(" iter INI %f %f %f \n",a[0],a[1],a[2]);  

  
  alamda=-1;
  mrqmin_d(x,y,sig,ndata, a,ia,ma,covar,alpha, chisq, funcs,&alamda);
  iter++;
  ochisq=*chisq*2;
  oldalamda=alamda;
  if(DEBUG) printf("iter 0 par %11.8g %11.8g %11.8g %g\n",a[0],a[1],a[2],*chisq);
  while((fabs(*chisq-ochisq)/ochisq>chisqTol  || alamda>oldalamda) ) {
    if(DEBUG) printf("iter  COND 1st %d   2nd %d   3rd %d   lam %g old %g\n",fabs(*chisq-ochisq)/ochisq>chisqTol,alamda>=0.001,alamda>=oldalamda,alamda,oldalamda);
    ochisq=*chisq;
    oldalamda=alamda;
    mrqmin_d(x,y,sig,ndata, a,ia,ma,covar,alpha, chisq, funcs,&alamda);
    if(DEBUG) printf("iter %d par %11.8g %11.8g chi %g lam %g ",iter,a[0],a[1],*chisq,alamda);
    iter++;
    if(DEBUG) printf(" ochisq %g chisq %g cond %g\n",ochisq,*chisq,fabs(*chisq-ochisq)/ochisq);
  }
/*  if(DEBUG)    printf("sal ite\n"); */
  alamda=0;
  mrqmin_d(x,y,sig,ndata, a,ia,ma,covar,alpha, chisq, funcs,&alamda);
  iter++;
  if(DEBUG) printf("iter %d END par %11.8g %11.8g %11.8g lam %g chi %g ochi %g c %g \n",iter,a[0],a[1],a[2],alamda,*chisq,ochisq,fabs(*chisq-ochisq)/ochisq);   
  return(iter);

  free_matrix_d(alpha,ma,ma);
}


void mrqmin_d(double  x[],double y[],double sig[],int ndata,double a[],int ia[],int ma, double **covar,double **alpha, double *chisq, void (*funcs)(double, double [], double *, double [],int), double *alamda) {
  int j,k,l,m;
  static int mfit;
  static double ochisq, *atry, *beta, *da, **oneda;
  if(*alamda < 0.0 ) {
    atry=vector_d(ma);
    beta=vector_d(ma);
    da  =vector_d(ma);
    
    for(mfit=0,j=1;j<=ma;j++) 
      if(ia[j-1]) mfit++;
    
    oneda=matrix_d(mfit,1);
    
    *alamda=0.001;
    /*      for(j=0;j<4;j++) printf(" fsolpar[%d] %f",j,a[j]);   */
    mrqcof_d(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);
    /*      printf(" chisq %g\n",*chisq);  */
    ochisq=(*chisq);
    for(j=1;j<=ma;j++) atry[j-1]=a[j-1];
  }
  if(DEBUG2) {
  for(j=1;j<=mfit;j++) {   
    for(k=1;k<=mfit;k++) {   
      printf(" alpha %d %d %g \n",j-1,k-1,alpha[j-1][k-1]);   
    }   
  }
  }
  for(j=0,l=1;l<=ma;l++) {
    if (ia[l-1]) {
      for(j++,k=0,m=1;m<=ma;m++) {
	if(ia[m-1]) {
	  k++;
	  covar[j-1][k-1]=alpha[j-1][k-1];
	}
      }
      if(DEBUG2) printf("covar %d %f\n",j,covar[j-1][j-1]); 
      covar[j-1][j-1]=alpha[j-1][j-1]*(1.0+(*alamda));
      oneda[j-1+(1-1)*mfit][0]=beta[j-1];
      if(DEBUG2) printf("covar %d %f\n",j,covar[j-1][j-1]); 
    }
  }
  if(DEBUG2) {
    for(j=1;j<=mfit;j++) {   
      for(k=1;k<=mfit;k++) {   
	printf(" %g ",covar[j-1][k-1]);   
      }   
      printf(" = %g\n",oneda[j-1][0]);   
    }   
    printf("antes gsuus\n");   
  }
  if(DEBUG2) printf("El mfit es %d\n",mfit);
  if(*alamda>=10000)  
    for(j=1;j<=mfit;j++) oneda[j-1][0]=oneda[j-1][0]/covar[j-1][j-1];
  else   gaussj_d(covar,mfit,oneda,1);
  if(DEBUG2) printf("sale gauus\n"); 
  for(j=1;j<=mfit;j++) da[j-1]=oneda[j-1][0];
  if(DEBUG2) printf(" da calculado\n");
  if (*alamda == 0.0 ) {
    if(DEBUG2) printf("Calculo covarianzas\n");  
    covsrt_d(covar,ma,ia,mfit);
    covsrt_d(alpha,ma,ia,mfit);
    if(DEBUG2) printf("Cova \n");  
    free_matrix_d(oneda,mfit,1);
    if(DEBUG2) printf("y ahora da\n");
    free(da);
    if(DEBUG2) printf("fr2 \n");  
    free(beta);
    if(DEBUG2) printf("ya\n");  
    free(atry); 
    if(DEBUG2) printf("sale\n");  
    return;
  }
  for(j=0,l=1;l<=ma;l++) 
    if(ia[l-1]){
      j++;
      atry[l-1]=a[l-1]+da[j -1];
      if(DEBUG2) printf("atry[%d] %f da %g\n",l-1,atry[l-1],da[j-1]);  
    }
  if(DEBUG2) {
    for(j=0;j<4;j++) printf(" fsolpar[%d] %f",j,atry[j]);   
    printf("\n");
  }
  mrqcof_d(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
  if(DEBUG2) printf(" chisq %f old %f\n",*chisq,ochisq); 
  if(DEBUG2) printf("second\n"); 
  if(*chisq < ochisq) {
    if(DEBUG2) printf("Decrem\n"); 
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
    if(DEBUG2)  printf("Increme\n"); 
    if(*alamda>MAXALAMDA);
    else *alamda *= 10.;
    *chisq=ochisq;
  }
  
  
}




void mrqcof_d(double x[],double y[],double sig[],int ndata,double a[],int ia[],int ma,double **alpha,double beta[],double *chisq,void (*funcs)(double ,double [],double *,double [],int)) {
  
  int i,j,k,l,m,mfit=0;
  double ymod,wt,sig2i,dy,*dyda;

  if(DEBUG2) {
    printf("mrqcof called with param: ");
    for(j=1;j<=ma;j++) printf(" %f ",a[j-1]);
    printf("\n");
  }

  dyda=vector_d(ma);

  for(j=1;j<=ma;j++)
    if(ia[j-1]) mfit++;
  if(DEBUG2)  printf(" mfit %d\n",mfit);
  for(j=1;j<=mfit;j++) {
    beta[j-1]=0.0;
    for(k=1;k<=j;k++) {
      alpha[j-1][k-1]=0.0;
      if(DEBUG2)       printf(" %d-%d al %g b %g\n",j,k,alpha[j-1][k-1],beta[j-1]); 
    }
  }
  *chisq=0.0;
  for(i=1;i<=ndata;i++) {
    if(DEBUG2) printf(" i %d x %f y %f sig %f\n",i,x[i-1],y[i-1],sig[i-1]);
    (*funcs)(x[i-1],a,&ymod,dyda,ma);
    if(DEBUG2) {
      printf(" ymod %f ",ymod);  
      for(j=1;j<=ma;j++) printf(" dyfun %g",dyda[j-1]);   
      printf("\n");  
    }
    sig2i=1.0/(sig[i-1]*sig[i-1]);
    dy=y[i-1]-ymod;
    if(DEBUG2) printf(" i %d y %f sigy %f dy %f sig2i %f\n",i,ymod,sig[i-1],dy,sig2i);  
    for(j=0,l=1;l<=ma;l++) {
      if(ia[l-1]) {
	wt=dyda[l-1]*sig2i;
	for(j++,k=0,m=1;m<=l;m++) {
	  if(ia[m-1]) {
	    k++;
	    alpha[j-1][k-1] += wt * dyda[m-1];
	    if(DEBUG2) printf(" al %g wt %g dyda %g\n",alpha[j-1][k-1],wt,dyda[m-1]); 
	  }
	}
	beta[j-1] += dy*wt;
	if(DEBUG2) printf(" beta %g\n",beta[j-1]);
      }
    }
    *chisq += dy*dy*sig2i;
    if(DEBUG2) printf(" x %f y %f ymod %f sig %f chi2 %g\n",x[i-1],y[i-1],ymod,sig[i-1],*chisq); 
  }
  for(j=2;j<=mfit;j++)
    for(k=1;k<j;k++) alpha[k-1][j-1]=alpha[j-1][k-1];
  if(DEBUG2) printf("and final chisq %f\n",*chisq);
  free(dyda);
  if(DEBUG2) printf(" Despues free\n");
}
	       



void covsrt_d(double **covar,int ma,int ia[],int mfit) {
  int i,j,k;
/*   double swap; */
  for(i=mfit+1;i<=ma;i++)
    for(j=1;j<=i;j++) covar[i-1][j-1]=covar[j-1][i-1]=0.0;
  k=mfit;
  for(j=ma;j>=1;j--) {
    if(ia[j-1]) {
      for( i=1;i<=ma;i++) SWAP_D(covar[i-1][k-1],covar[i-1][j-1]);
      for( i=1;i<=ma;i++) SWAP_D(covar[k-1][i-1],covar[j-1][i-1]);
      k--;
    }
  }
}
