#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_machine.h>

#define DEBUG 0
#define DEBUG2 0
#define TOL 1.e-12
#define MAXIT 1000

void Covars_g_g(int n,float *x,float *errx,double mean,double sigma,float covar[2][2]);

int ML_g_g(int n,float *x,float *errx,float *mean,float *sigma,float *errmean,float *errsigma,float *covarmeansigma) {
  
  double sold,snew;
  double m;
  int i;
  int iter=0;

  int sigma_method=1;

  float meanerrx,stderrx; 

  float fnul;

  double tmp1,tmp2,tmp3;
  double tmp4,tmp5,tmp6;
  double tmp7,tmp8;
  double tmp9,tmp10;
  double tmp11,tmp12;
  double tmp13,tmp14;
  float covar[2][2];

  m=StErrWeightMedia(n,x,errx,&fnul); 
  snew=fnul;
  meanerrx=StMedia(n,errx,&stderrx);

  do {
    if(DEBUG2) printf(" Antes sold %g snew %g \n",sold,snew);
    sold=snew;
    if(DEBUG2) printf(" Despues sold %g snew %g \n",sold,snew);
    tmp1=0;
    tmp2=0;
    tmp3=0;
    tmp4=0;
    tmp5=0;
    tmp6=0;
    tmp7=0;
    tmp8=0;
    tmp9=0;
    tmp10=0;
    tmp11=0;
    tmp12=0;
    tmp13=0;
    tmp14=0;
    for(i=0;i<n;i++) {
      tmp1+=x[i]/errx[i]/errx[i]/(1+sold*sold/errx[i]/errx[i]);
      tmp2+=1./errx[i]/errx[i]/(1+sold*sold/errx[i]/errx[i]);
      tmp4+=x[i]/(1+errx[i]*errx[i]/sold/sold);
      tmp5+=1./(1+sold*sold/errx[i]/errx[i]);
      tmp7+=x[i]/(sold*sold+errx[i]*errx[i]);
      tmp8+=1/(sold*sold+errx[i]*errx[i]);
      tmp10+=1./(1.+errx[i]*errx[i]/sold/sold);
      tmp12+=1./(errx[i]*errx[i] + sold*sold);
      tmp13+=1./(errx[i]*errx[i]);
      tmp14+=1./(errx[i]*errx[i]*errx[i]*errx[i]);
    }
    for(i=0;i<n;i++) {
      tmp6+=(x[i]/errx[i]/errx[i]+m/sold/sold)*(x[i]/errx[i]/errx[i]+m/sold/sold)/(1./errx[i]/errx[i]+1./sold/sold)/(1./errx[i]/errx[i]+1./sold/sold); 
      tmp6+=(-2*m*x[i]/errx[i]/errx[i]-2*m*m/sold/sold)/(1./errx[i]/errx[i]+1./sold/sold); 
      tmp3+=(x[i]-m)*(x[i]-m)/sold/sold/(1+errx[i]*errx[i]/sold/sold)/(1+errx[i]*errx[i]/sold/sold);
      tmp9+=(x[i]-m)*(x[i]-m)/(1+errx[i]*errx[i]/sold/sold)/(1+errx[i]*errx[i]/sold/sold);
      tmp11+=(x[i]-m)*(x[i]-m)/(errx[i]*errx[i]+sold*sold)/(errx[i]*errx[i]+sold*sold);
      if(DEBUG2) printf(" tmp9 %g errx[i] %g sold %g x[i]-m %g fact %g\n",tmp9,errx[i],sold, x[i]-m,(x[i]-m)*(x[i]-m)/(1+errx[i]*errx[i]/sold/sold)/(1+errx[i]*errx[i]/sold/sold));
    }
    
    /* Aqui pongo varias formulas, todas v�lidas, pero la mejor es la de tmp7,8,9,10 */
    
/*     snew=sqrt((n*m*m+tmp6)/((float)n-tmp5));    */
/*     snew=sqrt(tmp3/tmp2); */
/*     snew=sqrt(sold*sold*tmp11/tmp12);  */
/*     snew=sqrt((tmp13-tmp11)/tmp14);  */

/*     sigma_method=1; */
    if(sigma_method==1) snew=sqrt(tmp9/tmp10);   
    else                snew=sqrt((tmp13-tmp11)/tmp14);
/*     m=tmp4/((float)n-tmp5);   */
/*     m=tmp1/tmp2; */
    m=tmp7/tmp8;
    
    if(DEBUG2) printf(" New iter m= %.10g   m_A= %.10g   m_B= %.10g  m_C= %.10g\n",m,tmp7/tmp8,tmp1/tmp2,tmp4/((float)n-tmp5));
    if(DEBUG2) printf(" %4d   %1d s= %.10g  s2_A= %.10g  s2_B= %.10g s2_C= %.10g  s2_D= %.10g s2_E= %.10g\n",iter,sigma_method,snew,tmp9/tmp10,tmp3/tmp2,(n*m*m+tmp6)/((float)n-tmp5),sold*sold*tmp11/tmp12,(tmp13-tmp11)/tmp14);

    /* Correccion si el primer metodo estima valores muy bajos */
    if(snew<meanerrx/50.) sigma_method=5;
/*     if(snew>meanerrx/5.)    sigma_method=1; */
    if(DEBUG2) printf(" snew %g 10.*meanerrx %g\n",snew,10.*meanerrx);

    iter++;
  } while (fabs(sold-snew)/fabs(snew)> TOL && iter<MAXIT);
  
  *mean=m;
  *sigma=snew;

  Covars_g_g(n,x,errx,m,snew,covar);
  *errmean=sqrt(covar[0][0]);
  *errsigma=sqrt(covar[1][1]);
  *covarmeansigma=covar[0][1];

  if(iter>=MAXIT) return(2);
  return(0);

}  


void  Covars_g_g(int n,float *x,float *errx,double mean,double sigma,float covar[2][2]) {
  
  int i;
  double tmp1,tmp2,tmp3,tmp4;
  double hessian[2][2];
  double hdet;
  tmp1=0;
  tmp2=0;
  tmp3=0;
  tmp4=0;
  for(i=0;i<n;i++) {
    tmp1+=-1/(sigma*sigma+errx[i]*errx[i]);
    tmp2+=+2.*sigma*sigma*(+sigma*sigma+errx[i]*errx[i]-2*(x[i]-mean)*(x[i]-mean))/(sigma*sigma+errx[i]*errx[i])/(sigma*sigma+errx[i]*errx[i])/(sigma*sigma+errx[i]*errx[i]);
    tmp3+=+1./(sigma*sigma+errx[i]*errx[i])-1.*(x[i]-mean)*(x[i]-mean)/(sigma*sigma+errx[i]*errx[i])/(sigma*sigma+errx[i]*errx[i]);
    tmp4+=-2*sigma*(x[i]-mean)/(sigma*sigma+errx[i]*errx[i])/(sigma*sigma+errx[i]*errx[i]);
  }

  if(DEBUG) printf(" tmp3 %g tmp2  %g sigma %f mean %f\n",tmp3,tmp2,sigma,mean);

  hessian[0][0]=tmp1;
  hessian[1][1]=tmp2;
  hessian[0][1]=tmp4;
  hessian[1][0]=tmp4;
  if(DEBUG) printf(" h[0][0] %f \n",hessian[0][0]);
  if(DEBUG) printf(" h[1][1] %f \n",hessian[1][1]);
  if(DEBUG) printf(" h[0][1] %f \n",hessian[0][1]);
  if(DEBUG) printf(" h[1][0] %f \n",hessian[1][0]);

  hdet=hessian[0][0]*hessian[1][1]-hessian[1][0]*hessian[1][0];

  if(DEBUG) printf(" hdet %f\n",hdet);

  covar[0][0]=-hessian[1][1]/hdet;
  covar[1][1]=-hessian[0][0]/hdet;
  covar[0][1]=+hessian[0][1]/hdet;
  covar[1][0]=+hessian[1][0]/hdet;

  if(DEBUG) printf(" covar[0][0] %f \n",covar[0][0]);

}

