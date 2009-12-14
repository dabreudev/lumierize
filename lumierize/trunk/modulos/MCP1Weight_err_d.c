#include "modulos.h"

void MCP1Weight_err_d(int n, double *x, double *y, double *sig,double *a, double *b, double *erra, double *errb,double *covab, double *chi2, float *q)

{
  int i;
  double sx=0.,sy=0.,sxx,sxy,ss;
  double wt,t,sxoss,st2=0.0;
/*    sx= StSuma1(n,x,1); */
/*    sy= StSuma1(n,y,1); */
/*    sxx=StSuma1(n,x,2); */
/*    sxy=StSuma2(n,x,1,y,1); */
/*    d=n*sxx-sx*sx; */
  
  *b=0;
  sx=0;sy=0;ss=0;sxx=0;sxy=0;
  for(i=0; i<n; i++) {
    wt=       1.0/sig[i]/sig[i];
    ss+=      wt;
    sx+=     x[i]*wt;
    sy+=     y[i]*wt;
    sxx+=     x[i]*x[i]*wt; 
    sxy+=     x[i]*y[i]*wt; 
/*      printf(" wt %f ss %f sx %f sy %f\n",wt,ss,sx,sy); */
  }
  sxoss=sx/ss;
  for(i=0; i<n; i++) {
    t=(x[i]-sxoss)/sig[i];
    st2 += t*t;
    *b += t*y[i]/sig[i];
/*      printf(" t %f st2 %f b %f\n",t,st2,*b); */

  }
  
  *b /= st2;
  *a=(sy-sx*(*b))/ss;
  *erra=sqrt((1.0+sx*sx/(ss*st2))/ss);
  *errb=sqrt(1.0/st2);
  *covab=-sx/(ss*st2);
  *chi2=0.0;
  
  for(i=0; i<n; i++) 
    *chi2 +=  (y[i]-(*a + *b * x[i])) *  (y[i]-(*a + *b * x[i])) / sig[i] / sig[i];
  if(n>2) *q=gammq(0.5*(n-2.0),0.5*(*chi2)); 
  else    *q=0;





  /*    sigdat=sqrt((*chi2)/(n-2.0)); */      /*Esto hay que aplicarlo si desconoces los errores */
  /*    *erra *= sigdat;  */
  /*    *errb *= sigdat;  */



}




