#include "modulos.h"

void MCP1_err_d(int n, double *x, double *y,double *a, double *b, double *erra, double *errb,double *covab, double *chi2, float *q)

{
  int i;
  double sx=0.0,sy=0.0,sxx,sxy,ss;
  double t,sxoss,st2=0.0,sigdat;
/*    sx= StSuma1(n,x,1); */
/*    sy= StSuma1(n,y,1); */
/*    sxx=StSuma1(n,x,2); */
/*    sxy=StSuma2(n,x,1,y,1); */
/*    d=n*sxx-sx*sx; */
  
  *b=0;
  sx=0;sy=0;ss=0;sxx=0;sxy=0;
  for(i=0; i<n; i++) {
    sx+=     x[i];
    sy+=     y[i];
    sxx+=     x[i]*x[i]; 
    sxy+=     x[i]*y[i]; 
/*      printf(" wt %f ss %f sx %f sy %f\n",wt,ss,sx,sy); */
  }
  ss=n;

  sxoss=sx/ss;
  for(i=0; i<n; i++) {
    t=(x[i]-sxoss);
    st2 += t*t;
    *b += t*y[i];
/*     printf(" t %f st2 %f b %f\n",t,st2,*b);  */

  }
  
  *b /= st2;

  *a=(sy-sx*(*b))/ss;
  *erra=sqrt((1.0+sx*sx/(ss*st2))/ss);
  *errb=sqrt(1.0/st2);
  *covab=-sx/(ss*st2);
  *chi2=0.0;

  for(i=0; i<n; i++) 
    *chi2 +=  (y[i]-(*a + *b * x[i])) *  (y[i]-(*a + *b * x[i]));
/*    printf(" los %f %f de gam1q\n",0.5*(n-2.0),0.5*(*chi2)); */
  if(n>2)  *q=gammq(0.5*(n-2.0),0.5*(*chi2)); 
  else     *q=0;
  if(n>2)  sigdat=sqrt((*chi2)/(n-2.0));       /*Esto hay que aplicarlo si desconoces los errores */
  else     sigdat=0;
  *erra *= sigdat;   
  *errb *= sigdat;   
  *covab *= sigdat*sigdat;  

/*    printf(" Unos a %f b %f erra %f errb %f covab %f\n",*a,*b,*erra,*errb,*covab); */

/*    d=ss*sxx-sx*sx; */
/*    *a=(sy*sxx-sxy*sx)/d; */
/*    *b=(ss*sxy-sx*sy)/d; */
/*    *erra=sqrt(sxx/d); */
/*    *errb=sqrt(ss/d); */
/*    *covab=-sx/d; */

/*    printf(" Otros a %f b %f erra %f errb %f covab %f\n",*a,*b,*erra,*errb,*covab); */

}




