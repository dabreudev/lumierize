#include "modulos.h"

int compare_quartil_d(const void *X1,const void *X2);

int compare_quartil_d(const void *X1,const void *X2)
{  
  double *x1,*x2; 
  x1=(double *)X1;
  x2=(double *)X2;
  if(*x1  < *x2) return(-1);
  if(*x1  > *x2) return(+1);
  if(*x1 == *x2) return(0);
  return(0);
}



void Quartil_d(int n,double *x,double *first, double *median,double *third)
{

  double *y;

  y=vector_d(n);
 
 
  memcpy(y,x,n*sizeof(double));
  qsort(y,n,sizeof(double),compare_quartil_d);
  *first=y[(int)(n/4)];
  *median=y[(int)(n/2)];
  *third=y[(int)(3*n/4)];
  free(y);
}
  
/* //Haciendo simulaciones, para una gaussiana: stddev= (third-first)/1.35 */
/* // no he visto que se corresponda con ningun numero raro. */
