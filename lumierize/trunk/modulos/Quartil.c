#include "modulos.h"

int compare_quartil_f(const void *X1,const void *X2);

int compare_quartil_f(const void *X1,const void *X2)
{  
  float *x1,*x2; 
  x1=(float *)X1;
  x2=(float *)X2;
  if(*x1  < *x2) return(-1);
  if(*x1  > *x2) return(+1);
  if(*x1 == *x2) return(0);
  return(0);
}



void Quartil(int n,float *x,float *first, float *median,float *third)
{

  float *y;

  y=vector_f(n);
 
 
  memcpy(y,x,n*sizeof(float));
/*   for (i=0; i<n; i++) printf("QUAR %f   %f\n",x[i],y[i]);   */
  qsort(y,n,sizeof(float),compare_quartil_f);
/*   for (i=0; i<n; i++) printf("QUAR %f   %f\n",x[i],y[i]);   */

  *first=y[(int)(n/4)];
  *median=y[(int)(n/2)];
  *third=y[(int)(3*n/4)];
  free(y);
}
  
/* //Haciendo simulaciones, para una gaussiana: stddev= (third-first)/1.35 */
/* // no he visto que se corresponda con ningun numero raro. */
