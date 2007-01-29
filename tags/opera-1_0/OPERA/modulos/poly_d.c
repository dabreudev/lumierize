#include "modulos.h"

#define DEBUG 0

/*
  x       IN:  Valor de la variable x donde calcular el polinomio
  g       IN:  Grado del polinomio a ajustar
  *c      OUT: Matriz de los coeficientes del polinomio ajustado

  Dado el polinomio de coeficientes c, calcula el valor de y:
  y = c[0] + c[1] x + c[2] x**2 + ... + c[g] x**g
  
*/


double poly_d(double x,int g, double *c)
     
{
  double y=0;
  int j;
/*   //printf("x %f\n",x); */
  for(j=0; j<=g; j++) {
/*     //printf(" y %f\n",y); */
    y+= c[j]*powf(x,(float)j);
/*     //printf(" %d  y %f %f  cj %f\n",j,y,c[j]*powf(x,(float)j),c[j]); */
  }
/*   //exit(1); */
  return(y);  
}
