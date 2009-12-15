#include "modulos.h"

/*
  n       IN:  No. de datos
  *x      IN:  Matriz con las coordenadas X
  *y      IN:  Matriz con las coordenadas Y
  g       IN:  Grado del polinomio a ajustar
  *c      OUT: Matriz de los coeficientes del polinomio ajustado
  
  Dadas las matrices *x, *y con n datos, calcula el polinomio de grado g que
  mejor se ajusta a los datos mediante el algoritmo de mínimos cuadrados. 
  
  El polinomio ajustado es
  y = c[0] + c[1] x + c[2] x**2 + ... + c[g] x**g
  
  La subrutina devuelve la desviación est´ndar del ajuste 
  
  NOTA: La matriz de los coeficientes *c debe estar dimensionada a g+1 antes de la llamada a la subrutina 
  
  El programa aborta si: 
  
  g < 2 
  Aparece un error al dimensionar las matrices para el Sistema de Ecuaciones
  Aparece algún error al resolver el Sistema de Ecuaciones
  
*/


float MCPN(int n, float *x, float *y, int g, float *c)
     
{
  float *a,*b,d=0,s;
  int i,j;
  
  if(g<1) {
    printf("mcpn: ERROR. g=%d (g debe ser >=2)\n",g);
    exit(1);
  }
  
  if(g==1)  {
    MCP1(n, x, y, &c[0], &c[1]);
    goto desviacion;
  }
  
  
/*    Creo las matrices del sistema de ecuaciones */
/*    ------------------------------------------- */
  if( (a=malloc((g+1)*(g+1)*sizeof(float))) == NULL ) {
    printf("mcpn: ERROR. No puedo dimensionar la matriz a (g=%d)\n",g);
    exit(1);
  }
  if( (b=malloc((g+1)*sizeof(float))) == NULL ) {
    printf("mcpn: ERROR. No puedo dimensionar la matriz b (g=%d)\n",g);
    exit(1);
  }
  
  for(i=0; i<=g; i++) {
    b[i]=StSuma2(n,x,i,y,1);
    for(j=0; j<=g; j++) a[j+i*(g+1)]=StSuma1(n,x,i+j);
  }
  
/*    Y resuelvo el sistema por Gauss */
/*    ------------------------------- */
  if(!SELGauss(g+1,a,b,c)) {
    printf(" MCPN: No puedo resolver el sistema\n");
    exit(1);
  }
  
  free(a);
  free(b);

  desviacion:
  
/*    Calculo la desviacion estandar del ajuste */
/*    ----------------------------------------- */
  for(i=0; i<n; i++) {
    s=c[0];
    for(j=0; j<=g; j++) s += c[j]*powf(x[i],(float)j);
    s=y[i]-s;
    d += s*s;
  }
  d=sqrtf(d/(n-(g+1.0)));
  
  return(d);
  
}
