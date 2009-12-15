#include "modulos.h"


void FitPol_deg2(float x1, float y1,float x2,float y2,float x3,float y3,struct Pol_deg2 *pol)  {

  float **a,**b;
  
  a=matrix_f(3,3);  
  b=matrix_f(3,1);
  

//  a[0][0]=x1*x1; a[0][1]=x2*x2; a[0][2]=x3*x3;
// a[1][0]=x1   ; a[1][1]=x2   ; a[1][2]=x3   ;
//  a[2][0]=1    ; a[2][1]=1    ; a[2][2]=1    ;
  a[0][0]=x1*x1; a[1][0]=x2*x2; a[2][0]=x3*x3;
  a[0][1]=x1   ; a[1][1]=x2   ; a[2][1]=x3   ;
  a[0][2]=1    ; a[1][2]=1    ; a[2][2]=1    ;
  b[0][0]=y1   ; b[1][0]=y2   ; b[2][0]=y3   ;

/*    Y resuelvo el sistema por Gauss */
/*    ------------------------------- */

  gaussj(a,3,b,1);

  pol->a=b[0][0];
  pol->b=b[1][0];
  pol->c=b[2][0];
  
  free_matrix_f(a,3,3);
  free_matrix_f(b,3,1);


}
