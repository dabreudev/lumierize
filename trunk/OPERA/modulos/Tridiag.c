/*Este programa sirve para resolver un sistema de ecuaciones tridiagonal

El sistema es del tipo:
	b1*x1+c1*x2                  =r1
        a2*x1+b2*x2+c2*x3            =r2
	      a3*x2+b3*x3+c3*x4      =r3

			

			an*x(n-1)+bn*xn=rn

Todos los vectores a,b,c tienen dimension n
a1=0			!Cuidado!!
cn=0
La solucion es el vector u
 
Para compilar:
cc -c $s2/Proced/C/modulos/Tridiag.c


*/


#include "modulos.h"


int Tridiag(float *a,float *b,float *c,float *r,float *u,int n)
{
float *gam,bet;
int i;


if(n == 1) {
	 u[0]=r[0]/b[0];
	return(0);
	}

if((gam=malloc(n*sizeof(float))) == NULL) {
    printf("Tridiag: ERROR. No puedo dimensionar gam de %d bytes\n",
           n*sizeof(float));
    exit(1);
    }
if (b[0] == 0 ) return(1);
bet=b[0];
u[0]=r[0]/bet;
for(i=1;i<n;i++) {
	gam[i]=c[i-1]/bet;
	bet=b[i]-a[i]*gam[i];
	if (bet  == 0 ) return(1);
	u[i]=(r[i]-a[i]*u[i-1])/bet;

}
for(i=n-2;i>-1;i--) {   
	u[i]=u[i]-gam[i+1]*u[i+1];
	}
return(0);
}




