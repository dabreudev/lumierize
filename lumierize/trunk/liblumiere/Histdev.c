#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "alloc.h"
#include "random.h"

#define DEBUG 0

/* Pk debe estar dimensionado como k y xk como k+1 */
/* Pk es la probabilidad en el intervalo x[k],x[k+1] */


double Histdev(int k,double *Pk,double *xk)
{
  int i;
  double r;
  double x;
  static double *Pcumk;
  static double *Pkbuf;
  static long idum =-1;
  static int kbuf;
  int initcum=0;
  struct timeval tv;
  if(idum==-1) {
    gettimeofday(&tv,NULL);
    idum=tv.tv_usec;
    //idum=-(long)time(NULL)/2;
    Pcumk=vector_d(k);
    Pkbuf=vector_d(k);
    initcum=1;
    if(DEBUG) printf("ES IDEUM\n");
  }
  if(DEBUG) printf(" kbuf %d k %d\n",kbuf,k);
  if(kbuf!=k && initcum==0) {
    free(Pcumk);free(Pkbuf);
    Pcumk=vector_d(k);
    Pkbuf=vector_d(k);
    initcum=1;
    if(DEBUG)    printf(" ES AQUI\n");
  }
  kbuf=k;
  if(DEBUG)  printf(" initcum %d\n",initcum);
  if(DEBUG) for(i=0;i<k;i++) printf(" I %d Pkbuf %f Pk %f\n",i,Pkbuf[i],Pk[i]);
  if(memcmp(Pkbuf,Pk,k*sizeof(double))) initcum=1;
  if(initcum) {
    memcpy(Pkbuf,Pk,k*sizeof(double));
    Pcumk[0]=Pk[0]*(xk[1]-xk[0]);
    for(i=1;i<k;i++) {
      Pcumk[i]=Pk[i]*(xk[i+1]-xk[i])+Pcumk[i-1];
    }
    /* Por si acaso viene sin normalizar normalizo: */
    for(i=0;i<k;i++) {
      Pcumk[i]=Pcumk[i]/Pcumk[k-1];
      if(DEBUG)  printf(" k %d  P %f Pcumk %f [%f-%f]\n",i,Pk[i],Pcumk[i],xk[i],xk[i+1]);
    }
  }
  r=ran2(&idum); 
  if(DEBUG)  printf(" r %f\n",r);
  for(i=0;i<k;i++) {
    if(r<=Pcumk[i]) {
      if(DEBUG) printf(" Intervalo %d Pcum %f\n",i,Pcumk[i]);
      if(i==0) x=(r-0.)/(Pcumk[0]-0)*(xk[1]-xk[0])+xk[0];
      else     x=(r-Pcumk[i-1])/(Pcumk[i]-Pcumk[i-1])*(xk[i+1]-xk[i])+xk[i];
      if(DEBUG) printf(" Salida x= %f\n",x);
      return(x);
    }
  }
  if(DEBUG) printf(" Por aqui nunca debe pasar\n");
  return(xk[k]);
}
    
