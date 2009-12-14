#include "modulos.h"

#define NTOT 5000
#define NBIN 50
int main () 
{
  
	int i;
	float media,sigma;
  float x[NTOT];
	float xx[1000],yy[1000];
   for(i=0;i<1000;i++) {
    xx[i]=i/100.;
    yy[i]=NTOT*exp(-xx[i])*10/NBIN;
	printf(" xx %f yy %f\n",xx[i],yy[i]);
  }


  for(i=0;i<NTOT;i++) {
    x[i]=Expdev();
  }
  
  cpgopen("?");
  cpgswin(0.,10.,0.,yy[0]*1.3);
  //  cpgask(0);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpgline(1000,xx,yy);
  cpghist(NTOT,x,0.,10.,NBIN,1);
  cpgend(); 
  printf("\n\n");
  printf("Distribucion inicial de x:\n"); 
  media= StMedia(NTOT,x,&sigma); 
  printf(" La media es %f y la desviacion %f  S/N %f\n\n",media,sigma,media/sigma); 

  return 0;

}
