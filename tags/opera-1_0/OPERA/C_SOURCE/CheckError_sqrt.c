#include "modulos.h"

#define NTOT 200000
#define NTEO 1000
int main () 
{
  
  float x[NTOT];
  float y[NTOT];
  float z[NTOT];
  float r;
  int i;
  float media,sigma;
  float sx,sy;
/*   long idum=-(long)time(NULL)/2; */

  printf("Input sigma x: ");
  scanf("%f",&sx);
  printf("Input sigma y : ");
  scanf("%f",&sy);
  srandom((unsigned int)time(NULL)/2); 
  //cpgopen("/xserve");
  cpgopen("sqrt.ps/cps");
  cpgswin(-sx*3 ,sx*4,0.,NTOT/8.);
  //  cpgask(0);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);

  for(i=0;i<NTOT;i++) {
    x[i]=Gasdev()*sx;
    y[i]=Gasdev()*sy;
    z[i]=sqrt(x[i]*x[i]+y[i]*y[i]);
  }
  
  //  cpgline(200,x,y);
  cpghist(NTOT,x,-2*sx,2*sx,25,1);
  cpgsci(2);
  cpghist(NTOT,y,-2*sy,2*sy,25,1);
  cpgsci(3);
  cpghist(NTOT,z,0.,2*(sx+sy),25,1);
  cpgsci(4);
  for(i=0;i<NTEO;i++) {
    r=i*2*(sx+sy)/NTEO;
    cpgpt1(r,2*(sx+sy)*NTOT/25*r/sx/sx*exp(-r*r/sx/sx/2.),1);
  }

  cpgend(); 
  printf("\n\n");
  media= StMedia(NTOT,x,&sigma); 
  printf("BLANCO. Gaussiana.\n La media de x es %f y la desviacion %f \n\n",media,sigma); 

  media= StMedia(NTOT,y,&sigma); 
  printf("ROJO. Gaussiana.\n La media de y es %f y la desviacion %f \n\n",media,sigma); 

  media=StMedia(NTOT,z,&sigma); 
  printf("VERDE. z=sqrt(x*x+y*y).\n La media de z es %f y la desviacion %f \n\n",media,sigma); 

  printf(" AZUL: Curva teorica: z*exp(-z*z/2/sx/sy)*1/sx/sx\n La media es %f y la desviacion %f\n",sqrt(2*3.14156*sx)/2.,sx*sqrt((8-2*3.14156)/4));

  return 0;
}
