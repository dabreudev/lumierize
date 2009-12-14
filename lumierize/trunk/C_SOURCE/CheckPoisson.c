#include "modulos.h"
/* #define Lmax 7 */
#define NPT  100

#define NTOT 5000
int main () 
{
  
  float x[NTOT];
/*   float r; */
  int i,j;
  float media[NPT],sigma[NPT],coc[NPT];

  


  float lam=100;

  float lambda;
  
/*   float G; */

  float min,max;

  printf("Input parameter lambda: ");

  lam=readf(lam);

  for(j=0;j<NPT;j++) {
/*     j=94;  */
    lambda=lam/(NPT-1)*j+1;

    for(i=0;i<NTOT;i++) {
      x[i]=Poidev(lambda);
/*       //printf(" x[i] %f \n",x[i]); */
    }

    MinMax(NTOT,x,&min,&max);

    media[j]= StMedia(NTOT,x,sigma+j); 
    coc[j]=sigma[j]/sqrt(media[j]);
    printf(" %d lambda %f La media es %f y la desviacion %f  coc %f min %f max %f\n",j,lambda, media[j],sigma[j],sigma[j]/sqrt(media[j]),min,max); 
    media[j]=lambda;
/*      j=NPT+1;   */
  }
  cpgopen("?");
  cpgswin(0.,lam,0.9,1.5);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  
  cpgline(NPT,media,coc);

  return 0;
}
