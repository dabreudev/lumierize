#include "modulos.h"

#define NTOT 15000
int main () 
{
  
  float x[NTOT];
  float y1[NTOT];
  float y2[NTOT];
  float y3[NTOT];
  float r;
  int i;
  float media,sigma;



  float lambda;

  float QE;

  printf("Input parameter lambda: ");
  scanf("%f",&lambda);
  printf("Input quantumm eficciency: ");
  scanf("%f",&QE);

  for(i=0;i<NTOT;i++) {
    x[i]=Poidev(lambda);
    r=random()/2147483647.;
    y1[i]=x[i]*QE;               /* // Dividiendo a pelo  */
    y2[i]=Bnldev(QE,x[i]);       /* // Combinacion de Poisson + Binomial */
    y3[i]=Poidev(lambda*QE);     /* // Poisson para electrones */
/*     //    printf(" Valor : %f  ",y2[i]); */
  }
  
  cpgopen("fot.ps/cps");
  cpgswin(lambda*QE-3*sqrt(lambda*QE),lambda*QE+3*sqrt(lambda*QE),0.,NTOT/2.);
/*   //  cpgask(0); */
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
/*   //  cpgline(200,x,y); */
  cpghist(NTOT,y1,lambda*QE-3*sqrt(lambda*QE),lambda*QE+3*sqrt(lambda*QE),25,1);
  cpgsci(2);
  cpghist(NTOT,y2,lambda*QE-3*sqrt(lambda*QE),lambda*QE+3*sqrt(lambda*QE),25,1);
  cpgsci(3);
  cpghist(NTOT,y3,lambda*QE-3*sqrt(lambda*QE),lambda*QE+3*sqrt(lambda*QE),25,1);

  cpgend(); 
  printf("\n\n");
  printf("Distribucion inicial de fotones x:\n"); 
  media= StMedia(NTOT,x,&sigma); 
  printf(" La media es %f y la desviacion %f  S/N %f\n\n",media,sigma,media/sigma); 

  printf(" BLANCO: Distribucion de y=x*QE. Multiplicacin normal:\n"); 
  media= StMedia(NTOT,y1,&sigma); 
  printf(" La media es %f y la desviacion %f  S/N %f\n\n",media,sigma,media/sigma); 

  printf(" ROJO: Distribucion de la combinacion de Poisson(x)+Binomial(QE):\n");
  media=StMedia(NTOT,y2,&sigma); 
  printf(" La media es %f y la desviacion %f  S/N %f\n\n",media,sigma,media/sigma); 

  printf(" VERDE: Distribucion de Poisson(y) para electrones:\n");
  media=StMedia(NTOT,y3,&sigma); 
  printf(" La media es %f y la desviacion %f  S/N %f\n\n",media,sigma,media/sigma); 
  return(0);
}
