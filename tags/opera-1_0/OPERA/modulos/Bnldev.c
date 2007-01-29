#include "modulos.h"
/* Subrutina de desviaciones binomiales. Cogida del Numerical Recipes.
   Suficientemente testada en todos los casos */

double Bnldev(double pp,int n)
{
  double p;
  double am;
  double bin;
  double g;
  double t,y,sq;
  double en,em;
  /*double pi=3.141592654;*/
  double pc,plog,pclog;
  int j;
  static long idum =-1;
  struct timeval tv;
  if(idum==-1) 
  {
    gettimeofday(&tv,NULL);
    idum=tv.tv_usec;
  }
  //if(idum==-1) idum=-(long)time(NULL)/2;



/*   printf(" DENTROOOOOOO      pp %f       n %d\n",pp,n); */
  if(pp<=0.5) p=pp;
  else p=1.-pp;
/*   printf(" p = %f\n",p); */
  am=n*p;
  
  if(n<25) {
/*     printf(" 111111111\n"); */
    bin=0;
/*     for(j=0;j<n;j++) if((random()/2147483647.)<p) bin=bin+1; */
    for(j=0;j<n;j++) if((ran2(&idum))<p) bin=bin+1;
  }
  else if (am<1.) {
    /* Quito esto porque equivale a una poissioniana
       
       g=exp(-am);
       t=1.;
       for(j=0;j<n+1;j++) {
       t=t*random()/2147483647.;
       if(t<g) goto L1;
       }
       j=n;
       L1:  bin=j;
    */
    bin=Poidev(am);
  }
  else {
/*     printf(" 3333333 pi %f\n",pi); */
    en=n;
    g=gammln(en+1.);
    pc=1.-p;
    plog=log(p);
/*     printf(" p = %f plog = %f\n",p,plog); */
    pclog=log(pc);
    sq=sqrt(2.*am*pc);
    do {
      do {
	y=tan(M_PI*ran2(&idum));
	/*    y=tan(M_PI*random()/2147483647.); */
	em=sq*y+am;
      } while(em<0 || em>=(en+1.));
/*       printf(" em %f\n",em); */
      em=(int)em;
/*       printf(" em %f\n",em); */
      t=1.2*sq*(1.+y*y)*exp(g-gammln(em+1.)-gammln(en-em+1.)+em*plog+(en-em)*pclog);
      /*     if((random()/2147483647.) > t) goto L2; */
    } while((ran2(&idum)) > t);
/*     } while((random()/2147483647) > t); */
    bin=em;
  }
  if(!(p==pp)) bin=n-bin;
  return(bin);
}


    
