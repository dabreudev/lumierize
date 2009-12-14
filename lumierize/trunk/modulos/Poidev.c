#include "modulos.h"


/* La siguiente subrutina parece que hace tonterias para xm>1e7. Las desviaciones
   estandar empiezan a ser anormalmente altas */
	 
double Poidev(double xm)
{
/*   //float pi=3.141592654; */
  double g,t,y,sq,alxm;
  double em;
  static long idum =-1;
  struct timeval tv;
  if(idum==-1) 
  {
    gettimeofday(&tv,NULL);
    idum=tv.tv_usec;
  }
  //if(idum==-1) idum=-(long)time(NULL)/2;


  if(xm<12.0) {
    g=exp(-xm);
    em=-1;
    t=1.;
    do  {
      ++em;
/*       t=t*random()/2147483647.; */
      t=t*ran2(&idum);
    } while(t>g);
  }
  else {
    sq=sqrt(2.*xm);
    alxm=log(xm);
    g=xm*alxm-gammln(xm+1.);
/*     printf(" sq %f alxm %f g %f\n",sq,alxm,g); */
    do {
      do {
	/*       y=tan(M_PI*random()/2147483647.); */
	y=tan(M_PI*ran2(&idum));
	em=sq*y+xm;
/* 	printf(" em %f ",em); */
      } while (em<0.);
      em=floor(em);
/*       printf(" EM %f\n",em); */
      t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
/*       printf(" t %f\n",t); */
    } while(ran2(&idum)>t);
/*     } while(random()/2147483647.>t); */
  }
  return(em);
}








