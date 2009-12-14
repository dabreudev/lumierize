#include "modulos.h"


double Gasdev(void)
{
  double v1,v2;
  double r=2;
  static long idum =-1;
  struct timeval tv;
  if(idum==-1) 
  {
    gettimeofday(&tv,NULL);
    idum=tv.tv_usec;
  }
  //if(idum==-1) idum=-(long)time(NULL)/2;
/*   printf(" idum %d   \n",idum); */
/*   exit(1); */
  while(r>=1 || r == 0.0) {

/*     v1=2.*random()/2147483647.-1.; */
/*     v2=2.*random()/2147483647.-1.; */
     v1=2.*ran2(&idum)-1.; 
     v2=2.*ran2(&idum)-1.; 
    r=v1*v1+v2*v2;

  }

  return(v2*sqrt(-2.*log(r)/r));
}
    
