#include <stdlib.h>
#include "random.h"

    
double Powdev(double xmin, double xmax,double alfa) {
  
  float A;
  float tmp;
  float dum;
  static long idum =-1;
  struct timeval tv;
  if(idum==-1) 
  {
    gettimeofday(&tv,NULL);
    idum=tv.tv_usec;
  }
  //if(idum==-1) idum=-(long)time(NULL)/2;
  do dum=ran2(&idum);
  while (dum == 0.0);

  A=(1+alfa)/(pow(xmax,1+alfa)-pow(xmin,1+alfa));

  tmp=dum*(1+alfa)/A+pow(xmin,1+alfa);
  tmp=pow(tmp,1./(1+alfa));

  return(tmp);
  
}
