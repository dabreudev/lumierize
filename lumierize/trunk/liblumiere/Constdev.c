#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "random.h"

double Constdev(double xmin, double xmax) {
  static long idum =-1;
  struct timeval tv;
  if(idum==-1) 
  {
    gettimeofday(&tv,NULL);
    idum=tv.tv_usec;
  }
  //if(idum==-1) idum=-(long)time(NULL)/2;
  return(xmin+ran2(&idum)*(xmax-xmin));
}
