#include "modulos.h"


double Expdev(void)
{
  double dum;
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
  return -log(dum);
}
    
