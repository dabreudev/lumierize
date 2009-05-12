#include "modulos.h"

void cpgsvpXY(float lxmin,  float lxmax,  float lymin,  float lymax,
              float xmin,   float xmax,   float ymin,   float ymax)

{

  float vpxmin,vpxmax,vpymin,vpymax;
  
  cpgqvpXY(lxmin,lxmax,lymin,lymax,xmin,xmax,ymin,ymax,
           &vpxmin,&vpxmax,&vpymin,&vpymax);
/*   //printf("cpgsvpXY: VPORT(%f,%f,%f,%f)\n",vpxmin,vpxmax,vpymin,vpymax); */
  cpgsvp(vpxmin,vpxmax,vpymin,vpymax);
}
