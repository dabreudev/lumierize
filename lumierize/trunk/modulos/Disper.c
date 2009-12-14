#include "modulos.h"


int Disper( float cosx, float cosy, float cocin, float alfa, float *cosxs, float *cosys ) 
{ 
  /* cosx: Es el coseno director del rayo de entrada segun el eje x 
     cosy: Es el coseno director del rayo de entrada segun el eje y 
     cocin= n(vidrio)/n(aire): Cociente entre indices de refraccion 
     alfa= Angulo del prisma. En radiames. 
     cosxs= Es el coseno director del rayo de salida segun el eje x 
     cosys= Es el coseno director del rayo de salida segun el eje y 
     */ 
  float pi=4*atan(1.); 
  double cosz,phiz,deltasup,deltain,deltasa,phi; 
  if(alfa==0.) {
    *cosxs=cosx;
    *cosys=cosy;
    return(0);
  }
  cosz= sqrt( 1- cosx*cosx -cosy*cosy); 
  phiz= acos(cosz); 
  deltasup=asin(1/cocin*sin(phiz)); 
  if (cosx==0) phi=pi/2; 
  else phi=atan(cosy/cosx); 
  if (cosx < 0)  phi=pi+phi ; 
  deltain=acos(cos(phi)*sin(deltasup)*sin(alfa)+cos(deltasup)*cos(alfa)); 
  if (cocin*sin(deltain)>1)  return(1);  
  else deltasa=asin(cocin*sin(deltain)); 
  *cosxs=((cos(deltasa)*cos(deltain)-cos(deltain-deltasa))*cos(phi)*sin(deltasup)+(cos(deltain-deltasa)*cos(deltain)-cos(deltasa))*sin(alfa))/(cos(deltain)*cos(deltain)-1); 
  *cosys=((cos(deltasa)*cos(deltain)-cos(deltain-deltasa))*sin(phi)*sin(deltasup))/(cos(deltain)*cos(deltain)-1); 
  if (*cosxs>=1. || *cosxs<=-1. || *cosys<=-1. || *cosys>=1.) return(1); 
/*   printf("cosx %f phi %f deltasup %f \n deltain %f deltasa %f\n cosxs %f cosys %f\n",cosx,phi,deltasup,deltain,deltasa,*cosxs,*cosys); */
  return(0); 
} 
