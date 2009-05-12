#include "modulos.h"

void cpgelip(float xc,float yc,float a,float b,float t)

{
  int i;
  float ct,st,p,cp,sp;
  float xx,yy,x[46],y[46];

  ct=cos(t/180*M_PI);
  st=sin(t/180*M_PI);

  for(i=0;i<46;i++) {
    p=2*M_PI*i/45.0;
    cp=cos(p);
    sp=sin(p);
    xx=a*cp;
    yy=b*sp;
    x[i]=xx*ct+yy*st+xc;
    y[i]=xx*st-yy*ct+yc;
    }

  cpgline(46,x,y);
}
