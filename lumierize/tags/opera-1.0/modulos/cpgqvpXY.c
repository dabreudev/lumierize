#include "modulos.h"
/* //#include "/usr/local/pgplot/cpgplot.h" */

void cpgqvpXY(float lxmin,  float lxmax,  float lymin,  float lymax,
              float xmin,   float xmax,   float ymin,   float ymax,
              float *vpxmin,float *vpxmax,float *vpymin,float *vpymax)

{

  float gxmin=0,gxmax=0,gymin=0,gymax=0,dvp;

  cpgqvsz(2,&gxmin,&gxmax,&gymin,&gymax);

  if(gxmax == gxmin || gymax == gymin) return;
  
  *vpxmin=lxmin;
  *vpxmax=lxmax;
  *vpymin=lymin;
  *vpymax=lymin+(lxmax-lxmin)/(xmax-xmin)*(ymax-ymin)*
               (gxmax-gxmin)/(gymax-gymin);

  if(*vpymax > lymax) {
    *vpymax=lymax;
    *vpxmax=lxmin+(lymax-lymin)/(ymax-ymin)*(xmax-xmin)*
                 (gymax-gymin)/(gxmax-gxmin);
    dvp=*vpxmax-*vpxmin;
    dvp=(lxmax-lxmin-dvp)/2.0;
    *vpxmin += dvp;
    *vpxmax += dvp;
    }
  else {
    dvp=*vpymax-*vpymin;
    dvp=(lymax-lymin-dvp)/2.0;
    *vpymin += dvp;
    *vpymax += dvp;
    }
}
