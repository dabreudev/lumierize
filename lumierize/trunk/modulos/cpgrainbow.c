#include "modulos.h"

void cpgrainbow(int *cimin,int *cimax)

{
  int i,ncol;
  float r[255],g[255],b[255];
  float dum,ang;

  for(i=0; i<255; i++) {
    r[i]=0.0;
    g[i]=0.0;
    b[i]=0.0;
    }

  cpgqcol(cimin,cimax);
  if(*cimin != 0) {
    printf("cpgrainbow: ERROR. cimin != 0 (%d)\n",*cimin);
    exit(1);
    }
  if(*cimax > 255) {
    printf("cpgrainbow: ERROR. cimin > 255 (%d)\n",*cimax);
    exit(1);
    }
  *cimin = 16;
  ncol = *cimax - *cimin + 1;

  for(i=0; i<ncol; i++) {
    ang=1.0+(5.50-1.0)*i/(ncol-1.);
    if(ang < M_PI &&
       (dum=sin( ang+M_PI/2.0)) > 0 ) r[i]=dum;
    if((dum=sin(-ang))          > 0 ) r[i]+=dum;
    if((dum=sin( ang-M_PI/2.0)) > 0 ) g[i]=dum;
    if((dum=sin( ang))          > 0 ) b[i]=dum;
    }
  for(i=0; i<ncol; i++) cpgscr( i + *cimin,r[i],g[i],b[i]);
  
}

