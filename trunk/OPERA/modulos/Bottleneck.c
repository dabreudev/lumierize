#include "modulos.h"



void Bottleneck(float *esp,int nldo,int ntot,float *evec) {
  
  float *classes;
  float norma=0;
  int ii;
  int i=0,j;
  int nclasses;

/*   float *pc,*plc; */
  
  if((classes=malloc(nldo*ntot*sizeof(float)))==NULL) printf("I cannot dimension classes of %d elements \n",nldo*ntot);

  /* Aqui viene la normalizacion */
  for(j=0;j<ntot;j++) {
    for(i=0 ;i<nldo;i++) {
      norma+=esp[j+nldo*i];
    }
    for(ii=0 ;ii<nldo;ii++) {
      esp[j+nldo*i]=esp[j+nldo*i]/norma*nldo;
    }
  }
  /*Acaba la normalizacion */
  
  nclasses=ntot;
  for(j=0;j<ntot;j++) {
    classes[j+nldo*i]=esp[j+nldo*i];
  }
  
  for(nclasses=ntot-1;nclasses>0;nclasses--) {
    
    
    
  }
  

}
