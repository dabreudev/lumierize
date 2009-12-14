#include "modulos.h"


int    Pdev(double prob) {

  float p;
  p=Constdev(0.,1.);
  if(p>prob) return(0);
  else       return(1);
}
