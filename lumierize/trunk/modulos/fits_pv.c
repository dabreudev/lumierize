#include "modulos.h"
/* //#include "Util.h" */

int fits_pv(float *ima, struct headfits *h, float x, float y, float *pv)

{
  float pixx, pixy;

/* // Calculamos el punto en coordenadas pixel */
/* // ---------------------------------------- */
  pixx = (x - h->crval1) / h->cdelt1 + 1.0; 
  pixy = (y - h->crval2) / h->cdelt2 + 1.0; 

/* // Interpolamos */
/* // ------------ */
  if(!(I_Intp1p(h->naxis1, h->naxis2, ima, pixx, pixy,pv))) {
    fprintf(stderr,"fits_pv: ERROR. Coord. fuera de rango\n");
    fprintf(stderr,"         NAXIS1=%6d NAXIS2=%6d\n",h->naxis1,h->naxis2);
    fprintf(stderr,"         CRVAL1=%f CRVAL2=%f\n",h->crval1,h->crval2);
    fprintf(stderr,"         CDELT1=%f CDELT2=%f\n",h->cdelt1,h->cdelt2);
    fprintf(stderr,"         x=%f  y=%f\n",x,y);
    fprintf(stderr,"         pixx=%f  pixy=%f\n",pixx,pixy);
    return(0);
    } 

  return(1);
} 

