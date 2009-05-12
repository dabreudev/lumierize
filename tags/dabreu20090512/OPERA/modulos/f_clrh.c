#include "modulos.h"

void f_clrh(struct headfits *h)

{  

  h->simple  ='F';
  h->bitpix  =0;
  h->naxis   =0;
  h->naxis1  =0;
  h->naxis2  =0;
  h->crval1  =0.0;
  h->crval2  =0.0;
  h->bscale  =0.0;
  h->bzero   =0.0;
  h->datamax =0.0;
  h->datamin =0.0;
  h->cdelt1  =0.0;
  h->cdelt2  =0.0;
  sprintf(h->date,"\?\?\/\?\?\/\?\?");
  sprintf(h->instrume,"UNKNOWN");
  sprintf(h->bunit,"UNKNOWN");
  sprintf(h->object,"UNKNOWN");
  sprintf(h->origin,"UNKNOWN");
  sprintf(h->ctype1,"UNKNOWN");
  sprintf(h->ctype2,"UNKNOWN");
  sprintf(h->filename,"UNKNOWN");
  h->mosaic.nmosx=0;
  h->mosaic.nmosy=0;
  h->mosaic.cajax=0;
  h->mosaic.cajay=0;
  f_clrl(h->history);
  f_clrl(h->comment);
}
