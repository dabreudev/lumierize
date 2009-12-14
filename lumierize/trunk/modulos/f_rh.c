#include "modulos.h"

int f_rh(FILE *fp,struct headfits *h)

{

  f_clrh(h);

  f_kfc(fp,"SIMPLE  =",&(h->simple_c));
  f_kfs(fp,"DATE    =",h->date);
  f_kfi(fp,"BITPIX  =",&(h->bitpix));
  f_kfi(fp,"NAXIS   =",&(h->naxis));
  f_kfi(fp,"NAXIS1  =",&(h->naxis1));
  f_kfi(fp,"NAXIS2  =",&(h->naxis2));
  f_kff(fp,"CRVAL1  =",&(h->crval1));
  f_kff(fp,"CRVAL2  =",&(h->crval2));
  if(h->bitpix>0) {
    f_kff(fp,"BSCALE  =",&(h->bscale));
    f_kff(fp,"BZERO   =",&(h->bzero));
    }
  f_kff(fp,"DATAMAX =",&(h->datamax));
  f_kff(fp,"DATAMIN =",&(h->datamin));
  f_kff(fp,"CDELT1  =",&(h->cdelt1));
  f_kff(fp,"CDELT2  =",&(h->cdelt2));
  f_kfs(fp,"INSTRUME=",h->instrume);
  f_kfs(fp,"OBJECT  =",h->object);
  f_kfs(fp,"BUNIT   =",h->bunit);
  f_kfs(fp,"ORIGIN  =",h->origin);
  f_kfs(fp,"CTYPE1  =",h->ctype1);
  f_kfs(fp,"CTYPE2  =",h->ctype2);
  f_kfs(fp,"FILENAME=",h->filename);
  if(!strncmp(h->object,"MOSAICO",7)) {
    f_kfi(fp,"NMOSX   =",&(h->mosaic.nmosx));
    f_kfi(fp,"NMOSY   =",&(h->mosaic.nmosy));
    f_kfi(fp,"CAJAX   =",&(h->mosaic.cajax));
    f_kfi(fp,"CAJAY   =",&(h->mosaic.cajay));
    }
  f_kfl(fp,"COMMENT ",h->comment);
  f_kfl(fp,"HISTORY ",h->history);

  rewind(fp);
  return(f_testh(h));
}
