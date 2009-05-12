#include "modulos.h"

int fitsrh(char *ima, struct headfits *h)

{
  int status=0;
  char comment[100];
  fitsfile *fits;

  if( ffopen(&fits, ima, READONLY, &status)) fits_report_error(stderr,status);

  f_clrh(h);

  ffgky(fits,TINT    ,"BITPIX" ,&(h->bitpix)   ,comment,&status);status=0;
  ffgky(fits,TLOGICAL,"SIMPLE" ,&(h->simple)   ,comment,&status);status=0;
  ffgky(fits,TSTRING ,"DATE"   ,&(h->date)     ,comment,&status);status=0;
  ffgky(fits,TINT    ,"NAXIS"  ,&(h->naxis)    ,comment,&status);status=0;
  ffgky(fits,TINT    ,"NAXIS1" ,&(h->naxis1)   ,comment,&status);status=0;
  ffgky(fits,TINT    ,"NAXIS2" ,&(h->naxis2)   ,comment,&status);status=0;
  ffgky(fits,TFLOAT  ,"CRVAL1" ,&(h->crval1)   ,comment,&status);status=0;
  ffgky(fits,TFLOAT  ,"CRVAL2" ,&(h->crval2)   ,comment,&status);status=0;
  ffgky(fits,TINT    ,"NAXIS1" ,&(h->naxis1)   ,comment,&status);status=0;
  if(h->bitpix>0) {
    ffgky(fits,TFLOAT  ,"BSCALE" ,&(h->bscale)   ,comment,&status);status=0;
    ffgky(fits,TFLOAT  ,"BZERO"  ,&(h->bzero)    ,comment,&status);status=0;
  }
  ffgky(fits,TFLOAT  ,"DATAMAX" ,&(h->datamax)  ,comment,&status);status=0;
  ffgky(fits,TFLOAT  ,"DATAMIN" ,&(h->datamin)  ,comment,&status);status=0;
  ffgky(fits,TFLOAT  ,"CDELT1"  ,&(h->cdelt1)   ,comment,&status);status=0;
  ffgky(fits,TFLOAT  ,"CDELT2"  ,&(h->cdelt2)   ,comment,&status);status=0;
  ffgky(fits,TSTRING ,"INSTRUME",&(h->instrume) ,comment,&status);status=0;
  ffgky(fits,TSTRING ,"OBJECT"  ,&(h->object)   ,comment,&status);status=0;
  ffgky(fits,TSTRING ,"BUNIT"   ,&(h->bunit)    ,comment,&status);status=0;
  ffgky(fits,TSTRING ,"ORIGIN"  ,&(h->origin)   ,comment,&status);status=0;
  ffgky(fits,TSTRING ,"CTYPE1"  ,&(h->ctype1)   ,comment,&status);status=0;
  ffgky(fits,TSTRING ,"CTYPE2"  ,&(h->ctype2)   ,comment,&status);status=0;
  ffgky(fits,TSTRING ,"FILENAME",&(h->filename) ,comment,&status);status=0;
  ffgky(fits,TSTRING ,"COMMENT" ,&(h->comment)  ,comment,&status);status=0;
  ffgky(fits,TSTRING ,"HISTORY" ,&(h->history)  ,comment,&status);status=0;

  if(h->simple) h->simple_c='T';

  f_testh(h);
  
  if(fits_close_file(fits,&status)) fits_report_error(stderr,status);
  
  return(status);
}
