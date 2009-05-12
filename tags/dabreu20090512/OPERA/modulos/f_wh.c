#include "modulos.h"

int f_wh(FILE *fp, struct headfits *h)

{
  int i;

  f_kwc(fp,"SIMPLE  =",h->simple,"FITS Standard");
  f_kwi(fp,"BITPIX  =",h->bitpix,"FITS bits/pixel");
  f_kwi(fp,"NAXIS   =",h->naxis,"No. de  ejes");
  f_kwi(fp,"NAXIS1  =",h->naxis1,"No. de pixels del eje X");
  f_kwi(fp,"NAXIS2  =",h->naxis2,"No. de pixels del eje Y");
  if(h->bitpix > 0) {
    f_kwf(fp,"BSCALE  =",h->bscale,"Valor real = Valor FITS * BSCALE + BZERO");
    f_kwf(fp,"BZERO   =",h->bzero," ");
    }
  f_kws(fp,"BUNIT   =",h->bunit,"Unidades de los pixels");
  f_kws(fp,"OBJECT  =",h->object,"Objeto");
  f_kws(fp,"ORIGIN  =",h->origin,"Origen");
  f_kws(fp,"INSTRUME=",h->instrume,"Instrumento");
  f_kwf(fp,"CRVAL1  =",h->crval1,"Coordenada X del primer pixel");
  f_kwf(fp,"CRVAL2  =",h->crval2,"Coordenada Y del primer pixel");
  f_kwf(fp,"CDELT1  =",h->cdelt1,"Tamano del pixel X");
  f_kwf(fp,"CDELT2  =",h->cdelt2,"Tamano del pixel Y");
  f_kws(fp,"CTYPE1  =",h->ctype1,"Unidades del eje X");
  f_kws(fp,"CTYPE2  =",h->ctype2,"Unidades del eje Y");
  f_kws(fp,"DATE    =",h->date,"Fecha");
  f_kws(fp,"FILENAME=",h->filename,"Nombre del fichero");
  f_kwf(fp,"DATAMAX =",h->datamax,"Maximo");
  f_kwf(fp,"DATAMIN =",h->datamin,"Minimo");

  if(!strncmp(h->object,"MOSAICO",7)) {
    f_kwi(fp,"NMOSX   =",h->mosaic.nmosx,"No. de mosaicos en el eje X");
    f_kwi(fp,"NMOSY   =",h->mosaic.nmosy,"No. de mosaicos en el eje Y");
    f_kwi(fp,"CAJAX   =",h->mosaic.cajax,
                     "Tamano X de la caja del mosaico (pix)");
    f_kwi(fp,"CAJAY   =",h->mosaic.cajay,
                     "Tamano Y de la caja del mosaico (pix)");
    }

  for(i=0; i<MAX_KEYLAB; i++) 
   if(h->history[i][0] != ' ') fprintf(fp,"HISTORY   %-70s",h->history[i]);
  
  for(i=0; i<MAX_KEYLAB; i++) 
   if(h->comment[i][0] != ' ') fprintf(fp,"COMMENT   %-70s",h->comment[i]);

  fprintf(fp,"%-80s","END");
  f_fill(fp);
  return(ftell(fp)/2880);

}
