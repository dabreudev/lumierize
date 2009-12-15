#include "modulos.h"

int fitsinfo(char *ima,FILE *stream)

{
  struct headfits h;

  if(fitsrh(ima,&h)) return(0);
  
  if(stream != NULL ) {
    fprintf(stream,"\n");
    fprintf(stream,"Descriptores de la imagen %s:\n\n",ima);
    fprintf(stream,"Dimension            %6d x %6d pix\n",h.naxis1,h.naxis2);
    fprintf(stream,"Coord. 1er pixel     %f x %f\n",h.crval1,h.crval2);
    fprintf(stream,"Tamano del pixel     %f x %f\n",h.cdelt1,h.cdelt2);
    fprintf(stream,"Unidades del pixel   %s\n",h.bunit);
    fprintf(stream,"Valores Bscale,Bzero %f,%f\n",h.bscale,h.bzero);
    fprintf(stream,"Intrumento          >%s<\n",h.instrume);
    fprintf(stream,"Objeto              >%s<\n",h.object);
    fprintf(stream,"Fecha               >%s<\n",h.date);

    if(!strncmp(h.object,"MOSAICO",7)) {
      fprintf(stream,"No. de objetos      %d x %d\n",
                       h.mosaic.nmosx,h.mosaic.nmosy);
      fprintf(stream,"Tamano de la caja   %d x %d\n",
                       h.mosaic.cajax,h.mosaic.cajay);
      }
    }
  return(1);
}
