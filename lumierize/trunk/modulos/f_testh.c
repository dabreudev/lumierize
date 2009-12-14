#include "modulos.h"

int f_testh(struct headfits *h)

{
  if(!h->simple) {
    fprintf(stderr,"f_rh: ERROR. La imagen FITS no es standard\n");
    fprintf(stderr,"              SIMPLE =%c no puedo tratarla.\n",h->simple);
    exit(1);
  }
  
  if(h->naxis != 2) {
    fprintf(stderr,"f_rh: ERROR. NAXIS es distinto de 2.\n");
    fprintf(stderr,"              NAXIS = %d no puedo tratarlo.\n",h->naxis);
    exit(1);
  }
  
  if(h->bitpix != 8 && h->bitpix != 16 && h->bitpix !=32 && h->bitpix != -32 && h->bitpix != -64) {
    fprintf(stderr,"f_rh: ERROR. No puedo tratar BITPIX\n");
    fprintf(stderr,"              BITPIX = %d\n",h->bitpix);
    exit(1);
    }

  if((h->bscale == 0) && (h->bitpix > 0)) {
    printf("f_rh: WARNING. Tomo valores FITS bscale=1 bzero=0\n");
    h->bscale=1.0;
    h->bzero=0;
    }

  if((h->cdelt1) == 0) {
    printf("f_rh: WARNING. CDELT1=0. Tomo valor CDELT1=1\n");
    h->cdelt1=1.0;
    }

  if((h->cdelt2) == 0) {
    printf("f_rh: WARNING. CDELT2=0. Tomo valor CDELT2=1\n");
    h->cdelt2=1.0;
    }

  return(1);
}
