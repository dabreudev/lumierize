#include "modulos.h"

void f_rf(FILE *fp,struct headfits *h,float *ima)

{
  union ifc {
    int   i;
    float a;
    char  ch[4];
    } fits;

  int x,y,bpp,i;

  char bloq[2880],byte[4];

/* // Salta los bloques de la cabecera */
/* // -------------------------------- */
  rewind(fp);
  do {
    fread(bloq,2880,1,fp);
    } while (strstr(bloq,"END       ") == NULL);


/* // Lee los datos de 1, 2 o 4 bytes */
/* // ------------------------------- */

  bpp=(h->bitpix)/8;

  for(y=0;y < h->naxis2;y++)
    for(x=0;x < h->naxis1;x++) {
     fread(byte,1,abs(bpp),fp);
     fits.i=0;
     for(i=0; i<abs(bpp); i++) fits.ch[i]=byte[abs(bpp)-1-i];  /*  Vuelvo los bits */
     if((bpp==2) && (fits.i & 0x00008000)) fits.i = fits.i | 0xffff0000;
     if(bpp>0) *ima = fits.i*(h->bscale)+(h->bzero);
     else      *ima = fits.a;
     ima++;
     }

  rewind(fp);
}
