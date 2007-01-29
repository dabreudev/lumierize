#include "modulos.h"

void f_rfb(FILE *fp, struct headfits *h, float *ima,
                  int x1, int x2, int y1, int y2)

{

  int nh,dimx,dimy,bpp,x,y,seekin,i;
  char byte[4];
  union  ifc {
    int   i;
    float a;
    char  ch[4];
    } fits;

  nh=f_bhn(fp);
  bpp=h->bitpix/8;
  dimx=x2-x1+1;
  dimy=y2-y1+1;


  seekin=2880*nh+((y1-1)*h->naxis1+(x1-1))*abs(bpp);
  for(y=0; y<dimy; y++) {
    fseek(fp,seekin,0);
    for(x=0;x<dimx;x++) {
      fread(byte,1,abs(bpp),fp);
      fits.i=0;
      for(i=0; i<abs(bpp); i++) fits.ch[i]=byte[abs(bpp)-1-i];  /* // Vuelvo los bits */
      if((bpp==2) && (fits.i & 0x00008000)) fits.i = fits.i | 0xffff0000;
      if(bpp>0) *ima = fits.i*(h->bscale)+(h->bzero);
      else      *ima = fits.a;
      ima++;
      }
    seekin += h->naxis1*abs(bpp);
    }
  rewind(fp);
}
