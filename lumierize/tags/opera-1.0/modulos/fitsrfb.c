#include "modulos.h"

float *fitsrfb(char *name, struct headfits *h,
                      int *x1, int *x2, int *y1, int *y2)

{

  int dimx,dimy;
  float *ima;
  long blc[2],trc[2],incr[2],naxes[2];
  int nfound,status=0,anynull;
  fitsfile *fits;
  

  if( ffopen(&fits, name, READONLY, &status)) fits_report_error(stderr,status);
  if(fits_read_keys_lng(fits, "NAXIS", 1, 2, naxes, &nfound, &status)) fits_report_error(stderr,status); 



/* // Test */
/* // ---- */
  if(*x1 < 1) {
    printf("fitsrfb: WARNING. Coloco x1 a 1\n");
    printf("                         x1 = %d  a 1\n",*x1);
    *x1=1;
    }
  if(*y1 < 1) {
    printf("fitsrfb: WARNING. Coloco y1 a 1\n");
    printf("                         y1 = %d  a 1\n",*y1);
    *y1=1;
    }
  if(*x2 > h->naxis1) {
    printf("fitsrfb: WARNING. Coloco x2 a NAXIS1\n");
    printf("                         x2 = %d  a %d (NAXIS1)\n",*x2,h->naxis1);
    *x2=h->naxis1;
    }
  if(*y2 > h->naxis2) {
    printf("fitsrfb: WARNING. Coloco y2 a NAXIS2\n");
    printf("                         y2 = %d  a %d (NAXIS2)\n",*y2,h->naxis2);
    *y2=h->naxis2;
    }

  dimx=*x2-*x1+1;
  dimy=*y2-*y1+1;
  if( (dimx<1) || (dimy<1) ) {
    printf("fitsrfb: ERROR: Las dim. deben ser positivas \n");
    printf("x1=%6d, x2=%6d :  Dimx=%6d\n",*x1,*x2,dimx);
    printf("y1=%6d, y2=%6d :  Dimy=%6d\n",*y1,*y2,dimy);
    return(NULL);
    }

  if( (ima=malloc(dimx*dimy*sizeof(float))) == NULL) {
    printf("fitsrfb: ERROR. No puedo dimensionar la matriz ima ");
    printf("de %d bytes",dimx*dimy*sizeof(float));
    exit(1);
    }

  blc[0]=*x1;blc[1]=*y1;trc[0]=*x2;trc[1]=*y2;incr[0]=1;incr[1]=1;
/*  printf(" x1 %d x2 %d y1 %d y2 %d\n",blc[0],blc[1],trc[0],trc[1]);
  printf(" axes %d %d \n",naxes[0],naxes[1]); */
  if(fits_read_subset_flt(fits,1,2,naxes,blc,trc,incr,0.,ima,&anynull,&status)) fits_report_error(stderr,status);
 /* f_rfb(fp,h,ima,*x1,*x2,*y1,*y2); */ /* Esta era la forma anituga */
/*  printf(" Leida \n"); */
  fits_close_file(fits,&status);
  return(ima);
}
