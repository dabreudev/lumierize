#include "modulos.h"
/* //#include "Util.h" */

int fitswf(char *name,struct headfits *h,float *ima,char autoescala)

{
    union ifc {
    int   i;
    float a;
    char  ch[4];
    } fits;

  int i,j,bpp;
  long int bb1,bb2,min=0,max=0,f;
  float amin,amax;
  FILE *fp;
  char byte[4];

/* // Test varios */
/* // ----------- */
  if(!f_testh(h)) return(0);

  if( (fp = fopen(name,"wb")) == NULL ) {
    fprintf(stderr,"fitswfd: WARNING. No puedo abrir el fichero %s\n",name);
    return(0);
    }

/* // Actualizo los descriptores datamin,datamax */
/* //                            bscale, bzero (si se pide) */
/* //                            filename */
/* // ---------------------------------------------------------- */
  MinMax(h->naxis1*h->naxis2,ima,&amin,&amax);
  h->datamin=amin;
  h->datamax=amax;
  printf("fitswf: MESSAGE.\n");
  printf("                 DATAMAX= %f\n",h->datamax);
  printf("                 DATAMIN= %f\n",h->datamin);

  strcpy(h->filename,name);

  if((autoescala == 'T') && (h->bitpix > 0)) {
    if(h->bitpix == 8) {
      h->bzero = h->datamin;
      h->bscale = (h->datamax - h->datamin)/255;
      }
    else {
      bb1=powf(2.0,(float)(h->bitpix-1));
      bb2=powf(2.0,(float)h->bitpix)-1;
      h->bzero = (bb1*(h->datamax+h->datamin)-h->datamin)/bb2;
      h->bscale = (h->datamax-h->datamin)/bb2;
      }
    printf("fitswf: MESSAGE. He cambiado los valores BSCALE y BZERO:\n");
    printf("                 BSCALE = %f\n",h->bscale);
    printf("                 BZERO  = %f\n",h->bzero);
    }


/* // Escribo la cabecera */
/* // ------------------- */
  f_wh(fp,h);

  bpp=h->bitpix/8;
  if (bpp > 0) {
    if(bpp==1) {
      min=0;
      max=255;
    }
    else {
      min=-powf(2.0,(float)(h->bitpix-1));
      max=-min-1;
    }
}

  /*  Bucle de escritura de los datos */
  /*  ------------------------------- */

  for(i=0; i < h->naxis1 * h->naxis2; i++) {
    if(bpp < 0) fits.a=*ima;
    else {
/* //Esto lo he cambiado para que funcione al compilar en Linux: */
/* //      f=nintf((*ima-h->bzero)/h->bscale); */
      f=rintf((*ima-h->bzero)/h->bscale);
      if(f < min) f=min;
      if(f > max) f=max;
      fits.i=f;
      }
    for(j=0; j<abs(bpp); j++) byte[j]=fits.ch[abs(bpp)-1-j];/* // Vuelvo los bits */
    fwrite(byte,1,abs(bpp),fp);
    ima++;
    }

/* // Relleno el utltimo bloque */
/* // ------------------------- */
  f_fill(fp);
  
  fclose(fp);
  return(1);
}
