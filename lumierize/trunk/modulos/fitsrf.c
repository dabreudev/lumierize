#include "modulos.h"

float *fitsrf(char *name,struct headfits *h)

{
  float *ima;
  int dat;
  FILE *fp;

  if((fp=fopen(name,"rb")) == NULL ) {
    fprintf(stderr,"fitsrf: WARNING. No puedo abrir la imagen %s\n",name);
    return(NULL);
    }

  if(!f_rh(fp,h)) {
    fclose(fp);
    return(NULL);
    }

  dat=h->naxis1*h->naxis2;
  if((ima=malloc(dat*sizeof(float))) == NULL) {
    fprintf(stderr,"fitsrf: ERROR. No dispongo de memoria ");
    fprintf(stderr,"para %d bytes\n",sizeof(float)*dat);
    exit(1);
    }

  f_rf(fp,h,ima);
  fclose(fp);
  return(ima);
}
