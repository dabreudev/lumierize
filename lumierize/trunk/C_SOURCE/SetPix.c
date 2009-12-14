#include "modulos.h"


void Copy(struct image I, struct image *O, float value);

int main() {


  struct image input;
  struct image output;
  float value;

  input.aloc_flag=0; 
  output.aloc_flag=0;
  printf(" Input image: ");
  reads("",input.file);
  printf(" Output image: ");
  reads("",output.file);
  printf(" Value to set all pixels: ");
  value=readf(0.);
  ReadImage(&input);
  printf(" Image read\n");



  Copy(input, &output, value);
  printf(" Image copied\n");
  SaveImage(output);
  printf(" and saved \n");
  return 0;
}


void Copy(struct image I, struct image *O, float value) {
  int i,j;

/*   cpgopen("?"); */
/*   PlotImage(&I); */

  if(O->aloc_flag==1) free(O->array);
  O->nx=I.nx;
  O->ny=I.ny;
  O->npixels=I.npixels;
  O->naxes[0]=I.naxes[0];
  O->naxes[1]=I.naxes[1];
  O->datamin=I.datamin;
  O->datamax=I.datamax;
  O->array=vector_f(I.npixels);

/*   printf(" npixe %d nx %d ny %d\n",I.npixels,I.nx,I.ny); */

  
  for(j=0;j<I.ny;j++) {
    for(i=0;i<I.nx;i++) {
/*       printf(" %d %d \n",i+j*I.nx,(I.nx-i-1)+j*I.nx); */
      O->array[i+j*I.nx]=value;
    }
  }

/*   cpgopen("?"); */
/*   PlotImage(&(*O)); */
}
