#include "modulos.h"


void Copy(struct spectrum I, struct spectrum *O);

int main() {


  struct spectrum input;
  struct spectrum output;

  input.aloc_flag=0;  input.alocldo_flag=0;
  output.aloc_flag=0;  output.alocldo_flag=0;
  printf(" Input spectrum: ");
  reads("",input.file);
  printf(" Output spectrum: ");
  reads("",output.file);
  ReadSpec(&input);
  printf(" Spectrum read\n");



  Copy(input, &output);
  printf(" Spectrum copied\n");
  SaveSpec(output);
  printf(" and saved \n");
  return 0;
}


void Copy(struct spectrum I, struct spectrum *O) {
  int i;

/*   cpgopen("?"); */
/*   PlotSpectrum(&I); */

  if(O->aloc_flag==1) free(O->spec);
  if(O->alocldo_flag==1) free(O->ldo);
  O->aloc_flag=1;  O->alocldo_flag=1;
  O->nx=I.nx;
  O->npixels=I.npixels;
  O->naxes[0]=I.naxes[0];
  O->naxes[1]=I.naxes[1];
  O->datamin=I.datamin;
  O->datamax=I.datamax;
  O->ldomin=I.ldomin;
  O->deltaldo=I.deltaldo;
  O->spec=vector_f(I.npixels);
  O->ldo=vector_f(I.npixels);

/*   printf(" npixe %d nx %d ny %d\n",I.npixels,I.nx,I.ny); */

  
  for(i=0;i<I.nx;i++) {
/*       printf(" %d %d \n",i+j*I.nx,(I.nx-i-1)+j*I.nx); */
      O->spec[i]=I.spec[i];
  }

/*   cpgopen("?"); */
/*   PlotSpectrum(&(*O)); */
}
