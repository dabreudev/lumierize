#include "modulos.h"
 
#define NBUFFER 10

/* Parametros del programa */

/* Variables comunes */ 
int nspec;
struct spectrum spec[NBUFFER];
struct spectrum errspec[NBUFFER];


void LoadParam_kbd(void);
void Plot_main(void);
void InitSpec(int ispec);


int main(void) {

  LoadParam_kbd();
  Plot_main();

  return(0);
}


void LoadParam_kbd(void) {

  int ispec;

  nspec=1;
  ispec=0;
  InitSpec(ispec);
  cpgopen("?");
}

 
void InitSpec(int ispec) {

  printf(" Input FITS file with spectrum \n");
  reads("",spec[ispec].file);
  printf(" Input FITS file with error spectrum \n");
  reads("NONE",errspec[ispec].file);
  spec[ispec].alocldo_flag=0;
  spec[ispec].aloc_flag=0;
  ReadSpec(&(spec[ispec]));
  if(strcmp(errspec[ispec].file,"NONE")) ReadSpec(&(errspec[ispec]));

}



void Plot_main(void) {
  
  int zoomflag;
  float xbmin,xbmax,ybmin,ybmax;
  char cnul;

  int ispec;
  char opt='E';
  
  do{
    
    for(ispec=0;ispec<nspec;ispec++) {
      printf(" Voy a pintar %d \n",ispec);
      if(strcmp(errspec[ispec].file,"NONE")) {
	if(zoomflag) PlotSpec_err_zoom(spec[ispec],errspec[ispec],xbmin,xbmax,ybmin,ybmax);
	else PlotSpec_err(spec[ispec],errspec[ispec]);
      }
      else    {
	if(zoomflag) PlotSpec_zoom(spec[ispec],xbmin,xbmax,ybmin,ybmax);
	else PlotSpec(spec[ispec]);
      }
    }
    printf(" R Read spectrum in buffer\n");
    printf(" L List buffers\n");
    printf(" Z Zoom\n");
    printf(" E Exit\n");
    opt=readc(opt); 
    switch (opt) { 
    case 'R':
    case 'r':
      if(nspec==NBUFFER) ispec=nspec-1;
      else ispec=nspec;
      ispec=readi(ispec+1)-1;
      break;
    case 'L':
    case 'l':
      for(ispec=0;ispec<nspec;ispec++) {
	printf(" Buffer %2d     Spectrum: %s\n",ispec+1,spec[ispec].file);
	printf("         Error Spectrum: %s\n",errspec[ispec].file);
      }
      break;
    case 'Z':
    case 'z':
      cpgsci(2);
      printf(" Press bottom left square with mouse...\n");
      cpgcurs(&xbmin,&ybmin,&cnul);
      printf(" Press top right square with mouse...\n");
      cpgband(2,1,xbmin,ybmin,&xbmax,&ybmax,&cnul);
      cpgsci(1);
      zoomflag=1;
    }
  }while(opt!='E' && opt!='e');

  cpgclos();
  
  exit(1);

}
  
