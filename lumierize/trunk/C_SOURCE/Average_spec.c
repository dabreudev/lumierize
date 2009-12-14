#include "modulos.h"



int main() {
  
  int i;
  float *spec;
  float *errspec;
  float xdata[2];
  float errxdata[2];

  char fileout[101],errfileout[101];
  
  struct spectrum spec1,spec2,specout;
  struct spectrum errspec1,errspec2,errspecout;
  
  printf(" Input name for first spectrum ");
  reads("",spec1.file);
  printf(" Input name for first error spectrum ");
  reads("",errspec1.file);
  printf(" Input name for second spectrum ");
  reads("",spec2.file);
  printf(" Input name for second error spectrum ");
  reads("",errspec2.file);
  printf(" Input name for output spectrum ");
  reads("",fileout);
  printf(" Input name for output error spectrum ");
  reads("",errfileout);
  spec1.alocldo_flag=0;
  spec1.aloc_flag=0;
  spec2.alocldo_flag=0;
  spec2.aloc_flag=0;
  specout.alocldo_flag=0;
  specout.aloc_flag=0;
  errspec1.alocldo_flag=0;
  errspec1.aloc_flag=0;
  errspec2.alocldo_flag=0;
  errspec2.aloc_flag=0;
  errspecout.alocldo_flag=0;
  errspecout.aloc_flag=0;
    


  ReadSpec(&spec1);
  ReadSpec(&spec2);
  ReadSpec(&errspec1);
  ReadSpec(&errspec2);


  if(spec1.nx==spec2.nx && spec1.ldomin==spec2.ldomin && spec1.deltaldo==spec2.deltaldo)
    ;
  else {
    printf(" nx1= %d  nx2 = %d   ldomin1 = %f ldomin2 = %f   deltaldo1 = %f deltaldo2 = %f\n",spec1.nx,spec2.nx,spec1.ldomin,spec2.ldomin,spec1.deltaldo,spec2.deltaldo);
    printf(" Spectra are not equal in spectral range. Exiting\n");
    
    exit(1);
  }
  spec=vector_f(spec1.nx);
  errspec=vector_f(spec1.nx);
  
  for(i=0;i<spec1.nx;i++) {
    xdata[0]=spec1.spec[i];
    xdata[1]=spec2.spec[i];
    errxdata[0]=errspec1.spec[i];
    errxdata[1]=errspec2.spec[i];
    if(errxdata[0]==0 || errxdata[1]==0)    spec[i]=StMedia(2,xdata,errspec+i);
    else     spec[i]=StErrWeightMedia(2,xdata,errxdata,errspec+i);
    
    printf(" x1 %f e1 %f x2 %f e2 %f  s %f e %f\n",spec1.spec[i],errspec1.spec[i],spec2.spec[i],errspec2.spec[i],spec[i],errspec[i]);
  }
  printf(" ANTES\n");
  
  FillSpec(&specout,spec1.nx,spec,spec1.ldomin,spec1.deltaldo*spec1.nx+spec1.ldomin);
  FillSpec(&errspecout,spec1.nx,errspec,spec1.ldomin,spec1.deltaldo*spec1.nx+spec1.ldomin);

  printf(" BIEN\n");

  strcpy(specout.file,fileout);
  strcpy(errspecout.file,errfileout);

  printf(" AUI\n");
  
  cpgopen("?");
  PlotSpec_pix(spec1);
  PlotSpec_pix_ov(spec2);
  PlotSpec_pix_ov(specout);

  SaveSpec(specout);
  SaveSpec(errspecout);

  return(0);
}
