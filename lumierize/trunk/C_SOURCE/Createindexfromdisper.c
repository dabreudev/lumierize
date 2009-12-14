#include "modulos.h"



int main() {

  struct disper_prism DP;
  char indexfile[101];
  float alpha,focal;
  float ldomin,ldomax,ldozero,ldo;
  float nref,nzero;
  int npt=1000;
  int i;
  FILE *of;

  printf(" Introduce prism dispersion A ");
  DP.A=readf(1);
  printf(" Introduce prism dispersion B ");
  DP.B=readf(1);
  printf(" Introduce prism dispersion C ");
  DP.C=readf(1);
  printf(" Introduce prism dispersion tampix ");
  DP.tampix=readf(1);
  printf(" Introduce alpha for this prism ");
  alpha=readf(1)/180*M_PI;
  printf(" Introduce focal lenght for the telescope");
  focal=readf(1);



  printf(" Introduce name of refraction index file");
  reads("",indexfile);

  printf(" Introduce minimum wavelenght ");
  ldomin=readf(3000);
  printf(" Introduce maximum wavelenght ");
  ldomax=readf(10000);
  printf(" Introduce wavelenght with zero deviation ");
  ldozero=readf(6500);
  printf(" Introduce index at that wavelenght" );
  nzero=readf(1.2);

  if((of=fopen(indexfile,"w")) ==NULL) {
    printf(" I could not open %s \n",indexfile);
    exit(1);
  }


  for(i=0;i<npt;i++) {
    ldo=ldomin+i*(ldomax-ldomin)/(npt-1);
    nref=nzero+((pr_ldo2pix(ldo,DP)-pr_ldo2pix(ldozero,DP))*DP.tampix*1e-6/focal)/alpha ;
    fprintf(of," %10.2f %10.6f \n",ldo,nref);
  }

  fclose(of);
  return(0);

}
