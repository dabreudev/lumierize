#include "modulos.h"


float Band2ZP(char photband[51]);


int  main()
{
  
  
  char skyfile[51];
  char filterfile[51];
  char outfile[51];
  char photband[51];
  float zeropoint;
  FILE *fs,*ff,*fo;
  float *sky,*ldos;
  float *filter,*ldof;
  int nsky;
  int nfilter;
  int i;

  float step;
  float ldo;

  float ldoefecfil;
  float deltafilt;
  float filtint=0;
  float filtldoint=0;
  float filtldoefec;


  int  ninteg;
  
  float ldostart=4250;
  float ldoend=7420;
  float skyflux;
  float skyflux_A;
  float magsky;

  float finalmag;

  printf(" Input sky spectrum in W/m2/Ang/arcsec2: ");
  getline(skyfile,100);

  printf(" Input filter response (1 normalized) : ");
  getline(filterfile,100);

  printf(" Input which band does it correspond: ");
  getline(photband,100);


  zeropoint=Band2ZP(photband);
  
/*   strcpy(skyfile,"cielo.leinert.dat"); */
  fs=fopen(skyfile,"ro");
  ff=fopen(filterfile,"ro");
  
  nsky=FileNLin(skyfile);
  nfilter=FileNLin(filterfile);
  
  
  if((sky=malloc(nsky*sizeof(float)))==NULL) exit(1);
  if((ldos=malloc(nsky*sizeof(float)))==NULL) exit(1);

  if((filter=malloc(nfilter*sizeof(float)))==NULL) exit(1);
  if((ldof=malloc(nfilter*sizeof(float)))==NULL) exit(1);
  
  for(i=0;i<nsky;i++) {
    fscanf(fs,"%f %f\n",ldos+i,sky+i);
  }

 
  for(i=0;i<nfilter;i++) {
    fscanf(ff,"%f %f\n",ldof+i,filter+i);
  }
  
  
  MinMax(nfilter,ldof,&ldostart,&ldoend);
  step=(ldoend-ldostart)/nfilter;


  if((ldoend-ldostart)/nsky<step) step=(ldoend-ldostart)/nsky;

  ninteg=0.5*(ldoend-ldostart)/step;

  step=(ldoend-ldostart)/ninteg;

  printf(" lstart %f lend %f step %f nint %d\n",ldostart,ldoend,step,ninteg);

  // UNA vez leido el fichero, lo que hago es sumar para cada final de filtro
  
  //Asigno el principio del filtro a ldostart
  skyflux=0;
  for(i=0;i<ninteg;i++) {
    ldo=ldostart+i*step;
    skyflux+=Lagr2(ldos,sky,nsky,ldo)*Lagr2(ldof,filter,nfilter,ldo);
/*     printf(" %i %f %g %g     %g %g\n",i,ldo,skyflux,Lagr2(ldos,sky,nsky,ldo)*Lagr2(ldof,filter,nfilter,ldo),Lagr2(ldos,sky,nsky,ldo),Lagr2(ldof,filter,nfilter,ldo));   */
  }

  for(i=0;i<nfilter-1;i++) {
    filtint+=filter[i]*(ldof[i+1]-ldof[i]);
    filtldoint+=filter[i]*ldof[i]*(ldof[i+1]-ldof[i]);
  }
  ldoefecfil=filtldoint/filtint;
  filtldoefec=0;
  for(i=0;i<nfilter;i++) {
    if(ldof[i]>ldoefecfil) {
      filtldoefec=filter[i];
      break;
    }
  }

  deltafilt=filtint/filtldoefec;
  deltafilt=filtint;
  

  printf(" Filter effective wavelenght: %f\n",ldoefecfil);
  printf(" Filter response at effective wavelenght: %f\n",filtldoefec);
  printf(" Filter effective passband: %f\n",deltafilt);


  skyflux=skyflux*(ldoend-ldostart)/ninteg;
  skyflux_A=skyflux/deltafilt;
  printf(" Integrated flux in filter Flux %g\n",skyflux);
  printf(" Mean flux per wavelnght unit %g\n",skyflux_A);  
  
  
/*   cpgopen("?"); */
/*   cpgswin(ldos[istart+1],ldoend,0,skyflux[iend]*1.2); */
/*   cpgbox("BTNS",0.0,0,"BTNS",0.0,0); */
/*   cpgmove(0.,0.); */
/*   cpgsci(3); */
/*   cpgline((iend-istart),ldos+istart+1,skyflux+istart+1); */
/*   cpgsci(1); */
/*   cpgswin(ldos[istart+1],ldoend,1.e-22,6.e-20); */
//  cpgbox("BTNS",0.0,0,"CTMS",0.0,0);
/*   cpgmove(0.,0.); */
/*   cpgline((iend-istart),ldos+istart+1,sky+istart+1); */
/*   zstart=(ldos[istart+1]/6563.3)-1; */
/*   zend=(ldoend/6563.3)-1; */
/*   cpgswin(zstart,zend,0.,1.); */
/*   cpgbox("CTMS",0.0,0,"",0.0,0); */
//  cpgswin(ldos[istart+1],ldoend,0,skyflux[iend]*1.2); //Para el que viene
/*   cpgswin(ldos[istart+1],ldoend,2.0,-1.1); //Para el que viene */

  // Ahora hago lo mismo pero para el filtro desde 7220, pa tras.
/*   ldoend=7220; */
/*   for(i=istart+1;i<nsky;i++) { */
/*     if(ldos[i]>ldoend) { */
//      iend=i;
/*       break; */
/*     } */
/*   } */
  //Calculo magnilimite
/*   for(i=istart+1;i<iend+1;i++) { */
/*     ntot=i-istart; */
/*     maglim[i]=-2.5*log10(sqrt(skyflux[i]/4.7e-18));  //el 4.7 es el flujo con filtro R */
/*   } */


  magsky=-2.5*log10(skyflux_A)-zeropoint;
  printf("Magnitude in that band: %f (using zeropoint %f)\n",magsky,zeropoint);
  printf("If input spectra were erg/cm2: m = %f\n",magsky+7.5);

  printf(" Input magnitude to which you want to scale this spectrum (0=exit): ");
  finalmag=readf(0);
  if(finalmag!=0) {
    printf(" Input output file: ");
    getline(outfile,100);
/*     printf(" Paa qui\n"); */
    if((fo=fopen(outfile,"w"))==NULL) {
      printf(" Could not open %s\n",outfile);
      exit(1);
    }
/*     printf(" AIUI\n"); */
    for(i=0;i<nsky;i++) {
      fprintf(fo," %f %g\n",ldos[i],sky[i]*pow(10.,-0.4*(finalmag-magsky)));
    }
  }
  fclose(fo);

  return 0;
}




	
float Band2ZP(char photband[51]) {
  float zeropoint;
  if(!strcmp(photband,"U")) zeropoint=28.36;
  else if(!strcmp(photband,"B")) zeropoint=27.97;
  else if(!strcmp(photband,"V")) zeropoint=28.52;
  else if(!strcmp(photband,"R")) zeropoint=29.39;
  else if(!strcmp(photband,"I")) zeropoint=29.66;
  else if(!strcmp(photband,"J")) zeropoint=31.25;
  else if(!strcmp(photband,"H")) zeropoint=32.35;
  else if(!strcmp(photband,"K")) zeropoint=33.50;
  else {
    printf(" ERROR: Photometric band %-s not defined\n",photband);
    printf(" Photometric bands defined are: U B V R I J H K. \n EXITING\n");
    exit(1);
  }
  
  return(zeropoint);
}
