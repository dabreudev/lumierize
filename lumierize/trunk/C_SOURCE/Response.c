#include "modulos.h"
  
#define PLOT 0
#define DEBUG 0
#define DEBUG2 0
#define DEBUGINT 0
//#define MALLOC_CHECK_ 2 
 
/*Parametros del programa */

char specfile[101];
int speccol;
int errspeccol;

char outresponse[101];
char outstdresponse[101];

/* Variables globales */


char *speclist;
char *errspeclist;

int nspectra;
int nselspectra;
int npix;

struct spectrum *spectra;
struct spectrum *errspectra;
struct spectrum *selspectra;
struct spectrum *selerrspectra;


struct spectrum mediumspec;
struct spectrum stdspec;

void InitSpec(void);
void LoadParam_kbd(void);
void CreateListSpec(void);
void ReadAllSpectra(void);
void CopyAllToSelect(void);
void SelectStdSpectra(float nreject);
void SelectStdErrOffsetSpectra(float nreject);
void FinalClean(void);
void SaveResponse(void);


int main() {

  InitSpec();
  LoadParam_kbd();
  printf(" Phase 1\n");
  CreateListSpec();
  printf(" Phase 2\n");
  ReadAllSpectra();
  printf(" Phase 3\n");
  CopyAllToSelect();
  printf(" Phase 4\n");
  SelectStdSpectra(4.5);
  printf(" Phase 5\n");
  SelectStdErrOffsetSpectra(3.);
  printf(" Phase 6\n");
  FinalClean();
  printf(" Phase 7\n");
  SaveResponse();

  return(0);
}

void LoadParam_kbd(void) {

  printf(" Input file with spectra ");
  reads("",specfile);
  printf(" Input column with spectra name: ");
  speccol=readi(1);
  printf(" Input column with error spectra name: ");
  errspeccol=readi(2);
  printf(" Input file with output response ");
  reads("response.fits",outresponse);
  printf(" Input file with output error in response ");
  reads("errresponse.fits",outstdresponse);

}


void CreateListSpec(void) {

  int i;
  int nspec;
  int *islog;
  int *ieslog;

  char *slist;
  char *eslist;

  nspec=FileNLin(specfile);

  printf(" Number of spectra: %d\n",nspec);

  if(nspec<20) {
    printf(" No sense to compute response function with less than 20 spectra. Exiting\n");
    exit(1);
  }

  islog=vector_i(nspec);
  ieslog=vector_i(nspec);
  if((   slist=malloc(nspec*101*sizeof(char)))==NULL) {
    printf("I cannot dimension slist of %d bytes \n",nspec*sizeof(char));
    exit(1);
  }
  if((     eslist=malloc(nspec*101*sizeof(char)))==NULL) {
    printf("I cannot dimension eslist of %d bytes \n",nspec*sizeof(char));
    exit(1);
  }
  if((speclist=malloc(nspec*101*sizeof(char)))==NULL) {
    printf("I cannot dimension speclist of %d bytes \n",nspec*sizeof(char));
    exit(1);
  }
  if((errspeclist=malloc(nspec*101*sizeof(char)))==NULL) {
    printf("I cannot dimension errspeclist of %d bytes \n",nspec*sizeof(char));
    exit(1);
  }

  printf(" Antes leer\n"); 

  ReadCharcol(specfile,speccol   ,slist ,islog ,101,&nspec);
  printf(" Termino la primera\n");    
  ReadCharcol(specfile,errspeccol,eslist,ieslog,101,&nspec);
  nspectra=0;
  printf(" Termino ller\n"); 
 
  for(i=0;i<nspec;i++) {
/*     printf(" Voy por %d i\n",i);  */
    if(islog[i] && ieslog[i]) {
/*       printf(" Incluyo %d \n",nspectra); */
/*       printf(" Son %s %s \n",slist+i*101,eslist+i*101); */
      strcpy(speclist   +nspectra*101,slist +i*101); 
/*       printf(" La primera\n"); */
      strcpy(errspeclist+nspectra*101,eslist+i*101); 
/*       printf(" vale %s \n",speclist+nspectra*101);  */
      nspectra++;
    }
  }

  free(slist);
  free(eslist);
  free(islog);
  free(ieslog);
} 

void ReadAllSpectra(void) {

  int i;

  printf(" Entro en readsl nspectra %d\n",nspectra);
  npix=0;

  if((spectra   =malloc(nspectra*sizeof(struct spectrum)))==NULL) {
    printf("I cannot dimension spectra of %d bytes \n",nspectra*sizeof(struct spectrum));
    exit(1);
  }
  if((errspectra=malloc(nspectra*sizeof(struct spectrum)))==NULL) {
    printf("I cannot dimension errspectra of %d bytes \n",nspectra*sizeof(struct spectrum));
    exit(1);
  }


  if(PLOT) cpgopen("?");  

  for(i=0;i<nspectra;i++) {
    printf("%06d/%06d",i,nspectra);
/*     printf(" Me voy por %d \n",i); */
    strcpy(spectra[i].file   ,speclist+i*101   );
    strcpy(errspectra[i].file,errspeclist+i*101);
/*     printf(" Antes ller,\n"); */
    spectra[i].aloc_flag=0; 
    spectra[i].alocldo_flag=0;
    errspectra[i].aloc_flag=0;
    errspectra[i].alocldo_flag=0;
    ReadSpec(   spectra+i);
    ReadSpec(errspectra+i);
/*     printf(" TERmin oller\n"); */
/*     PlotSpec_pix_err(spectra[i],errspectra[i]); */
/*     i=readi(i); */
      
/*     printf(" DATAMAX %g MIN %g\n",spectra[i].datamax,spectra[i].datamin); */
    
    if(npix==0) npix=spectra[i].npixels;
    else if(npix!= spectra[i].npixels || npix!= errspectra[i].npixels)  {
      printf(" No matching pixels for all the spectra. Exiting\n");
      exit(1);
    }
  }

  free(   speclist);
  free(errspeclist);
}

void CopyAllToSelect(void) {

  int i;
  static int aloc=0;

  if(aloc) {
    free(selspectra);
    free(selerrspectra);
  }

  if((selspectra   =malloc(nspectra*sizeof(struct spectrum)))==NULL) {
    printf("I cannot dimension selspectra of %d bytes \n",nspectra*sizeof(struct spectrum));
    exit(1);
  }
  if((selerrspectra=malloc(nspectra*sizeof(struct spectrum)))==NULL) {
    printf("I cannot dimension selerrspectra of %d bytes \n",nspectra*sizeof(struct spectrum));
    exit(1);
  }
  aloc=1;

  for(i=0;i<nspectra;i++) {
    selspectra[i].nx=spectra[i].nx;
    selspectra[i].npixels=spectra[i].npixels;
    selspectra[i].datamin=spectra[i].datamin;
    selspectra[i].datamax=spectra[i].datamax;
    selspectra[i].ldomin=spectra[i].ldomin;
    selspectra[i].deltaldo=spectra[i].deltaldo;
    selspectra[i].spec=vector_f(selspectra[i].nx);
    selspectra[i].ldo=vector_f(selspectra[i].nx);
    selspectra[i].aloc_flag=1;
    selspectra[i].alocldo_flag=1;
    memcpy(selspectra[i].spec,spectra[i].spec,spectra[i].nx*sizeof(float));
    memcpy(selspectra[i].ldo,spectra[i].ldo,spectra[i].nx*sizeof(float));
    memcpy(selspectra[i].file,spectra[i].file,200*sizeof(char));
    selerrspectra[i].nx=errspectra[i].nx;
    selerrspectra[i].npixels=errspectra[i].npixels;
    selerrspectra[i].datamin=errspectra[i].datamin;
    selerrspectra[i].datamax=errspectra[i].datamax;
    selerrspectra[i].ldomin=errspectra[i].ldomin;
    selerrspectra[i].deltaldo=errspectra[i].deltaldo;
    selerrspectra[i].spec=vector_f(selerrspectra[i].nx);
    selerrspectra[i].ldo=vector_f(selerrspectra[i].nx);
    selerrspectra[i].aloc_flag=1;
    selerrspectra[i].alocldo_flag=1;
    memcpy(selerrspectra[i].spec,errspectra[i].spec,errspectra[i].nx*sizeof(float));
    memcpy(selerrspectra[i].ldo,errspectra[i].ldo,errspectra[i].nx*sizeof(float));
    memcpy(selerrspectra[i].file,errspectra[i].file,200*sizeof(char));
  }
  nselspectra=nspectra;
}

void SelectStdSpectra(float nreject) {
   
  float *pixvalue;
  float *medium;
  float *std;

  float first,median,third;
  float mmin,mmax;

  int i,j;

  medium  =vector_f(npix);
  std     =vector_f(npix);
  pixvalue=vector_f(nselspectra);


  if(DEBUG) printf(" Primer caluclo de medium\n");
  if(DEBUG) printf(" Uso %d espectros\n",nselspectra);

      

  for(j=0;j<npix;j++) {
    for(i=0;i<nselspectra;i++) {
      pixvalue[i]=(selspectra[i]).spec[j]/selspectra[i].datamax;
      if(DEBUG2) printf(" pix %d esp %d  va %g    %g   %g\n",j,i,pixvalue[i],(selspectra[i]).spec[j],selspectra[i].datamax); 
    }
    Quartil(nselspectra,pixvalue,&first,&median,&third);
    medium[j]=median;
    std[j]=(third-first)/1.35;
    if(DEBUG) printf(" %d medium %f std %f  f %g m %g t %g\n",j,medium[j],std[j],first,median,third); 
  }

  /* Voy a normalizar el espectro medio. Asi me aseguro que se va a uno el maximo */
  MinMax(npix,medium,&mmin,&mmax);
  for(j=0;j<npix;j++) {
    medium[j]=medium[j]/mmax;
    std[j]=std[j]/mmax;
  }

  mediumspec.aloc_flag=0;
  mediumspec.alocldo_flag=0;
  stdspec.aloc_flag=0;
  stdspec.alocldo_flag=0;
  FillSpec(&mediumspec,npix,medium,0.,(float)npix);
  FillSpec(  & stdspec,npix,std   ,0.,(float)npix);


  if(PLOT) {
    PlotSpec_pix(mediumspec);
    PlotSpec_pix_err(mediumspec,stdspec);
    cpgopen("?");
  }

  if(PLOT) {
    cpgswin(0.,(float)npix,0.,1.5);
    cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
    for(i=0;i<nselspectra;i++) {
      cpgswin(0.,(float)npix,0.,selspectra[i].datamax*1.5);
/*       PlotSpec_pix_err(selspectra[i],selerrspectra[i]); */
/*       j=readi(j); */
    }
  }

  if(DEBUG) i=readi(i);

      

  
/*   SaveSpec(mediumspec);   */

  for(j=0;j<npix;j++) {
    if(DEBUG) printf(" Me voy por pixel %d\n",j);
    for(i=0;i<nselspectra;i++) {
      if((fabs((selspectra[i]).spec[j]/(selspectra[i]).datamax-medium[j])/std[j])>nreject ||
	 !(finite((selspectra[i]).spec[j])) ) {
 	if(DEBUG) printf(" Deleting %d \n",i); 
   	if(PLOT) PlotSpec(selspectra[i]);   
 	if(DEBUG) printf(" es %g  me %f std %f  re %f fi %d\n",(selspectra[i]).spec[j]/(selspectra[i]).datamax,medium[j],std[j],(fabs((selspectra[i]).spec[j]/(selspectra[i]).datamax-medium[j])/std[j]),(finite((selspectra[i]).spec[j]))); 
 	if(DEBUG2) j=readi(j);  
 	
/* 	printf(" Antes ºnselspecrtra %d\n",nselspectra); */
 	memmove(  selspectra+i ,  selspectra+i+1   ,((nselspectra)-i-1)*sizeof(struct spectrum));
/* 	printf("Despues 1st move\n"); */
	memmove(selerrspectra+i,selerrspectra+i+1  ,((nselspectra)-i-1)*sizeof(struct spectrum));
/* 	printf("Despues move\n"); */
	nselspectra--;
	i--;
      }
    }
  }

  if(DEBUG) printf(" Segundo caluclo de medium\n");

  if(DEBUG) printf(" Uso %d espectros\n",nselspectra);

  for(j=0;j<npix;j++) {
    for(i=0;i<nselspectra;i++) {
      pixvalue[i]=(selspectra[i]).spec[j]/selspectra[i].datamax;
/*       printf(" pix va %g    %g   %g\n",pixvalue[i],(selspectra[i]).spec[j],selspectra[i].datamax); */
    }
    Quartil(nselspectra,pixvalue,&first,&median,&third);
    medium[j]=median;
    std[j]=(third-first)/1.35;
    if(DEBUG) printf(" %d medium %f std %f  f %g m %g t %g\n",j,medium[j],std[j],first,median,third); 
  }
  /* Voy a normalizar el espectro medio. Asi me aseguro que se va a uno el maximo */
  MinMax(npix,medium,&mmin,&mmax);
  for(j=0;j<npix;j++) {
    medium[j]=medium[j]/mmax;
    std[j]=std[j]/mmax;
  }
  FillSpec(&mediumspec,npix,medium,0.,(float)npix);
  FillSpec(  & stdspec,npix,std   ,0.,(float)npix);

  if(PLOT) {
    PlotSpec_pix(mediumspec);
    PlotSpec_pix_err(mediumspec,stdspec);
    cpgopen("?");
  }

  if(DEBUGINT) j=readi(j);

  free(pixvalue);
  free(std);
}


void SelectStdErrOffsetSpectra(float nreject) {

  float *newmedium;
  float *newstd;
  float shift;

  float *shifted;
  float *pixvalue;
  float *errpixvalue;
  float *medium;
  float *std;
  float *meanerr;
  float *stderr;
  float *meanlogerr;
  float *stdlogerr;

  float first,median,third;
  float mmin,mmax;

  int i,j;

  newmedium=vector_f(npix);
  newstd=vector_f(npix);
  shifted=vector_f(npix);
  meanerr=vector_f(npix);
  stderr=vector_f(npix);

  CopyAllToSelect();

  /* Aqui hago un shift de todos los espectros */


  for(i=0;i<nselspectra;i++) {
    printf(" Computing shift for spectra %d\n",i);
    shift=offsetpix(npix,mediumspec.spec,selspectra[i].spec);
/*     if(PLOT) PlotSpec_pix(mediumspec);  */
/*     if(PLOT) PlotSpec_pix_ov(selspectra[i]);  */
    /*     cpgswin(selspectra[i].ldomin,selspectra[i].ldomin+selspectra[i].nx*selspectra[i].deltaldo,selspectra[i].datamin-(selspectra[i].datamax-selspectra[i].datamin)*.2,selspectra[i].datamax+(selspectra[i].datamax-selspectra[i].datamin)*.45); */
    if(fabs(shift)>10) shift=0;
    shiftspec(shift,npix,selspectra[i].spec,shifted);
    memcpy(selspectra[i].spec,shifted,npix*sizeof(float));
    /*     PlotSpec_pix_ov(selspectra[i]); */
    if(DEBUG) printf(" La %d tiene offset %f \n",i,shift);  
    /*     i=readi(i);  */    
  }
  
/*   printf(" Salio de shiftspec\n"); */
  
/*   printf(" npix %d\n",npix);   */
  /* Ahora mismo estan todos los espectros shifted. Puedo volver a calcular
     el espectro medio */
  if(DEBUG) printf(" allocateo con nselspectra %d sizeof %d\n",nselspectra,sizeof(float));
  medium  =vector_f(npix);
  std     =vector_f(npix);
  meanlogerr=vector_f(npix);
  stdlogerr=vector_f(npix);
  pixvalue=vector_f(nselspectra);
  errpixvalue=vector_f(nselspectra);

  if(DEBUG) printf(" Todos los punteros\n");
/*   if(DEBUG) printf(" Punteros med %d pix %d errpix %d std %d meer %d stder %d menlo %d stdlo %d newme %d new s %d sh %d\n", */
/* 	 medium,pixvalue,errpixvalue,std,meanerr,stderr,meanlogerr,stdlogerr,newmedium,newstd,shifted); */

/*   printf(" MAS antes free errpixvalue\n"); */
/*   free(errpixvalue); */
/*   printf(" 3\n"); */


  for(j=0;j<npix;j++) {
    for(i=0;i<nselspectra;i++)   pixvalue[i]=(selspectra[i]).spec[j]/selspectra[i].datamax; 
    Quartil(nselspectra,pixvalue,&first,&median,&third);
    medium[j]=median;
    std[j]=(third-first)/1.35;
    for(i=0;i<nselspectra;i++)   pixvalue[i]=(selerrspectra[i]).spec[j]/selspectra[i].datamax;   
    Quartil(nselspectra,pixvalue,&first,&median,&third);
    meanerr[j]=median;
    stderr[j]=(third-first)/1.35;
    for(i=0;i<nselspectra;i++)   pixvalue[i]=log10(fabs((selerrspectra[i]).spec[j]/selspectra[i].datamax));  
    Quartil(nselspectra,pixvalue,&first,&median,&third); 
    meanlogerr[j]=median; 
    stdlogerr[j]=(third-first)/1.35;    
    printf(" Calculo medias j %d med %f meanerr %f meanlog %f\n",j,medium[j],meanerr[j],meanlogerr[j]);
  }



  if(DEBUG) {
    printf(" MEDIUM primero \n");
    for(j=0;j<npix;j++) {
      printf(" pix %d med %f sig %f\n",j,medium[j],stderr[j]);
    }
  }
  /* Voy a normalizar el espectro medio. Asi me aseguro que se va a uno el maximo */
  MinMax(npix,medium,&mmin,&mmax);
  for(j=0;j<npix;j++) {
    medium[j]=medium[j]/mmax;
    std[j]=std[j]/mmax;
  }

  FillSpec(&mediumspec,npix,medium,0.,(float)npix);
  FillSpec(  & stdspec,npix,std   ,0.,(float)npix);
  if(PLOT) PlotSpec_pix(mediumspec); 
  if(PLOT) PlotSpec_pix_err(mediumspec,stdspec); 
  if(DEBUG) printf(" He Calculado  la primera media dentro de Err\n");
  if(DEBUG) printf(" Esta media es con %d\n",nselspectra);
  if(DEBUGINT) j=readi(j);
        
  /* Elimino los que se van, igual que antes */


  
  for(j=0;j<npix;j++) {
    if(DEBUG) printf(" Para %d : mean %f std %f meanerr %f\n",j,medium[j],std[j],meanerr[j]);
    for(i=0;i<nselspectra;i++) {
      if((fabs((selspectra[i]).spec[j]/(selspectra[i]).datamax-medium[j])/std[j])>5.5 || 
	 fabs(((selspectra[i]).spec[j]/(selspectra[i]).datamax-medium[j])/selerrspectra[i].spec[j]*selspectra[i].datamax*meanerr[j]/std[j])>5 ||
	 fabs((log10((selerrspectra[i]).spec[j]/selspectra[i].datamax)-meanlogerr[j])/stdlogerr[j])>3 ||
	 !(finite((selspectra[i]).spec[j])) ||
	 !(finite((selerrspectra[i]).spec[j])) ) {
 	if(DEBUG) printf(" Deleting %d \n",i); 
   	if(PLOT) PlotSpec_pix_err(selspectra[i],selerrspectra[i]);   
 	if(DEBUG) printf(" es %g  me %f std %f  re %f (5.5) fi %d fi %d\n",(selspectra[i]).spec[j]/(selspectra[i]).datamax,medium[j],std[j],(fabs((selspectra[i]).spec[j]/(selspectra[i]).datamax-medium[j])/std[j]),(finite((selspectra[i]).spec[j])),(finite((selerrspectra[i]).spec[j]))); 
 	if(DEBUG) printf(" es %g  err %f me %f  std/meanerr %f re %f (5)\n",(selspectra[i]).spec[j]/(selspectra[i]).datamax,selerrspectra[i].spec[j]/selspectra[i].datamax,medium[j],std[j]/meanerr[j],fabs(((selspectra[i]).spec[j]/(selspectra[i]).datamax-medium[j])/selerrspectra[i].spec[j]*selspectra[i].datamax*meanerr[j]/std[j]));
 	if(DEBUG) printf(" log10erres %g  me %f std %f  logre %f (3)\n",log10((selerrspectra[i]).spec[j]/(selspectra[i]).datamax),meanlogerr[j],stdlogerr[j],	 fabs((log10((selerrspectra[i]).spec[j]/selspectra[i].datamax)-meanlogerr[j])/stdlogerr[j]));
 	if(DEBUG2) j=readi(j);  
 	memmove(  selspectra+i ,  selspectra+i+1   ,((nselspectra)-i-1)*sizeof(struct spectrum));
	memmove(selerrspectra+i,selerrspectra+i+1  ,((nselspectra)-i-1)*sizeof(struct spectrum));
	nselspectra--;
	i--;
      }
    }
  }

  /* Ahora calculo la media promediada */




 
  for(j=0;j<npix;j++) {
    for(i=0;i<nselspectra;i++) {
      pixvalue[i]=(selspectra[i]).spec[j]/selspectra[i].datamax;
      errpixvalue[i]=(selerrspectra[i]).spec[j]/selspectra[i].datamax/meanerr[j]*std[j];
      if(DEBUG) printf(" %d pixval %f  errpixval %f  \n",i,pixvalue[i],errpixvalue[i]);
    }
    medium[j]=StMedia(nselspectra,pixvalue,std+j);
    if(DEBUG) printf(" Con normal %f std %f\n",medium[j],std[j]);
    Quartil(nselspectra,pixvalue,&first,&median,&third);
    medium[j]=StErrWeightMedia(nselspectra,pixvalue,errpixvalue,std+j);
    StMedia(nselspectra,pixvalue,std+j);
    /* utilizo el promedio pesado para la media, pero la std es la normal, sin pesos ni nada */
    if(DEBUG) printf(" pixel %d Con err %f std %f\n",j,medium[j],std[j]);
  }

  if(DEBUG) {
    printf(" MEDIUM segundo \n");
    for(j=0;j<npix;j++) {
      printf(" pix %d med %f sig %f\n",j,medium[j],stderr[j]);
    }
  }
  /* En esta ya no normalizo el espectro medio. Me fio ya de casi todos los puntos */


  FillSpec(&mediumspec,npix,medium,0.,(float)npix);
  FillSpec(  & stdspec,npix,std   ,0.,(float)npix);
  if(PLOT) PlotSpec_pix_err(mediumspec,stdspec);  
  if(DEBUG) printf(" Calculo la penultima media promerdiada\n");
  if(DEBUG) printf(" Esta media es con %d\n",nselspectra);
  if(DEBUG) j=readi(j);
 

  /* Elimino spectros que se van  por encima */
  for(j=0;j<npix;j++) {
    for(i=0;i<nselspectra;i++) {
      if((fabs((selspectra[i]).spec[j]/(selspectra[i]).datamax-medium[j])/std[j])>nreject 
	 /* 	 || fabs(((selspectra[i]).spec[j]/(selspectra[i]).datamax-medium[j])/selerrspectra[i].spec[j]*selspectra[i].datamax*meanerr[j]/std[j])>nreject  */
	 ) {
 	if(DEBUG) printf(" Deleting %d  nselspectra %d\n",i,nselspectra); 
  	if(PLOT) PlotSpec_pix_err(selspectra[i],selerrspectra[i]);  
 	if(DEBUG) printf(" es %g  me %f std %f  re %f segerr %f re2 %f\n",(selspectra[i]).spec[j]/(selspectra[i]).datamax,medium[j],std[j],(fabs((selspectra[i]).spec[j]/(selspectra[i]).datamax-medium[j])/std[j]),1./(1./selerrspectra[i].spec[j]*selspectra[i].datamax*meanerr[j]/std[j]),fabs(((selspectra[i]).spec[j]/(selspectra[i]).datamax-medium[j])/selerrspectra[i].spec[j]*selspectra[i].datamax*meanerr[j]/std[j])); 
	if(DEBUG) j=readi(j);
/* 	printf(" Antes ºnselspecrtra %d\n",nselspectra); */
 	memmove(  selspectra+i ,  selspectra+i+1   ,((nselspectra)-i-1)*sizeof(struct spectrum));
/* 	printf("Despues 1st move\n"); */
	memmove(selerrspectra+i,selerrspectra+i+1  ,((nselspectra)-i-1)*sizeof(struct spectrum));
/* 	printf("Despues move\n"); */
	nselspectra--;
	i--;
      }
    }
  }



  /* Vuelvo a calcular la media promediada */

 
  for(j=0;j<npix;j++) {
    for(i=0;i<nselspectra;i++) {
      pixvalue[i]=(selspectra[i]).spec[j]/selspectra[i].datamax;
      errpixvalue[i]=(selerrspectra[i]).spec[j]/selspectra[i].datamax/meanerr[j]*std[j];
    }
    medium[j]=StErrWeightMedia(nselspectra,pixvalue,errpixvalue,std+j);
    std[j]=std[j]/sqrt(nselspectra); /* Para hallar el error en medium, ya no es la desviacion standard */
    Quartil(nselspectra,pixvalue,&first,&median,&third);
/*     printf(" %d medium %f std %f  f %g m %g t %g std_q %g\n",j,medium[j],std[j],first,median,third,(third-first)/1.35);  */
  }

  if(DEBUG) {
    printf(" MEDIUM tercero y final \n");
    for(j=0;j<npix;j++) {
      printf(" pix %d med %f sig %f\n",j,medium[j],stderr[j]);
    }
  }




  FillSpec(&mediumspec,npix,medium,0.,(float)npix);
  FillSpec(  & stdspec,npix,std   ,0.,(float)npix);
  if(PLOT) PlotSpec_pix_err(mediumspec,stdspec);  

  if(DEBUG) printf(" Ultima vez la media promerdiada\n");
  if(DEBUG) printf(" Esta media es con %d\n",nselspectra);
  if(DEBUGINT) j=readi(j);

  if(DEBUG) printf(" Este ya es el definitivo\n");



  

  if(DEBUGINT) printf(" Pulsa una tecla para ver todos los espectros ");

  if(DEBUGINT) {
    for(i=0;i<nselspectra;i++) {
      cpgswin(0.,(float)npix,(selspectra[i]).datamin-((selspectra[i]).datamax-(selspectra[i]).datamin)*.2,(selspectra[i]).datamax+((selspectra[i]).datamax-(selspectra[i]).datamin)*.45);
      PlotSpec_pix_ov(selspectra[i]);
      printf(" %s \n",selspectra[i].file);
    }
  }


  if(DEBUG) printf(" Antes free\n");
/*   if(DEBUG) printf(" Punteros med %d pix %d errpix %d std %d meer %d stder %d menlo %d stdlo %d newme %d new s %d sh %d\n", */
/* 	 medium,pixvalue,errpixvalue,std,meanerr,stderr,meanlogerr,stdlogerr,newmedium,newstd,shifted); */
  free(medium);
  if(DEBUG) printf(" 1\n");
  free(pixvalue);
  if(DEBUG) printf(" 2\n");
  if(DEBUG) printf(" antes free errpixvalue\n");
  free(errpixvalue);
  if(DEBUG) printf(" 3\n");
  free(std);
  if(DEBUG) printf(" 4\n");
  free(meanerr);
  if(DEBUG) printf(" 5\n");
  free(stderr);
  if(DEBUG) printf(" 6\n");
  free(meanlogerr);
  if(DEBUG)  printf(" 7\n");
  free(stdlogerr);
  if(DEBUG) printf(" 8\n");

  free(newmedium);
  printf(" 9\n");
  free(newstd);
  printf(" 10\n");
  free(shifted);

  printf(" DEepus free\n");

  printf(" Number of spectra used at last: %d\n",nselspectra);
  
}


void InitSpec(void) {
  mediumspec.aloc_flag=0;
  mediumspec.alocldo_flag=0;
  stdspec.aloc_flag=0;
  stdspec.alocldo_flag=0;
  strcpy(mediumspec.file,outresponse);
  strcpy(stdspec.file,outstdresponse);
}


void FinalClean(void) {

  int i;
  
  for(i=0;i<npix;i++) {
    if(mediumspec.spec[i]<0.03) {
      mediumspec.spec[i]=0.;
/*       stdspec.spec[i]=0.; */
    }
  }
  
}

 
void SaveResponse(void) {

  strcpy(mediumspec.file,outresponse);
  strcpy(stdspec.file,outstdresponse);
  SaveSpec(mediumspec);
  SaveSpec(stdspec);
}
