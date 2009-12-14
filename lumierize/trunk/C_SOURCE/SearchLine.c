#include "modulos.h"

#define  FTOL    0.000005 
#define  FTOL2   0.000005 
#define  FTOLERR    0.00012
#define  FTOLERR2    0.00012
#define NITER  300
#define NITER2 300
#define NITERERR 120
#define NITERERR2 150
#define NCONFL 15
#define NCONFL2 15


#define PLOT 0
#define DEBUG 0
#define DEBUG2 0
#define DEBUGINT 0
#define INTERACT 0
 
/*Parametros del programa */
char respfile[101];
char specfile[101];
int speccol;
int errspeccol;
int xcol;
int ycol;
int racol;
int deccol;
int skycol;
char outfile[101];
char selecfile[101];
float seeing,platescale;
float sigmax;
float minselew,maxselew,minselldo,maxselldo;
struct disper_prism DP;

/* Variables para la funcion respuesta */
struct spectrum response;
struct spectrum errfromresp;

/* Variables para los espectros */
char *speclist;
char *errspeclist;

int nspectra;
int npix;

struct spectrum *spectra;
struct spectrum *errspectra;
float *xspec,*yspec,*skyspec;
double *raspec,*decspec;

struct spectrum *specamo;
struct spectrum *errspecamo;

int nomatchingpixels_flag=0;
float lowpix,uppix;
float sumresp;
double ybegin,fluxorig,x0,dfdx,sigmaflux;



void ReadResp(char respfile[]);
float resp(float ldo);
void LoadParam_file(char file[100]);
void SaveParam(void);
void LoadParam_kbd(void);
void CreateListSpec(void);
void ReadAllSpectra(void);
int FitSpec(struct spectrum spec, struct spectrum errspec,float *ew, float *errew, float *lpix, float *errlpix, float *Kfit, float *errKfit);
int Minimize(struct spectrum spec, struct spectrum errspec, float *ew, float *errew, float *lpix, float *errlpix, float *Kfit, float *errKfit);
void CreateSintec(struct spectrum *sintec, float ew_pix, float lpix, float fwhm, float K);
double CalcML(struct spectrum es, struct spectrum erres, struct spectrum si);
double CalcML2(struct spectrum es, struct spectrum erres, struct spectrum si);
double Funk_int_ml(double x);
void ParseSpectra(int *selfit, float *ewfit, float *errewfit, float *ldofit, float *errldofit, float *Kfit, float *errKfit);
void PCA(void);
double amoe_fit_sl(int n, double *x, double *y, double *p);
double amoe_fit_sl2(int n, double *x, double *y, double *p);
void PrintOutFiles(int *selfit, float *ewfit, float *errewfit, float *ldofit, float *errldofit, float *Kfit, float *errKfit);

int main(int argc, char **argv) {

  int *selfit;
  float *ewfit,*errewfit,*ldofit,*errldofit,*Kfit,*errKfit;
  char touchchar[500];

  printf(" Welcome to SearchLine 3.5\n");

  srandom((unsigned int)time(NULL)/2); 
  if(argc<2) {
    LoadParam_kbd();
  }
  else LoadParam_file(argv[1]);
  if(DEBUG) {
    cpgopen("?");
  }

    
  ReadResp(respfile);
  CreateListSpec();
  ReadAllSpectra();
  selfit=vector_i(nspectra);
  ewfit=vector_f(nspectra);
  errewfit=vector_f(nspectra);
  ldofit=vector_f(nspectra);
  errldofit=vector_f(nspectra);
  Kfit=vector_f(nspectra);
  errKfit=vector_f(nspectra);

  printf("done\n");
  printf(" Searching for emission lines...\n");
  ParseSpectra(selfit,ewfit,errewfit,ldofit,errldofit,Kfit,errKfit);
  printf(" done\n");
  printf(" Writing results...\n");
  PrintOutFiles(selfit,ewfit,errewfit,ldofit,errldofit,Kfit,errKfit);
  printf(" done\n");
  sprintf(touchchar,"/bin/touch %s.sldone\n",specfile);
  system(touchchar);

  return(0);
}

void ParseSpectra(int *selfit, float *ewfit, float *errewfit, float *ldofit, float *errldofit, float *Kfit, float *errKfit) {
  int i; 
  int ipix;
  int seltot=0;
  float sigma;
  float *lpixfit,*errlpixfit;
  
  lpixfit=vector_f(nspectra);
  errlpixfit=vector_f(nspectra);

/*   PCA();  */
  lowpix=pr_ldo2pix(minselldo,DP);
  uppix=pr_ldo2pix(maxselldo,DP);

  errfromresp.aloc_flag=0;
  errfromresp.alocldo_flag=0;
  CopySpec(&errfromresp,response);
  for(ipix=0;ipix<response.nx;ipix++) {
    if(response.spec[ipix]<0.08) ;
    else {
      StMedia(3,response.spec+ipix-1,&sigma);
      errfromresp.spec[ipix]=sigma/5.;
    }
  }



  sumresp=StSuma1((int)uppix-(int)lowpix+1,response.spec+(int)lowpix,1);
  for(i=0;i<nspectra;i++) {
    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    printf("%06d/%06d   Selected: %06d",i,nspectra,seltot);
    fflush(NULL);
    selfit[i]=FitSpec(spectra[i],errspectra[i],ewfit+i,errewfit+i,lpixfit+i,errlpixfit+i,Kfit+i, errKfit+i); 
    ldofit[i]=pr_pix2ldo(lpixfit[i],DP);
    errldofit[i]=pr_dldp(lpixfit[i],DP)*errlpixfit[i];
    seltot+=selfit[i];
    if(DEBUG) i=readi(i);
  }

  free(lpixfit);
  free(errlpixfit);

}


void PrintOutFiles(int *selfit, float *ewfit, float *errewfit, float *ldofit, float *errldofit, float *Kfit, float *errKfit) {

  FILE *fo, *fc;
  int i;

  if((fo=fopen(outfile,"w"))==NULL) {
    printf("ERROR: Can't open output fitting file %s\n",outfile);
    exit(1);
  }
  if((fc=fopen(selecfile,"w"))==NULL) {
    printf("ERROR: Can't open candidate file %s\n",selecfile);
    exit(1);
  }

  fprintf(fo,"# Output Fitting file from SearchLine\n");
  fprintf(fo,"# Original spectra file: %s\n",specfile);
  fprintf(fo,"#%-101s  %-101s  %-15s %-15s %-15s\n","Spectrum","Err_spectrum","EW_fit","LDO_fit","K_lin");

  fprintf(fc,"# Output candidate file from SearchLine\n");
  fprintf(fc,"# Original spectra file: %s\n",specfile);
  fprintf(fc,"# Selection cuts: EW %f - %f   LDO %f -%f\n",minselew,maxselew,minselldo,maxselldo);
  fprintf(fc,"#%-101s  %-101s  %-15s %-15s %-15s %-15s %-15s\n","Spectrum","Err_spectrum","EW_fit","ERR_EW_fit","LDO_fit","ERR_LDO_fit","K_lin");
  
  for(i=0;i<nspectra;i++) {
    fprintf(fo,"%-101s  %-101s   %-15.5g %-15.5g %-15.10g\n",spectra[i].file,errspectra[i].file,ewfit[i],ldofit[i],Kfit[i]);
    if(selfit[i]) {
      fprintf(fc,"%-101s  %-101s   %-15.5g %-15.5g %-15.5g %-15.5g %-15.10g\n",spectra[i].file,errspectra[i].file,ewfit[i],errewfit[i],ldofit[i],errldofit[i],Kfit[i]);
    }
  }

  fclose(fo);
  fclose(fc);

}

void LoadParam_kbd(void) {


  printf(" Input file with spectra ");
  reads("",specfile);
  printf(" Input column with spectra name: ");
  speccol=readi(1);
  printf(" Input column with error spectra name: ");
  errspeccol=readi(2);
  printf(" Input column with X extraction position: ");
  xcol=readi(3);
  printf(" Input column with Y extraction position: ");
  ycol=readi(4);
  printf(" Input column with RA: ");
  racol=readi(5);
  printf(" Input column with DEC: ");
  deccol=readi(6);
  printf(" Input column with sky: ");
  skycol=readi(7);
  printf(" Ouput fitting file: ");
  reads(outfile,outfile);
  printf(" Output file with candidates: ");
  reads(selecfile,selecfile);
  printf(" FITS file with response function of the system (NONE for no response function: ");
  reads(respfile,respfile);
  printf(" Minimum EW selected: ");
  minselew=readf(minselew);
  printf(" Maximum EW selected: ");
  maxselew=readf(maxselew);
  printf(" Seeing: ");
  seeing=readf(seeing);
  printf(" Plate scale (\"/pix): ");
  platescale=readf(platescale);
  printf(" Input error in X spectra extraction: ");
  sigmax=readf(1.);
  printf(" Prism dispersion Constant A: ");
  DP.A=readf(1.);
  printf(" Prism dispersion Constant B: ");
  DP.B=readf(1.);
  printf(" Prism dispersion Constant C: ");
  DP.C=readf(1.);
  printf(" Input pixel size (microns) ");
  DP.tampix=readf(24.);
  printf(" Minimum wavelenght selected: ");
  minselldo=readf(minselldo);
  printf(" Maximum wavelenght selected: ");
  maxselldo=readf(maxselldo);
}


void LoadParam_file(char file[100])
{
  int status=0;
  char comment[51];
  fitsfile *parfile;

  printf("Using parameter file <<%s>>\n",file);
  if( ffopen2(&parfile,file, READONLY, &status)) fits_report_error(stderr,status);
  status=0;
  ffgky(parfile,TINT,"XCOL",&xcol,comment,&status);
  ffgky(parfile,TINT,"YCOL",&ycol,comment,&status);
  ffgky(parfile,TINT,"RACOL",&racol,comment,&status);
  ffgky(parfile,TINT,"DECCOL",&deccol,comment,&status);
  ffgky(parfile,TINT,"SKYCOL",&skycol,comment,&status);
  ffgky(parfile,TINT,"NAMECOL",&speccol,comment,&status);
  ffgky(parfile,TINT,"ENAMECOL",&errspeccol,comment,&status);
  ffgky(parfile,TSTRING,"OBJFILE",specfile,comment,&status);
  ffgky(parfile,TSTRING,"OUTFILE",outfile,comment,&status);
  ffgky(parfile,TSTRING,"SELFILE",selecfile,comment,&status);
  ffgky(parfile,TSTRING,"RESPFILE",respfile,comment,&status);
  ffgky(parfile,TFLOAT,"MINEW",&minselew,comment,&status);
  ffgky(parfile,TFLOAT,"MAXEW",&maxselew,comment,&status);
  ffgky(parfile,TFLOAT,"SEEING",&seeing,comment,&status);
  ffgky(parfile,TFLOAT,"SECPPIX",&platescale,comment,&status);
  ffgky(parfile,TFLOAT,"ERRXSP",&sigmax,comment,&status);
  ffgky(parfile,TFLOAT,"A_DISP",&DP.A,comment,&status);
  ffgky(parfile,TFLOAT,"B_DISP",&DP.B,comment,&status);
  ffgky(parfile,TFLOAT,"C_DISP",&DP.C,comment,&status);
  ffgky(parfile,TFLOAT,"PIXSIZE",&DP.tampix,comment,&status);
  ffgky(parfile,TFLOAT,"MINLDO",&minselldo,comment,&status);
  ffgky(parfile,TFLOAT,"MAXLDO",&maxselldo,comment,&status);
  if(status!=0) {
    fits_report_error(stderr,status);
    printf("ERROR: Reading parameter file.\nERROR: Not enough keywords or badformed parameter file\n");
    exit(1);
  }
  fits_close_file(parfile,&status);

}


void CreateListSpec(void) {

  int i;
  int nspec;
  int *islog;
  int *ieslog;
  int *ilog;

  nspec=FileNLin(specfile);

  printf(" Number of spectra: %d\n",nspec);

  islog=vector_i(nspec);
  ieslog=vector_i(nspec);
  ilog=vector_i(nspec);
  if((speclist=malloc(nspec*101*sizeof(char)))==NULL) {
    printf("I cannot dimension speclist of %d bytes \n",nspec*sizeof(char));
    exit(1);
  }
  if((errspeclist=malloc(nspec*101*sizeof(char)))==NULL) {
    printf("I cannot dimension errspeclist of %d bytes \n",nspec*sizeof(char));
    exit(1);
  }
  xspec=vector_f(nspec);
  yspec=vector_f(nspec);
  skyspec=vector_f(nspec);
  raspec=vector_d(nspec);
  decspec=vector_d(nspec);

  for(i=0;i<nspec;i++) {islog[i]=1;ieslog[i]=1;ilog[i]=1;}

  ReadCharcol(specfile,speccol   ,speclist ,islog ,101,&nspec);
  ReadCharcol(specfile,errspeccol,errspeclist,ieslog,101,&nspec);
  ReadNumcol(specfile,      xcol,xspec ,  ilog,&nspec);
  ReadNumcol(specfile,      ycol,yspec ,  ilog,&nspec);
  ReadNumcol(specfile,    skycol,skyspec ,  ilog,&nspec);
  ReadDoublecol(specfile,    racol,raspec ,  ilog,&nspec);
  ReadDoublecol(specfile,   deccol,decspec ,  ilog,&nspec);
  nspectra=0;
 
  nspectra=nspec;
  for(i=0;i<nspectra;i++) {
    if(islog[i] && ieslog[i]) {
      if(DEBUG2) printf(" Espectro %s\n",speclist+i*101);
    }
    else {
      if(DEBUG2) printf(" Borrando %d\n",i);
      DeleteRecord_s(speclist,nspectra,101,i);
      DeleteRecord_s(errspeclist,nspectra,101,i);
      DeleteRecord_f(xspec,nspectra,i); 
      DeleteRecord_f(yspec,nspectra,i); 
      DeleteRecord_f(skyspec,nspectra,i); 
      DeleteRecord_d(raspec,nspectra,i); 
      DeleteRecord_d(decspec,nspectra,i); 
      DeleteRecord_i(islog,nspectra,i); 
      DeleteRecord_i(ieslog,nspectra,i); 
      DeleteRecord_i(ilog,nspectra,i); 
      nspectra--;
      i--;
    }
  }

  free(islog);
  free(ieslog);
  free(ilog);
}  

void ReadAllSpectra(void) {

  int i;

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
 
  printf("%06d/%06d",0,nspectra);
  for(i=0;i<nspectra;i++) {
    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b");
    printf("%06d/%06d",i,nspectra);
    if(DEBUG2) printf(" spectra %s\n",speclist+i*101); 
    strcpy(spectra[i].file   ,speclist+i*101   );
    strcpy(errspectra[i].file,errspeclist+i*101);
    spectra[i].aloc_flag=0; 
    spectra[i].alocldo_flag=0;
    errspectra[i].aloc_flag=0;
    errspectra[i].alocldo_flag=0;
    ReadSpec_buf(   spectra+i);
    ReadSpec_buf(errspectra+i);

    if(npix==0) npix=spectra[i].npixels;
    else if(npix!= spectra[i].npixels || npix!= errspectra[i].npixels)  {
      nomatchingpixels_flag=1;
    }
  }
}

 
void ReadResp(char respfile[])
{
  int i;
  if(strcmp(respfile,"NONE")) {
    response.aloc_flag=0;
    response.alocldo_flag=0;
    strcpy(response.file,respfile);
    ReadSpec(&response);
  }
  if(DEBUG) {
    for(i=0;i<response.nx;i++) {
      printf(" %d ldoresp %f resp %f resp int. %f %f\n",i,response.ldo[i],response.spec[i],resp(response.spec[i]),Lagr4(response.ldo,response.spec,response.nx,response.ldo[i]));  
    }
  }
}


float resp(float ldo) {
  return(Lagr2(response.ldo,response.spec,response.nx,ldo));
}



void PCA(void) {
  float med[nspectra];
  float *cov,*cor;
  float *evec;
  float *esp;
  float eval[npix];
  int i,j;
  int nnllddoo;
  float *espec;
  float *lambda;
  float ldomin=0,ldomax=0;
  float datamin,datamax;
/*   float fnul; */
  char cnul='\0'; 
  float crpix=1.,crval=5870.4,cdelt=28.32;
  float proj1[nspectra],proj2[nspectra],proj3[nspectra];
  float xmin,xmax,ymin,ymax;
  float mindist=1.e15,dist;
  float xcur=0,ycur=0; 
  int jsel=0; 

  double norma;

  if(nomatchingpixels_flag) {
    printf(" All the spectra are not equal. Giving up PCA. Exiting\n");
    exit(1);
  }

  esp=vector_f(npix*nspectra);
  for(j=0 ;j<nspectra;j++) {
    printf(" Por espectro %d\n",j);
    norma=0;
    for(i=0 ;i<npix;i++) {
      norma+=spectra[j].spec[i];
    }
    for(i=0 ;i<npix;i++) {
      esp[i+npix*(j-1)]=spectra[j].spec[i]/norma*npix;
      esp[i+npix*(j-1)]=spectra[j].spec[i];
      printf(" esp %f \n",esp[i+npix*(j-1)]);
    }
  }
  

  printf(" Npix vale %d\n",npix);
  if((cov=malloc(npix*npix*sizeof(float)))==NULL) printf("I cannot dimension cov of %d elements \n",npix*npix);
  if((cor=malloc(npix*npix*sizeof(float)))==NULL) printf("I cannot dimension cor of %d elements \n",npix*npix);
  if((evec=malloc(npix*npix*sizeof(float)))==NULL) printf("I cannot dimension evec of %d elements \n",npix*npix);
  if((espec=malloc(npix*sizeof(float)))==NULL) printf("I cannot dimension espec of %d elements \n",npix);
  if((lambda=malloc(npix*sizeof(float)))==NULL) printf("I cannot dimension lambda of %d elements \n",npix);

  printf("Aqui\n");
  covcor(esp,nspectra,npix,med,cov,cor);
  printf("Aqui 2\n");
  if(DEBUG) {
    for(i=0;i<npix;i++) {
      for(j=0;j<i+1;j++) {
	printf(" %f ",cov[i+npix*j]);
      }
      printf(" \n");
    }
  }
  eigen(cov,npix,eval,evec);
/*   //eigen(cor,npix,eval,evec); */
  printf(" Npix vale %d\n",npix);
  nnllddoo=npix;
/*   //eigen(cov,50,eval,cov); */
  printf("Aqui 3\n");

  for(i=0;i<npix;i++) {
    printf(" autoval %d: %g \n",i,eval[i]);
    printf(" autovec %d: ",i);
    for(j=0;j<npix;j++) {
      espec[j]=evec[j+npix*i];
      lambda[j]=(j+1-crpix)*cdelt+crval;
      printf(" %f ",evec[j+npix*i]);
    }
    ldomin=(1-crpix)*cdelt+crval;
    ldomax=(npix-crpix)*cdelt+crval;
    datamin=1.0e30;
    datamax=-1.0e30;
    for (j=0 ;j<npix;j++) {
      if(espec[j]< datamin) datamin=espec[j];
      if(espec[j]> datamax) datamax=espec[j];
    }
    cpgpage();
    cpgswin(ldomin,ldomax,datamin-(datamax-datamin)*.2,datamax+(datamax-datamin)*.2);
    cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
    cpglab("Wavelength \\(2078)","Flux","");
    cpgline(npix,lambda,espec);
/*     //cpgcurs(&fnul,&fnul,&cnul); */

    
    printf("\n");
  }

  for(j=0;j<nspectra;j++) {
    proj1[j]=0;
    proj2[j]=0;
    proj3[j]=0;
    for(i=0;i<npix;i++) {
      proj1[j]+=evec[i+npix*0]*esp[i+npix*j];
      proj2[j]+=evec[i+npix*1]*esp[i+npix*j];
      proj3[j]+=evec[i+npix*2]*esp[i+npix*j];
    }
  }
  xmin=1.0e30;ymin=1.0e30;
  xmax=-1.0e30;ymax=-1.0e30;
  for (j=0 ;j<nspectra;j++) {
    if(proj1[j]< xmin) xmin=proj1[j];
    if(proj1[j]> xmax) xmax=proj1[j];
    if(proj2[j]< ymin) ymin=proj2[j];
    if(proj2[j]> ymax) ymax=proj2[j];
  }

  while(!(cnul=='m')) {
    mindist=1.e15;
    cpgpage();
    cpgswin(xmin*1.1,xmax*1.1,ymin*1.1,ymax*1.1);
    cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
    cpglab("Projection on eigenvector 1","Projection on eigenvector 2","");
    cpgpt(nspectra,proj1,proj2,1);
    cnul='m';
/*     //cpgcurs(&xcur,&ycur,&cnul); */
    for(j=0;j<nspectra;j++) {
      dist=((proj1[j]-xcur)*(proj1[j]-xcur)+(proj2[j]-ycur)*(proj2[j]-ycur));
      if(mindist*mindist>dist*dist) {
	mindist=dist;
	jsel=j;
      }
    }
    
    for(j=0;j<npix;j++) {
      espec[j]=esp[j+npix*jsel];
    }
    
    datamin=1.0e30;
    datamax=-1.0e30;
    for (j=0 ;j<npix;j++) {
      if(espec[j]< datamin) datamin=espec[j];
      if(espec[j]> datamax) datamax=espec[j];
    }
    cpgpage();
    cpgswin(ldomin,ldomax,datamin-(datamax-datamin)*.2,datamax+(datamax-datamin)*.2);
    cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
    cpglab("Wavelength \\(2078)","Flux","");
    cpgline(npix,lambda,espec);
    printf(" Proj1 %f Proj2 %f\n",proj1[jsel],proj2[jsel]);
    
    for(j=0;j<npix;j++) {
      espec[j]=evec[j+npix*0]*proj1[jsel];
    }
    cpgsci(2);
    cpgline(npix,lambda,espec);
    for(j=0;j<npix;j++) {
      espec[j]=evec[j+npix*0]*proj1[jsel]+evec[j+npix*1]*proj2[jsel];
    }
    cpgsci(3);
    cpgline(npix,lambda,espec);
    for(j=0;j<npix;j++) {
      espec[j]=evec[j+npix*0]*proj1[jsel]+evec[j+npix*1]*proj2[jsel]+evec[j+npix*2]*proj3[jsel];
    }
    cpgsci(4);
    cpgline(npix,lambda,espec);
    cpgsci(1);

/*     //cpgcurs(&fnul,&fnul,&cnul); */
  }
    
  if(!strcmp(respfile,"NONE")) {
    response.spec=vector_f(npix);
    response.ldo=vector_f(npix);
    response.nx=npix;
    response.naxes[1]=npix;
    response.naxes[0]=1;
    response.aloc_flag=1;
    response.alocldo_flag=1;
    for(i=0;i<response.nx;i++) {
      response.spec[i]=evec[i+npix*0];
      response.ldo[i]=lambda[i];
    }
  }
  printf(" Ha salido de PCA\n");

}

 

int FitSpec(struct spectrum spec, struct spectrum errspec,float *ew, float *errew, float *lpix, float *errlpix, float *K, float *errK) {

  int i;

/*   float first,third,min; */
/*   float modespec,moderesp; */
/*   float maxresp,maxspec; */
  struct spectrum specshifted;
  struct spectrum errspecshifted;
  float shift;
  float sumspec;
  int selected;

 

  if(DEBUG) printf(" sumspec %f sumresp %f\n",sumspec,sumresp);

  specshifted.aloc_flag=0;
  specshifted.alocldo_flag=0;
  CopySpec(&specshifted,spec);
  if(DEBUG) printf(" after first copy\n");
  errspecshifted.aloc_flag=0;
  errspecshifted.alocldo_flag=0;
  CopySpec(&errspecshifted,errspec);
  if(DEBUG) printf(" Before shift\n");
  shift=offsetpix(spec.nx,spec.spec,response.spec);
  if(DEBUG) printf(" shift %f\n",shift);
  if(fabs(shift<10) && fabs(shift)>sigmax/4) {
    shiftspec(-shift,spec.nx,spec.spec,specshifted.spec);
    shiftspec(-shift,errspec.nx,errspec.spec,errspecshifted.spec);
  }
  if(DEBUG) printf(" Despues del shift\n");
  sumspec=StSuma1((int)uppix-(int)lowpix+1,spec.spec+(int)lowpix,1);
  if(DEBUG) printf(" Calculo norma %f\n",sumspec); 
  for(i=0;i<specshifted.nx;i++) {
    specshifted.spec[i]*= sumresp/sumspec;
    errspecshifted.spec[i]*= sumresp/sumspec;
  }

  selected=Minimize(specshifted,errspecshifted,ew,errew,lpix,errlpix,K,errK);

  CloseSpec(&specshifted); 
  CloseSpec(&errspecshifted);
  return(selected);
}


int Minimize(struct spectrum spec, struct spectrum errspec, float *ew, float *errew, float *lpix, float *errlpix, float *K, float *errK) {

  int iew,iz,i;
  int new=15,nz=25;
  struct spectrum sintec;
  struct spectrum errspectot;
  struct spectrum specleft;
  struct spectrum specright;
  double **ML;
  float dpdl;
  float ldoline;
  double mlfit=1e39;
  int iewmin=0,izmin=0;
  double mlrandom,ml2random,ml0,sigmlrandom;
  double mlleft,mlright;
  int npixml;
  int iter;
  double par[3],sigpar[3],**covarpar;
  float sigmatol;
  int selected=0;
  int errstat;
  double *xdob,*ydob;

  xdob=vector_d(spec.nx);  ydob=vector_d(spec.nx);
  covarpar=matrix_d(3,3);

  sintec.aloc_flag=0;
  sintec.alocldo_flag=0;
  errspectot.aloc_flag=0;
  errspectot.alocldo_flag=0;
  ML=matrix_d(new,nz);
  CopySpec(&errspectot,errspec);
  for (i=0;i<spec.nx;i++) {
    errspectot.spec[i]=sqrt(errspec.spec[i]*errspec.spec[i]+errfromresp.spec[i]*errfromresp.spec[i]);
    errspectot.spec[i]=errspec.spec[i];
  }
    
  
  for(iew=0;iew<new;iew++) {
    *ew=exp(-0.3+(log(maxselew)+0.3)*iew/(new-1));
    for(iz=0;iz<nz;iz++) {
      *lpix=lowpix+(uppix-lowpix)*iz/(nz-1);
      
      ldoline=pr_pix2ldo(*lpix,DP);
      dpdl=pr_dpdl(ldoline,DP);
      if(DEBUG2) printf(" Ldo line %f EW (pix) %f\n",ldoline,*ew*dpdl); 
	     
      CreateSintec(&sintec,*ew*dpdl,*lpix,seeing/platescale,0.); 
      if(DEBUG2) PlotSpec_pix_err(spec,errspec);
      if(DEBUG2) PlotSpec_ov(sintec);
      ML[iew][iz]=CalcML(spec,errspectot,sintec);
      if(ML[iew][iz]<mlfit) {
	mlfit=ML[iew][iz];
	iewmin=iew;
	izmin=iz;
      }
      if(DEBUG2) printf(" EW %f lpix %f ML %f\n",*ew,*lpix,ML[iew][iz]);
      if(DEBUG2) iz=readi(iz);
    }
  }
  CreateSintec(&sintec,0.,*lpix,seeing/platescale,0.); 
  ml0=CalcML(spec,errspectot,sintec);
  mlrandom=0;
  /* Esto de aqui es el valor teorico de L para una distribucion aleatoria  
   Es la media del logarimo de de una gaussiana para una distribucin gaussiana.
  Y luego multiplicadas todas */
  npixml=0;
  for (i=0;i<spec.nx;i++) {
    if(response.spec[i]<0.08) ;
    else {
      npixml++;
      mlrandom -= -log(errspectot.spec[i]*sqrt(2*M_PI)) - 0.5;  
      ml2random -= log(errspectot.spec[i]*sqrt(2*M_PI)) +log(errspectot.spec[i]*sqrt(2*M_PI)) * log(errspectot.spec[i]*sqrt(2*M_PI)) + 0.75;
    }  
  }
  sigmlrandom=sqrt(npixml/2.);
  *ew=exp(-0.3+(log(maxselew)+0.3)*iewmin/(new-1));
  *lpix=lowpix+(uppix-lowpix)*izmin/(nz-1);
  CreateSintec(&sintec,*ew*dpdl,*lpix,seeing/platescale,0.); 
  if(DEBUG) PlotSpec_pix_err(spec,errspectot);
  if(DEBUG2) PlotSpec_pix_ov(sintec);
  
  if(DEBUG) printf(" AJUSTE EW %f lpix %f ML %f\n",*ew,*lpix,mlfit);
  if(DEBUG) printf(" mlrandom %f sigma %f\n",mlrandom,sigmlrandom);
  if(DEBUG) printf(" ml EW=0 : %f\n",ml0);
  if(fabs(mlfit-ml0)<2.0*sigmlrandom) {
    if(DEBUG) printf(" Compatible con EW=0\n");
  }
  else {
    specamo=&spec;
    errspecamo=&errspectot;
    par[0]=*ew;sigpar[0]=*ew/2;
    par[1]=*lpix;sigpar[1]=2.;
    par[2]=0.;sigpar[2]=0.1/npixml;
    iter=Amoeba_d(spec.nx,xdob,ydob,3,par,sigpar,FTOL,NITER,amoe_fit_sl);
    if(DEBUG) printf(" Sale de amoeba con ew=%f lpi=%f K=%f \n",par[0],par[1],par[2]);
    *ew=par[0];
    *lpix=par[1];
    *K=par[2];
    ldoline=pr_pix2ldo(*lpix,DP);
    dpdl=pr_dpdl(ldoline,DP);
    mlfit=amoe_fit_sl(spec.nx,xdob,ydob,par);
    CreateSintec(&sintec,*ew*dpdl,*lpix,seeing/platescale,*K); 
    if(DEBUG2) PlotSpec_pix_err(spec,errspectot);
    if(DEBUG) PlotSpec_pix_ov(sintec);
    specleft.aloc_flag=0;
    specleft.alocldo_flag=0;
    specright.aloc_flag=0;
    specright.alocldo_flag=0;
    CopySpec(&specleft,spec);
    CopySpec(&specright,spec);
    shiftspec(-0.5,spec.nx,spec.spec,specleft.spec);
    shiftspec(+0.5,spec.nx,spec.spec,specright.spec);
    mlleft=CalcML(specleft,errspectot,sintec);
    mlright=CalcML(specright,errspectot,sintec);
    if(DEBUG) {
      printf(" Amoeba mlfit %f\n",mlfit);
      printf(" ML2Fit %f\n",CalcML2(spec,errspectot,sintec));
      printf(" mlleft %f mlright %f\n",mlleft,mlright);
      mlleft=CalcML2(specleft,errspectot,sintec);
      mlright=CalcML2(specright,errspectot,sintec);
      printf(" CON ML2 mlleft %f mlright %f\n",mlleft,mlright);
      mlleft=CalcML(specleft,errspectot,sintec);
      mlright=CalcML(specright,errspectot,sintec);
    }
    /*     if(fabs(mlleft-mlfit)<3*sigmlrandom && fabs(mlright-mlfit)<3*sigmlrandom && fabs(mlrandom-mlfit)> 5*sigmlrandom) { */
    if(fabs(mlleft-mlfit)<0.2*fabs(mlrandom-mlfit) && fabs(mlright-mlfit)<0.2*fabs(mlrandom-mlfit) && fabs(mlrandom-mlfit)> 5*sigmlrandom) {
      if(DEBUG) printf(" La funcion response no es muy buena\n");
    }
    else {
      if(DEBUG) printf(" Buena funcion respuesta\n");
/*       if(fabs(mlrandom-mlfit)> 10*sigmlrandom) { */
/* 	if(DEBUG) printf(" No es un buen ajuste compatible con random\n"); */
/* 	sigmatol=4; */
/*       } */
/*       else { */
/* 	if(DEBUG) printf(" ES un buen ajuste en plan random\n"); */
/* 	sigmatol=3; */
/*       } */
      sigmatol=fabs(mlrandom-mlfit)/sigmlrandom;
      if(sigmatol<3) sigmatol=3;
      if(DEBUG) printf(" sigmatol %f \n",sigmatol);

      sigpar[0]=par[0]/30.;
      sigpar[1]=0.5;
      sigpar[2]=par[2]/5.;
      covarpar[0][0]=-1;
      covarpar[1][1]=-1;
      errstat=mlerr_amo_d(spec.nx,xdob,ydob,3,par,sigpar,FTOLERR,NITERERR,amoe_fit_sl,NCONFL,covarpar);
      if(covarpar[0][0]<0 || covarpar[1][1]<0) errstat=mlerr_amo_d(spec.nx,xdob,ydob,3,par,sigpar,FTOLERR,NITERERR,amoe_fit_sl,NCONFL,covarpar); 
      if(covarpar[0][0]<0 || covarpar[1][1]<0) errstat=mlerr_amo_d(spec.nx,xdob,ydob,3,par,sigpar,FTOLERR,NITERERR,amoe_fit_sl,NCONFL,covarpar);
      if(DEBUG) printf(" Esto antes ultimo %f %f\n",covarpar[0][0],covarpar[1][1]);
      if(covarpar[0][0]<0 || covarpar[1][1]<0|| (errstat==0)) {
	if(DEBUG) printf(" I was not able to compute errores\n");
      }
      else {
	if(DEBUG) printf(" Errores: ew %f +/- %f lpix %f +/- %f K %f +/- %f\n",par[0],sqrt(covarpar[0][0]),par[1],sqrt(covarpar[1][1]),par[2],sqrt(covarpar[2][2]));
	*errew=sqrt(covarpar[0][0]);
	*errlpix=sqrt(covarpar[1][1]);
	*errK=sqrt(covarpar[2][2]);
	
	if(fabs(*ew-minselew)>sigmatol**errew && *ew> minselew) {
	  if(DEBUG) printf(" Va por el buen camino. Segundo test\n");
          par[0]=*ew;sigpar[0]=*errew*2;
          par[1]=*lpix;sigpar[1]=*errlpix*2;
          par[2]=*K;sigpar[2]=*K/20;
          iter=Amoeba_d(spec.nx,xdob,ydob,3,par,sigpar,FTOL2,NITER2,amoe_fit_sl2);
          *ew=par[0];*lpix=par[1];*K=par[2];
          covarpar[0][0]=-1;
          covarpar[1][1]=-1;
          errstat=mlerr_amo_d(spec.nx,xdob,ydob,3,par,sigpar,FTOLERR2,NITERERR2,amoe_fit_sl2,NCONFL2,covarpar);
          if(covarpar[0][0]<0 || covarpar[1][1]<0) mlerr_amo_d(spec.nx,xdob,ydob,3,par,sigpar,FTOLERR2,NITERERR2,amoe_fit_sl2,NCONFL2,covarpar);
          if(covarpar[0][0]<0 || covarpar[1][1]<0) errstat=mlerr_amo_d(spec.nx,xdob,ydob,3,par,sigpar,FTOLERR2,NITERERR2,amoe_fit_sl2,NCONFL2,covarpar);
          if(covarpar[0][0]<0 || covarpar[1][1]<0|| (errstat==0)) {
	     if(DEBUG) printf(" I was not able to compute second errores\n");
          }
          else {
	    /* Aqui reestablezco el nivel para sigmatol */
	    sigmatol=4;
 	    if(DEBUG) printf(" Errores2: ew %f +/- %f lpix %f +/- %f K %f +/- %f\n",par[0],sqrt(covarpar[0][0]),par[1],sqrt(covarpar[1][1]),par[2],sqrt(covarpar[2][2]));
	    *errew=sqrt(covarpar[0][0]);
	    *errlpix=sqrt(covarpar[1][1]);
            ldoline=pr_pix2ldo(*lpix,DP);
            if(DEBUG) printf(" ldo %f\n",ldoline);
	    if(fabs(*ew-minselew)>sigmatol**errew && *ew> minselew && ldoline>minselldo && ldoline<maxselldo) {
	      if(DEBUG) printf(" SIIIIIII. Este es candidato\n");
	      if(sigmatol==5) selected=2;
	      else            selected=1;
	      if(DEBUG2) mlfit=readf(mlfit);
            }
          }
	}
	else {
          if(DEBUG) printf(" Se quedo fuera\n");
        }
      }
    }
    if(DEBUG) printf(" FInal \n");
    CloseSpec(&specleft);
    CloseSpec(&specright);
  }
  
  if(DEBUG) {
    mlfit=readf(mlfit); 
    if(selected) {
      PlotSpec_pix_err(spec,errspectot); 
      PlotSpec_pix_ov(sintec);
    }
  }

  free(xdob);free(ydob);
  free_matrix_d(covarpar,3,3);
  free_matrix_d(ML,new,nz);
  CloseSpec(&sintec);
  CloseSpec(&errspectot);
  return(selected);
}


double CalcML(struct spectrum es, struct spectrum erres, struct spectrum si) {

  int i;

  double ml=0;
  double K;

  K=-log(sqrt(2.*M_PI));
  for (i=0;i<es.nx;i++) {
/*     printf(" es %f si %f err %f  arg %f ga %f ml %f\n",es.spec[i],si.spec[i],erres.spec[i],(es.spec[i]-si.spec[i])/erres.spec[i],gaussian(es.spec[i],si.spec[i],erres.spec[i]),ml);  */
    if(response.spec[i]<0.08) ;
    else ml -= K - log(erres.spec[i]) - (es.spec[i]-si.spec[i])*(es.spec[i]-si.spec[i])/2./erres.spec[i]/erres.spec[i];
    /*     ml - = log(gaussian(es.spec[i],si.spec[i],erres.spec[i])); */
  }
  return(ml);
}

void CreateSintec(struct spectrum *sintec, float ew_pix, float lpix, float fwhm, float K) {


  float sigma;
  int i;
  float sumsintec;

  CopySpec(sintec,response);
  sigma=fwhm/2.35;
  for(i=0;i<response.nx;i++) {
    (*sintec).spec[i]=response.spec[i]*(1+K*i+ew_pix*gaussian((float)i,lpix,sigma));
  }
  sumsintec=StSuma1((int)uppix-(int)lowpix+1,(*sintec).spec+(int)lowpix,1);
  for(i=0;i<response.nx;i++) (*sintec).spec[i]*=sumresp/sumsintec;

}

double amoe_fit_sl(int n, double *x, double *y, double *p) {

  double ew,lpix,K;
  double ldoline,dpdl;
  struct spectrum sintec;
  double ml;
  
  ew=p[0];
  lpix=p[1];
  K=p[2];

  if(ew>maxselew*2) p[0]=maxselew;
  
  ldoline=pr_pix2ldo(lpix,DP);
  dpdl=pr_dpdl(ldoline,DP);

  sintec.aloc_flag=0;
  sintec.alocldo_flag=0;
  CreateSintec(&sintec,ew*dpdl,lpix,seeing/platescale,K); 
  
/*   printf(" aloc %d\n",sintec.aloc_flag); */
/*   printf(" alocldo %d\n",sintec.alocldo_flag); */
  
  ml=CalcML(*specamo,*errspecamo,sintec);
  
/*   printf(" AMo ew = %f  lpix = %f  K = %f   ml=%f\n",ew,lpix,K,ml);  */
  
  CloseSpec(&sintec);

  return(ml);
}

double amoe_fit_sl2(int n, double *x, double *y, double *p) {

  double ew,lpix,K;
  double ldoline,dpdl;
  struct spectrum sintec;
  double ml;
  
  ew=p[0];
  lpix=p[1];
  K=p[2];
  
  ldoline=pr_pix2ldo(lpix,DP);
  dpdl=pr_dpdl(ldoline,DP);

  sintec.aloc_flag=0;
  sintec.alocldo_flag=0;
  CreateSintec(&sintec,ew*dpdl,lpix,seeing/platescale,K); 
  
  ml=CalcML2(*specamo,*errspecamo,sintec);
  
/*   printf(" AMo ew = %f  lpix = %f  K = %f   ml=%f\n",ew,lpix,K,ml);  */
  
  CloseSpec(&sintec);

  return(ml);
}

double CalcML2(struct spectrum es, struct spectrum erres, struct spectrum si) {

  int i,j;


  double ml=0;
  double K;
  double sigmaf_x,xmin,xmax,x1,x2,f_x,ymin,ymax;
  double intg,gi;
  //double intg1,gi1,ml1=0;

  K=-log(sqrt(2.*M_PI));
  for (i=0;i<es.nx;i++) {
    /*     printf(" es %f si %f err %f  arg %f ga %f ml %f\n",es.spec[i],si.spec[i],erres.spec[i],(es.spec[i]-si.spec[i])/erres.spec[i],gaussian(es.spec[i],si.spec[i],erres.spec[i]),ml);  */
    if(response.spec[i]<0.08) ;
    else {
      intg=0;
      //intg1=0;
      sigmaflux=erres.spec[i];
      fluxorig=es.spec[i];
      for(j=-(int)(3*sigmax)-1;j<+(int)(3*sigmax)+1;j++) {
        ymin=si.spec[i+j];ymax=si.spec[i+j+1];
        dfdx=ymax-ymin;
	ybegin=ymin;
	x0=(double)(-j);
	if(dfdx==0) {
          xmin=x0-3*sigmax;
          xmax=x0+3*sigmax;
          x1=maxf(xmin,0.); if(x1>=1) x1=0.;
          x2=minf(xmax,1.); if(x2<=0) x2=1.;

	}
	else {
	  sigmaf_x=erres.spec[i]/fabs(dfdx);
          if(sigmaf_x<sigmax) {
	    f_x=(es.spec[i]-ybegin)/dfdx;
	    xmin=f_x-3*sigmaf_x;
	    xmax=f_x+3*sigmaf_x;
	    //printf("   PRIMERO x %f %f f_x %f s_x %f y %f %f fl %f efl %f\n",xmin,xmax,f_x,sigmaf_x,ymin,ymax,es.spec[i],erres.spec[i]);
	    x1=maxf(xmin,0.); if(x1>=1) x1=0.;
	    x2=minf(xmax,1.); if(x2<=0) x2=1.;
          }
          else {
	    xmin=x0-3*sigmax;
	    xmax=x0+3*sigmax;
	    //printf("   SEGUNDO x %f %f f_x %f s_x %f y %f %f fl %f efl %f\n",xmin,xmax,f_x,sigmaf_x,ymin,ymax,es.spec[i],erres.spec[i]);
	    x1=maxf(xmin,0.); if(x1>=1) x1=0.;
	    x2=minf(xmax,1.); if(x2<=0) x2=1.;
          }
	}
	gi=gaussintleg_d(Funk_int_ml,(double)x1,(double)x2,7);
        intg+=gi;
	//printf(" pix %d  j %d x %f-%f gi %f intg %f\n",i,j,x1,x2,gi,intg);
	//x1=0.;x2=1.;
	//dfdx=0.;
	//ybegin=si.spec[i];
	//gi1=gaussintleg_d(Funk_int_ml,(double)x1,(double)x2,7);
        //intg1+=gi1;
	//printf("  ml1pix %d  j %d x %f-%f gi1 %f intg1 %f\n",i,j,x1,x2,gi1,intg1);
      }
      ml -= log(intg);
      //ml1 -= K - log(erres.spec[i]) - (es.spec[i]-si.spec[i])*(es.spec[i]-si.spec[i])/2./erres.spec[i]/erres.spec[i];
    }
    //printf(" PIX %d int g %f mlold %f mlold(int) %f ml %f ml1 %f\n",i,log(intg),K - log(erres.spec[i]) - (es.spec[i]-si.spec[i])*(es.spec[i]-si.spec[i])/2./erres.spec[i]/erres.spec[i],log(intg1),ml,ml1);
    //printf(" f0 %f sig %f\n",fluxorig,sigmaflux);
  }
  return(ml);
}


double Funk_int_ml(double x) {
  double fluxsintec;
  fluxsintec=ybegin+dfdx*x;
  //printf("  x %f x0 %f fo %f fs %f fac1 %f fac2 %f\n",x,x0,fluxorig,fluxsintec,gaussian(fluxorig,fluxsintec,sigmaflux),gaussian(x,x0,sigmax));
  return(gaussian(fluxorig,fluxsintec,sigmaflux)*gaussian(x,x0,sigmax));
}
