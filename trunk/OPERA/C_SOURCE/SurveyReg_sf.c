#include "modulos.h"
#define DEBUG 0
#define DEBUG2 0

#define FTOL 1.e-6
#define NMON 1000

/* Common variables */
/* Parameter file */
char parfilename[100]="";
int interact=0;
int amoegauss=0;
int asstatus=1,mcstatus=1,mlstatus=1,sestatus=1,specstatus=1;

/* Parameters of the program */
char imagefile[51];
char specfile[51];
char respfile[51];
char filecross[51];
int  fluxcol;
int  errfluxcol;
int  magcol;
int  xcol,ycol;
char filereg[51];
char pgdevice[51];
float platescale;
float fluxsatura;

/*data matrix */
float *flux,*errflux,*magni,*xobj,*yobj;
int *ilog;
int ndat;

/*Flux calibration variables. mag = a + b log10(flux)*/
double a,b;
double erra,errb,covab;
float rms;
static double sigma=1.e29;
static float nsigrej=3;

/*Astrometric information */
double alfac, deltac;
double xdim, ydim;
double rot;
double epoch;
double aststdx,aststdy;

/*Limiting magnitude  Numbre counts distribution variables*/
float mmaximum;
float mminimum;
float mmode;
float mfermicut;
float deltafermi;
float gammapowerlaw;
float errmfermicut;
float errdeltafermi;
float errgammapowerlaw;
float *xglob,*yglob;
double *xglob_d,*yglob_d;
int nbin=50;
float deltabin;


/*Seeing variables */
float seeing;
float stdseeing;
float *xpglob,*cutyglob;
int nspacut;

void LoadParam_kbd(void);
void LoadParam_file(void);
void SaveParam(void);
int  AstromInfo(void);
void FluxCal_Basic(double *logcts,double *errlogcts,double *m,double *resi,int *nfot);
int  FluxCalibration(void);
int  MagLim(void);
int  Seeing(void);
int  SpecInfo(void);
void ReadCat(void);
void PrintLog(void);
void PlotMagCal(double *logcts,double *m,double *errlogcts,int nfot);
/* float Amoe_Funk(int n, float *x, float *y, float *p); */
float logFermi(float x,float mu,float T);
float PowLaw(float x,float A,float x0,float gam);
float DerPowLaw(float x,float A,float x0,float gam);
float logDerPowLaw(float x,float A,float x0,float gam);
float Fselerf(float x,float A,float B);
void funmaglim_d(double x,double *p,double *y,double *dyda,int n);
void funseeing(float x,float *p,float *y,float *dyda,int n);

void funlin(float x,float *y,int n);

int main(int argc, char **argv)
{

  char option;


/*    char filebuho[51]; */
/*    int fluxcol=4,magcol=30; */

/*    float *fluxdet,*magdet; */
/*    int fluxdetcol=4; */
/*    float fluxmin,fluxmax,magmin,magmax; */
/*    float a,b; */
/*    int i,ii; */
/*    int nobj,ntot,ndetect,ntotdetect; */
/*    float sigma; */


  if(argc<2) {
    LoadParam_kbd();
    interact=1;
  }
  else {
    strcpy(parfilename,argv[1]);
    LoadParam_file(); 
  }

  asstatus=AstromInfo();
  printf(" ra %f dec %f xdim %f ydim %f rot %f\n",alfac,deltac,xdim,ydim,rot); 
  ReadCat();
  
  printf(" >>>>> Flux calibration computation\n");

  if(ndat<10) {
    printf(" No reliable registering can be done with less than 10 points\n");
    a=0;erra=0;b=0;errb=0;covab=0;rms=0;mmaximum=0;mmode=0;mfermicut=0;errmfermicut=0;deltafermi=0;errdeltafermi=0;gammapowerlaw=0;errgammapowerlaw=0;seeing=0;stdseeing=0;
    asstatus=mcstatus=mlstatus=sestatus=0;
    PrintLog();
    if(interact) {
      printf(" Do you want to save parameters in a file? [y/n]: ");
      option=readc('y'); 
      if(option=='y')   SaveParam();
    }
    cpgend();
    return 1;
  }

  
  mcstatus=FluxCalibration(); 
  printf("m= (%f +- %f)  + (%f +- %f) * log10(Flux)\n",a,erra,b,errb);
  printf("Covar: %f\n",covab);
/*      exit(1);   */
/*   option=readc('y');   */
  printf(" >>>>> Selection function computation\n");
  mlstatus=MagLim(); 
/*   exit(1); */
/*   printf(" Donde estas \n"); */
/*   option=readc('y');   */
  
  printf(" >>>>> Seeing computation\n");
  sestatus=Seeing();  
  PrintLog();

  
  if(interact) {
    printf(" Do you want to save parameters in a file? [y/n]: ");
    option=readc('y'); 
    if(option=='y')   SaveParam();
  }

  printf("HE acabado\n");
  cpgend();

  return 0;
}



void LoadParam_kbd(void) {

  fluxcol=4;
  errfluxcol=23;
  magcol=46;
  xcol=14;
  ycol=15;

  printf(" Input FITS image to register: ");
  reads("",imagefile);
  printf(" Input file with spectra database: ");
  reads("",specfile);
  printf(" Input file with response function: ");
  reads("",respfile);
  printf(" Input catalogue with matched entries: ");
  reads("",filecross);
  printf(" Input column with x position in %s: ",filecross);
  xcol=readi(xcol);
  printf(" Input column with y position in %s: ",filecross);
  ycol=readi(ycol);
  printf(" Input column with flux in %s: ",filecross);
  fluxcol=readi(fluxcol);
  printf(" Input column with fluxerr  in %s: ",filecross);
  errfluxcol=readi(errfluxcol);
  printf(" Input column with magnitude to calibrate in %s: ",filecross);
  magcol=readi(magcol);
  printf(" Input saturation level (0=no saturation): ");
  fluxsatura=readf(0.);
  printf(" Input file where to register the survey image: ");
  reads("",filereg);
  cpgopen("?");
  cpgask(1);
}



void LoadParam_file(void) {

  int status=0;
  char comment[51];
  fitsfile *parfile;
  printf("Using parameter file <<%s>>\n",parfilename);
  if( ffopen2(&parfile,parfilename, READONLY, &status)) fits_report_error(stderr,status);

  ffgky(parfile,TSTRING,"IMAGE",imagefile,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"SPECFILE",specfile,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"RESPFILE",respfile,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"CATALOG",filecross,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"FLUXCOL",&fluxcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"EFLUXCOL",&errfluxcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"MAGCOL",&magcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"XCOL",&xcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"YCOL",&ycol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"SECPPIX",&platescale,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"SATURA",&fluxsatura,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"FILEREG",filereg,comment,&status); 
  fits_report_error(stderr,status); 
  ffgky(parfile,TSTRING,"DEVICE",pgdevice,comment,&status);
  fits_report_error(stderr,status);
  
  
  if(status!=0) {
    fits_report_error(stderr,status);
    printf("ERROR: Reading parameter file. Not enough keywords.\n");
    exit(1);
  }
  fits_close_file(parfile,&status);
  
  cpgopen(pgdevice);
  cpgask(0);
}

void SaveParam(void) {

  FILE *fp;
  int nc=0,nt;
  char ch51[51];
  printf("Name of parameter file: ");
  reads(parfilename,parfilename);
  if((fp=fopen(parfilename,"w")) ==NULL) {
    printf("ERROR: Can't open file %s\n",parfilename);
    return;
  }
  fprintf(fp,"COMMENT  Parameter file for SurveyReg                                          \n");
  nc++;
  sprintf(ch51,"'%s'",imagefile);
  fprintf(fp,"IMAGE   = %-51.51s / Image to reg.  \n",ch51);
  sprintf(ch51,"'%s'",specfile);
  fprintf(fp,"SPECFILE= %-51.51s /Spectra database\n",ch51);
  sprintf(ch51,"'%s'",respfile);
  fprintf(fp,"RESPFILE= %-51.51s /Response func.  \n",ch51);
  sprintf(ch51,"'%s'",filecross);
  fprintf(fp,"CATALOG = %-51.51s / File with objs.\n",ch51);
  fprintf(fp,"FLUXCOL =%21d / Column with CCD flux in file CATALOG          \n",fluxcol     );
  fprintf(fp,"EFLUXCOL=%21d / Column with CCD flux in file CATALOG          \n",errfluxcol     );
  fprintf(fp,"MAGCOL  =%21d / Column with calibration magnitude in CATALOG  \n",magcol      );
  fprintf(fp,"XCOL    =%21d / Column with X position in CATALOG             \n",xcol        );
  fprintf(fp,"YCOL    =%21d / Column with Y position in CATALOG             \n",ycol        );
  fprintf(fp,"SECPPIX =%21f / Platescale (arcsecs per pixel)                \n",platescale    );
  fprintf(fp,"SATURA  =%21f / Saturation level for flux (not to fit in calib\n",platescale    );
  sprintf(ch51,"'%s'",filereg );
  fprintf(fp,"FILEREG = %-51.51s / Register file  \n",ch51);
  sprintf(ch51,"'%s'","?" );
  fprintf(fp,"DEVICE  = %-51.51s / PGPLOT device  \n",ch51);
  nc+=12;
  fprintf(fp,"COMMENT  Blank lines to assure the last block of 2880 caracters is complete    \n");
  fprintf(fp,"COMMENT                                                                        \n");
  nc+=2;
  for(nt=nc;nt<(1+(int)(nc/36))*36-1;nt++) {
    fprintf(fp,"COMMENT                                                                        \n");
  }
  fprintf(fp,"END                                                                            \n");
  fclose(fp);
}







int  AstromInfo(void) {
  /* FITS variables for FITSIO*/
  fitsfile *fitsimage,*ftemp;
  long naxes[2];
  int nfound;
  int status=0;
  /*Header variables for WCSTOOLS */
  char *header;
  int lhead,nbfits;
  char comment[51];
  /* WCS variables for WCSTOOLS */
  struct WorldCoor *wcsim=0;      /* World coordinate system structure */
  int ncross;

  long rand_i;
  char tempfile[200];
  
  /* Open the FITS file and create a dummy file to read the header with fitsrhead */
  if( ffopen(&fitsimage, imagefile, READONLY, &status)) fits_report_error(stderr,status);
  if(fits_read_keys_lng(fitsimage, "NAXIS", 1, 2, naxes, &nfound, &status)) fits_report_error(stderr,status);
  if(ffgky(fitsimage,TDOUBLE,"STDX",&aststdx,comment,&status)) fits_report_error(stderr,status);
  if(ffgky(fitsimage,TDOUBLE,"STDY",&aststdy,comment,&status)) fits_report_error(stderr,status);
  if(ffgky(fitsimage,TINT,"NCROSS",&ncross,comment,&status)) fits_report_error(stderr,status);
  aststdx=aststdx/sqrt((float)ncross);
  aststdy=aststdy/sqrt((float)ncross);
  rand_i=random();
  sprintf(tempfile,"temp%012ld_header.fits",rand_i);
  unlink(tempfile);
  ffinit(&ftemp,tempfile,&status);
  fits_copy_header(fitsimage, ftemp, &status);
  ffclos(fitsimage,&status);
  ffclos(ftemp,&status);
  if ((header = fitsrhead (tempfile, &lhead, &nbfits)) == NULL) {
    fprintf (stderr, "Cannot read FITS header of file %s\n", imagefile);
    exit(1);
  }
  unlink(tempfile);

  /* Initialize WCS structure with header keywords */
  if((wcsim=wcsinit(header))==NULL) {
    printf(" No WCS information found in header\n");  
    alfac=0;
    deltac=0;
    rot=0;
    epoch=0;
    xdim=0;
    ydim=0;
    return(0);
  }
  
  pix2wcs(wcsim,naxes[0]/2., naxes[1]/2. ,&alfac, &deltac);
  xdim=(double)(naxes[0]*fabs((float)((*wcsim).xinc)*3600.));
  ydim=(double)(naxes[1]*fabs((float)((*wcsim).yinc)*3600.));
  rot=(*wcsim).rot;
  epoch=(*wcsim).epoch;
  return(1);
}







void ReadCat(void) {

  ndat=FileNLin(filecross);
  

  if((flux=malloc(ndat*sizeof(float)))==NULL) printf("I cannot dimension flux  of %d elements \n",ndat);
  if((errflux=malloc(ndat*sizeof(float)))==NULL) printf("I cannot dimension errflux  of %d elements \n",ndat);
  if((magni=malloc(ndat*sizeof(float)))==NULL) printf("I cannot dimension magni  of %d elements \n",ndat);
  if((xobj=malloc(ndat*sizeof(float)))==NULL) printf("I cannot dimension xobj  of %d elements \n",ndat);
  if((yobj=malloc(ndat*sizeof(float)))==NULL) printf("I cannot dimension yobj  of %d elements \n",ndat);
  if((ilog=malloc(ndat*sizeof(int)))==NULL) printf("I cannot dimension fluxcol  of %d elements \n",ndat);

  
  ReadNumcol(filecross,fluxcol,flux,ilog,&ndat);
  if(errfluxcol!=0)  ReadNumcol(filecross,errfluxcol,errflux,ilog,&ndat);
  ReadNumcol(filecross,magcol,magni,ilog,&ndat);
  ReadNumcol(filecross,xcol,xobj,ilog,&ndat);
  ReadNumcol(filecross,ycol,yobj,ilog,&ndat);


}







int  FluxCalibration(void) {
  double *logcts;
  double *errlogcts;
  double *m;
  double *resi;
  float ctsmin,ctsmax;
  float magmin,magmax;
  char title[200];
  int nfot;


  nfot=0;

  if((logcts=malloc(ndat*sizeof(double)))==NULL) printf("I cannot dimension logcts  of %d elements \n",ndat);
  if((errlogcts=malloc(ndat*sizeof(double)))==NULL) printf("I cannot dimension errlogcts  of %d elements \n",ndat);
  if((m=malloc(ndat*sizeof(double)))==NULL) printf("I cannot dimension m  of %d elements \n",ndat);
  if((resi=malloc(ndat*sizeof(double)))==NULL) printf("I cannot dimension resi  of %d elements \n",ndat);


  cpgsci(1);
  FluxCal_Basic(logcts,errlogcts,m,resi,&nfot);
  printf("m= (%f +- %f)  + (%f +- %f) * log10(Flux)",a,erra,b,errb);

  pgLimits_d(nfot,logcts,&ctsmin,&ctsmax);
  pgLimits_d(nfot,m,&magmin,&magmax);
  cpgswin(ctsmin,ctsmax,magmin,magmax);
/*    cpgswin(2.56,7.24,7.9,20.25); */
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);  
  cpgsch(1.);
  sprintf(title,"Calibration plot using file %s",filecross);
  cpglab("log10(Flux)","Calibration magnitude",title);

  if(fluxsatura==0) fluxsatura=pow(10.,ctsmax);


  PlotMagCal(logcts,m,errlogcts,nfot);
  cpgsci(3);
  FluxCal_Basic(logcts,errlogcts,m,resi,&nfot);
  printf("m= (%f +- %f)  + (%f +- %f) * log10(Flux)",a,erra,b,errb);
  cpgswin(ctsmin,ctsmax,magmin,magmax);
  PlotMagCal(logcts,m,errlogcts,nfot);
  cpgsci(4);
  FluxCal_Basic(logcts,errlogcts,m,resi,&nfot);
  printf("m= (%f +- %f)  + (%f +- %f) * log10(Flux)",a,erra,b,errb);
  cpgswin(ctsmin,ctsmax,magmin,magmax);
  PlotMagCal(logcts,m,errlogcts,nfot);
  cpgsci(2);
  FluxCal_Basic(logcts,errlogcts,m,resi,&nfot);
  printf("m= (%f +- %f)  + (%f +- %f) * log10(Flux)",a,erra,b,errb);
  cpgswin(ctsmin,ctsmax,magmin,magmax);
  PlotMagCal(logcts,m,errlogcts,nfot);
  cpgsci(1);
  sprintf(title,"m= (%f +- %f)  + (%f +- %f) * log10(Flux)",a,erra,b,errb);
  cpgtext(.1*(ctsmax-ctsmin)+ctsmin,.2*(magmax-magmin)+magmin,title);
  sprintf(title,"RMS fit: %f using %d points. Cov(a,b)=%g",sigma,nfot,covab);
  cpgtext(.1*(ctsmax-ctsmin)+ctsmin,.1*(magmax-magmin)+magmin,title);
  rms=sigma;

  if(nfot<2) {
    printf(" Not possible to fit a line in magnitude calibration\n");
    return(0);
  }
  free(logcts);free(errlogcts);free(m);free(resi);
  return(1);
}








void FluxCal_Basic(double *logcts,double *errlogcts,double *m,double *resi,int *nfot) {
  static int calflag=0;
  int i;
  double chi2;
  float q;
  double ap,bp,errap,errbp,covabp;
/*   float par[2]; */
/*   float covpar[4]; */
  double apww,bpww,errapww,errbpww,covabpww;

  static   double aresi,bresi;
  

  FILE *kk; 
  unlink("pajustar"); 
  kk=fopen("pajustar","w"); 
  
  if(calflag) StMedia_d(*nfot,resi,&sigma);
  printf(" sigma %f n %d\n",sigma,*nfot);

  *nfot=0;
  for(i=0;i<ndat;i++) {
    if(ilog[i]) {
/*       printf(" i %d fl %f sa %f m %f ssss %f n*s %f s %f\n",i,flux[i],fluxsatura,mag[i],fabs((double)(mag[i]-a-b*log10(flux[i]))),nsigrej*sigma,sigma);  */
      if( !calflag) {
        logcts[*nfot]=log10(flux[i]);
        if(errfluxcol==0) errlogcts[*nfot]=1;
        else              errlogcts[*nfot]=errflux[i]/flux[i]/log(10.);
        m[*nfot]=magni[i];
        (*nfot)++;
      }
      else if((flux[i]<fluxsatura  && fabs((double)(magni[i]-a-b*log10(flux[i])))<nsigrej*sigma)) {
	
	logcts[*nfot]=log10(flux[i]);
	if(errfluxcol==0) {
	  if(!calflag ) errlogcts[*nfot]=1;
/* 	  else         errlogcts[*nfot]=(fabs(logcts[*nfot]-(m[*nfot]-a)/b)); */
/* 	  else         errlogcts[*nfot]=(1/(logcts[*nfot])); */
	  else         errlogcts[*nfot]=aresi+bresi*logcts[*nfot];
          if(errlogcts[*nfot]<sigma) errlogcts[*nfot]=sigma;
	}
	else              errlogcts[*nfot]=errflux[i]/flux[i]/log(10.);
	m[*nfot]=magni[i];
	resi[*nfot]=fabs((float)(m[*nfot]-a-b*logcts[*nfot]));
        if(DEBUG2) printf(" Calculo resi %f\n",resi[*nfot]);
/*      	printf("  %f %f %f \n",m[*nfot],logcts[*nfot],errlogcts[*nfot]);       */
      	fprintf(kk," %f %f %f \n",m[*nfot],logcts[*nfot],errlogcts[*nfot]);   
	(*nfot)++;
      }
    }
  }
  if(DEBUG) printf(" Viy a hacer ajuste con %d puntos \n",*nfot);


/*   fclose(kk); */
  unlink("pajustar");
  printf(" NFOT %d\n",*nfot);

  if(calflag) {
    MCP1_d(*nfot,logcts,resi,&aresi,&bresi);
    printf(" RESI a %f  b %f \n",aresi,bresi); 
  }
  else {
    aresi=1;
    bresi=0;
  }
/*   cpgswin(3.,7.1,0.,2.); */
/*   PlotMagCal(logcts,resi,errlogcts,*nfot); */


  if(calflag) {
    MCP1_err_d(*nfot,m,logcts,&apww,&bpww,&errapww,&errbpww,&covabpww,&chi2,&q); 
    if(DEBUG)  printf(" a %f +- %f b %f +- %f cov %f\n",apww,errapww,bpww,errbpww,covabpww); 
    
    MCP1Weight_err_d(*nfot,m,logcts,errlogcts,&ap,&bp,&errap,&errbp,&covabp,&chi2,&q);  
  if(DEBUG)    printf(" a %f +- %f b %f +- %f cov %f\n",ap,errap,bp,errbp,covabp); 
/*     errap=errbp=covabp=0.0; */
  }
  else  MCP1_d(*nfot,m,logcts,&ap,&bp);
  

  /*    MCP1Weight_err(*nfot,logcts,m,errlogcts,&a,&b,&erra,&errb,&covab,&chi2,&q);     */
  /*    svdfit(m,logcts,errlogcts,*nfot,par,2,covpar,&chi2,&funlin); */

  if(errfluxcol==0) {
    errap=errapww;errbp=errbpww;covabp=covabpww;
  }
  else {
    errap=errapww;errbp=errbpww;covabp=covabpww;
  }
				  
  /*    printf(" Chi2 %f  Q %g\n",chi2,q); */
  /*    Ajusto al reves */
  /*    printf("REVESSS m= (%f +- %f)  + (%f +- %f) * log10(Flux)  Covar %g\n",ap,errap,bp,errbp,covabp);   */
  a=-ap/bp;  
  erra=sqrt( (errap/bp)*(errap/bp) + (ap/bp/bp*errbp)*(ap/bp/bp*errbp) -2*ap*covabp/(bp*bp*bp));  
  b=1./bp;  
  errb=errbp/bp/bp;  
  covab=-ap*errbp*errbp/bp/bp/bp/bp+covabp/bp/bp/bp;

  if(!calflag) {
    *nfot=0;
    for(i=0;i<ndat;i++) {
      if(ilog[i]) {
         resi[*nfot]=fabs((float)(m[*nfot]-a-b*logcts[*nfot]));
         if(DEBUG2) printf(" Calculo resi %f\n",resi[*nfot]);
         (*nfot)++;
      }  
    }
  }



  calflag=1;
  if (DEBUG) calflag=readi(calflag);
/*    exit(1); */
}






void PlotMagCal(double *logcts,double *m,double *errlogcts,int nfot) {
  int i;
  for(i=0;i<nfot;i++) {
     cpgpt1((double)logcts[i],(double)m[i],3); 
/*     if(errfluxcol!=0) cpgerr1(5,logcts[i],m[i],20*errlogcts[i],0.1); */
  }
  cpgmove(0.,a);
  cpgdraw(100.,a+b*100.);
/*    printf(" m = %f + %f *log10(flux)\n",a,b); */
}






int MagLim(void) {
  float *mcal;
  int nhist=0;
  int *histo,*cumhisto;
  float *x,*y,*sigy;
  double *x_d,*y_d,*sigy_d;
  int i;
  float mmin,mmax;
  int nhistmax = 0;
  int maxhist;
  float mhalfmaxd;
  int nhalfmaxd,halfmaxd;
  float mhalfmaxi;
  int nhalfmaxi,halfmaxi;
  float mthirdmax;
  int nthirdmax,thirdmax;
  float mtwothirdmax;
  int ntwothirdmax,twothirdmax;
  float par[2];
  int   ipar[2];
  float **covarpar;
  float chisq;
  double par_d[2];
  double **covarpar_d;
  double chisq_d;
  int iter;
/*    float *fsel; */
  float xdum;
  double sumf,sumy,fact;
  char title[200];

  covarpar=matrix_f(2,2);
  covarpar_d=matrix_d(2,2);
  mcal=vector_f(ndat);

  nhist=0;
  for(i=0;i<ndat;i++) {
    if(ilog[i]) {
      if(fabs((float)(magni[i]-a-b*log10(flux[i])))<nsigrej*sigma) {
	mcal[nhist]=a+b*log10(flux[i]);
  	if(DEBUG) printf(" f %f m %f %f x %f y %f\n",flux[i],magni[i],mcal[nhist],xobj[i],yobj[i]);  
	nhist++;
      }
    }
  }

  if(DEBUG) printf(" nhist %d\n",nhist);


  MinMax(nhist,mcal,&mminimum,&mmaximum);
  if(nhist<3 || mminimum==mmaximum ) {
    printf(" Not possible to fit selection function with this poor histogram\n");
    mmaximum=0;mmode=0;mfermicut=0;errmfermicut=0;deltafermi=0;errdeltafermi=0;gammapowerlaw=0;errgammapowerlaw=0;
    return 0;
  }


/*   mmaximum=mmax; */
/*   mminimum=mmin; */
  printf(" Maximum magnitude: %f\n",mmaximum);
  printf(" Minimum magnitude: %f\n",mminimum);
  mmax=mmaximum+(mmaximum-mminimum)*0.2;
  mmin=mminimum+0.6;

/*    mmax=27;  */
/*    mmin=15;  */


  histo=StHisto2(nhist,mcal,nbin,&mmin,&mmax);
  deltabin=(mmax-mmin)/nbin;
  maxhist=0;
  halfmaxd=0;
  for(i=0;i<nbin;i++ ){
    if(histo[i]>maxhist) {
      nhistmax=i;
      maxhist=histo[i];
    }
  }

  if(DEBUG) printf(" maxhist %d \n",maxhist);

  nbin=(int)(nbin*maxhist/100);
  if(nbin<20) nbin=20;
  if(nbin>200) nbin=200;
  printf(" NBIN %d\n",nbin);
  free(histo);
  histo=StHisto2(nhist,mcal,nbin,&mmin,&mmax);
  deltabin=(mmax-mmin)/nbin;
  maxhist=0;
  for(i=0;i<nbin;i++ ){
    if(histo[i]>maxhist) {
      nhistmax=i;
      maxhist=histo[i];
    }
  }
  for(i=nhistmax;i<nbin;i++) {
    if(histo[i]<maxhist/2.) break;
  }
  nhalfmaxd=i;
  halfmaxd=histo[i];
  if(halfmaxd==0) halfmaxd=1;
  for(i=0;i<nhistmax;i++) {
    if(histo[i]>maxhist/2.) break;
  }
  nhalfmaxi=i;
  halfmaxi=histo[i];
  for(i=0;i<nhistmax;i++) {
    if(histo[i]>maxhist/3.) break;
  }
  nthirdmax=i-1;
  thirdmax=histo[i-1];
  for(i=0;i<nhistmax;i++) {
    if(histo[i]>maxhist*2./3.) break;
  }
  ntwothirdmax=i+1;
  twothirdmax=histo[i+1];

  if(thirdmax==0) thirdmax=1;

  if(twothirdmax<thirdmax) {
    twothirdmax=maxhist;
    ntwothirdmax=nhistmax;
  }




/*   printf(" hi %d-%d his %f\n",nhalfmaxi,nhistmax,halfmaxi); */
  

  if(DEBUG) printf(" 22 maxhist %d \n",maxhist);


  cpgpage();
  cpgswin(mmin,mmax,0.,maxhist*1.2);
  cpgbox("BCTNS",0.0,0,"BTNS",0.0,0);  
  cpghist(nhist,mcal,mmin,mmax,nbin,1); 
  sprintf(title,"Selection function for image %s",imagefile);
  cpglab("Magnitude","Number of objects",title);

  mmode=mmin+(nhistmax+0.5)*(mmax-mmin)/nbin;
  mhalfmaxd=mmin+(nhalfmaxd+0.5)*(mmax-mmin)/nbin;
  mhalfmaxi=mmin+(nhalfmaxi+0.5)*(mmax-mmin)/nbin;
  mthirdmax=mmin+(nthirdmax+0.5)*(mmax-mmin)/nbin;
  mtwothirdmax=mmin+(ntwothirdmax+0.5)*(mmax-mmin)/nbin;
  printf(" Mode mag: %f \n",mmode);

  cumhisto=StCumHisto(nhist,mcal,nbin,&mmin,&mmax);
  cpgswin(mmin,mmax,0.,2.*(nhist));
  cpgbox("BCTNS",0.0,0,"CTMS",0.0,0);  
  cpgmove(mmin+0.5*(mmax-mmin)/nbin,(cumhisto[0]));
  for(i=1;i<nbin;i++ )  cpgdraw(mmin+(i+0.5)*(mmax-mmin)/nbin,(cumhisto[i]));

/*   printf(" %d %d %d %d %d %d\n",x,y,sigy,x_d,y_d,sigy_d); */
  
  if((x     =malloc(nbin*sizeof(float)))==NULL) printf("I cannot dimension x of %d elements \n",nbin);
  if((y     =malloc(nbin*sizeof(float)))==NULL) printf("I cannot dimension y of %d elements \n",nbin);
  if((sigy  =malloc(nbin*sizeof(float)))==NULL) printf("I cannot dimension sigy of %d elements \n",nbin);

  if((x_d   =malloc(nbin*sizeof(double)))==NULL) printf("I cannot dimension x_d of %d elements \n",nbin);
/*   printf(" %d %d %d %d %d %d\n",x,y,sigy,x_d,y_d,sigy_d); */
/*   printf(" IN KAKAK\n"); */
/*   free(x_d);   */
/*   printf(" IN 666 ka\n"); */
/*   printf(" %d %d %d %d %d %d\n",x,y,sigy,x_d,y_d,sigy_d); */

  if((y_d   =malloc(nbin*sizeof(double)))==NULL) printf("I cannot dimension y_d of %d elements \n",nbin);
  if((sigy_d=malloc(nbin*sizeof(double)))==NULL) printf("I cannot dimension sigy_d of %d elements \n",nbin);
/*   printf(" %d %d %d %d %d %d\n",x,y,sigy,x_d,y_d,sigy_d); */
/*   printf(" SA KAKAK\n"); */
/*   free(x_d);    */
/*   printf(" SA 666 ka\n"); */
  
/*   printf(" x_d %d x_df %d\n",x_d,x_d+nbin-1); */

  for(i=0;i<nbin;i++ ){
    x[i]=(i+0.5)*(mmax-mmin)/nbin+mmin;
    x_d[i]=(i+0.5)*(mmax-mmin)/nbin+mmin;
    y[i]=(float)(histo[i]);
    if(i!=nbin-1) y[i+1]=(float)(histo[i+1]);
    y_d[i]=(double)(histo[i]);
    if(i!=nbin-1) y_d[i+1]=(double)(histo[i+1]);

    if(i==0)               sigy[i]=sqrt((fabs(y[i])+fabs(y[i+1]))/2.); 
    if(i==nbin-1)          sigy[i]=sqrt((fabs(y[i-1])+fabs(y[i]))/2.); 
    if(i!=0 && i!=nbin-1)  sigy[i]=sqrt((fabs(y[i-1])+fabs(y[i])+fabs(y[i+1]))/3.); 
    if(sigy[i]<1)          sigy[i]=1;

    if(i==0)               sigy_d[i]=sqrt((fabs(y_d[i])+fabs(y_d[i+1]))/2.); 
    if(i==nbin-1)          sigy_d[i]=sqrt((fabs(y_d[i-1])+fabs(y_d[i]))/2.); 
    if(i!=0 && i!=nbin-1)  sigy_d[i]=sqrt((fabs(y_d[i-1])+fabs(y_d[i])+fabs(y_d[i+1]))/3.); 
    if(sigy_d[i]<1)        sigy_d[i]=1;
  }

/*   printf(" x_d %d\n",x_d); */
/*   for(i=0;i<nbin;i++ )         printf(" x %f y %f s %f\n",x[i],y[i],sigy[i]);  */

/*   printf(" 444 KAKAK\n"); */
/*   free(x_d);   */
/*   printf(" 444 666 ka\n"); */

  if(DEBUG) printf(" 33 maxhist %d \n",maxhist);


  cpgsci(2);

  ipar[0]=1;ipar[1]=1;
  printf(" %d %d %f %f\n",twothirdmax,thirdmax,mtwothirdmax,mthirdmax);
  printf(" %d %d %f %f\n",maxhist,halfmaxd,mmode,mhalfmaxd);
/*   par[0]=0.5; */
  par[0]=(mmode);
  if(DEBUG) printf(" mhm %f mmo %f ma %d hal %d\n",mhalfmaxd,mmode,maxhist,halfmaxd);  
  if(halfmaxd==maxhist) par[1]=(mhalfmaxd-mmode)/log(2);
  else                 par[1]=(mhalfmaxd-mmode)/log(2*(float)maxhist/(float)halfmaxd-1);
/*   printf(" TENTA %f %f %f\n",par[0],par[1],par[2]); */

  par_d[0]=par[0];par_d[1]=par[1];

  xglob=x;yglob=y; 
  xglob_d=x_d;yglob_d=y_d; 

  sumy=sumf=0.0;
  for(i=0;i<nbin;i++) {
    sumf+=Fermi(x[i],par[1],par[2]);
    sumy+=y[i]; 
  }  
  fact=sumy/sumf;
  cpgswin(mmin,mmax,0.,maxhist*1.2);
  
  cpgsci(3);
  for(i=0;i<nbin;i++ ){
/*     cpgdraw(x[i],DerPowLaw(x[i],fact,0.,par[0])*Fermi(x[i],par[1],par[2]));  */
  }
  cpgsci(5);
  par_d[0]=par[0];par_d[1]=par[1];;
  
  /*    exit(1); */
  if(DEBUG) printf(" 44 maxhist %d \n",maxhist);


  printf(" Calling MRQ...\n");
  if(DEBUG) {
    for(i=0;i<nbin;i++) printf(" x %f y %f s %f\n",x_d[i],y_d[i],sigy_d[i]);    
  }
      
  if(DEBUG)  printf(" Llamo a Mrq\n"); 
  iter=Mrq_d(x_d,y_d,sigy_d,nbin,par_d,ipar,2,covarpar_d,&chisq_d,funmaglim_d);
  if(DEBUG)  printf(" iter %d END %f %f  chi %f\n",iter,par_d[0],par_d[1],chisq_d);  
  
/*   goto salta;  */


  sumy=sumf=0.0;
  for(i=0;i<nbin;i++) {
    sumf+=Fermi(x_d[i],par_d[0],par_d[1]);
    sumy+=y_d[i]; 
  }  
  fact=sumy/sumf;
  if(DEBUG) printf(" gact %f sumy %f sumf %f\n",fact,sumy,sumf);
  mmax=mmaximum;
  if(fact==0) mmin=mminimum-(mmaximum-mminimum)*0.1;
  else     mmin=mminimum+0.6;

  if(DEBUG) printf(" par0 %f\n",par[0]);
  if(DEBUG) printf(" pard0 %f\n",par_d[0]);
  if(DEBUG) printf(" Aqui mmin %g maxhist %d fact %g\n",mmin,maxhist,fact);
  free(histo);
  histo=StHisto2(nhist,mcal,nbin,&mmin,&mmax);
  deltabin=(mmax-mmin)/nbin;
  maxhist=0;
  halfmaxd=0;
  for(i=0;i<nbin;i++ ){
    if(histo[i]>maxhist) {
      nhistmax=i;
      maxhist=histo[i];
    }
  }
  nbin=(int)(nbin*maxhist/100.);
  if(nbin<20) nbin=20;
  if(nbin>200) nbin=200;
  printf(" NBIN %d\n",nbin);
  printf(" mmax %f mmin %f\n",mmax,mmin);
  free(histo);
  if(DEBUG) {
    for(i=0;i<nhist;i++) printf(" mcal %f\n",mcal[i]);  
  }

  histo=StHisto2(nhist,mcal,nbin,&mmin,&mmax);
  deltabin=(mmax-mmin)/nbin;
  maxhist=0;
  for(i=0;i<nbin;i++ ){
    if(histo[i]>maxhist) {
      nhistmax=i;
      maxhist=histo[i];
    }
  }
  cpgpage();
  cpgsci(1);
  cpgswin(mmin,mmax,0.,maxhist*1.2);
  cpgbox("BCTNS",0.0,0,"BTNS",0.0,0);  
  cpghist(nhist,mcal,mmin,mmax,nbin,1); 
  cpglab("Magnitude","Number of objects",title);
  mmode=mmin+(nhistmax+0.5)*(mmax-mmin)/nbin;
  printf(" Mode mag: %f \n",mmode);
  free(cumhisto);
  if(DEBUG) printf(" Despues\n");
  cumhisto=StCumHisto(nhist,mcal,nbin,&mmin,&mmax);
  cpgswin(mmin,mmax,0.,2.*(nhist));
  cpgbox("BCTNS",0.0,0,"CTMS",0.0,0);  
  cpgmove(mmin+0.5*(mmax-mmin)/nbin,(cumhisto[0]));
  for(i=1;i<nbin;i++ ) cpgdraw(mmin+(i+0.5)*(mmax-mmin)/nbin,(cumhisto[i]));
  xglob_d=NULL;
  yglob_d=NULL;
/*   printf(" Antes free\n");  */
  free(x_d);free(y_d);free(sigy_d); 
/*   printf(" DEspeues fre\n");  */
  if((x_d   =malloc(nbin*sizeof(double)))==NULL) printf("I cannot dimension x_d of %d elements \n",nbin);
  if((y_d   =malloc(nbin*sizeof(double)))==NULL) printf("I cannot dimension y_d of %d elements \n",nbin);
  if((sigy_d=malloc(nbin*sizeof(double)))==NULL) printf("I cannot dimension sigy_d of %d elements \n",nbin);
/*    printf(" DEsp allo\n");  */
  xglob_d=x_d;yglob_d=y_d; 
  for(i=0;i<nbin;i++ ) {
    x_d[i]=(i+0.5)*(mmax-mmin)/nbin+mmin;
    y_d[i]=(double)(histo[i]);
    if(i!=nbin-1) y_d[i+1]=(double)(histo[i+1]);
    if(i==0)               sigy_d[i]=sqrt((fabs(y_d[i])+fabs(y_d[i+1]))/2.); 
    if(i==nbin-1)          sigy_d[i]=sqrt((fabs(y_d[i-1])+fabs(y_d[i]))/2.); 
    if(i!=0 && i!=nbin-1)  sigy_d[i]=sqrt((fabs(y_d[i-1])+fabs(y_d[i])+fabs(y_d[i+1]))/3.); 
    if(sigy_d[i]<1)        sigy_d[i]=1;
  }

  printf(" Calling MRQ...\n");
/*   printf(" iter INI %f %f %f \n",par_d[0],par_d[1],par_d[2]);  */
  cpgswin(mmin,mmax,0.,maxhist*1.2);
  if(DEBUG) {
    for(i=0;i<nbin;i++) printf(" x %f y %f s %f\n",x_d[i],y_d[i],sigy_d[i]);    
  }

  iter=Mrq_d(x_d,y_d,sigy_d,nbin,par_d,ipar,2,covarpar_d,&chisq_d,funmaglim_d);
/*   printf(" iter %d END %f %f %f chi %f\n",iter,par_d[0],par_d[1],par_d[2],chisq_d);  */
 
  


/*  salta: */




  par[0]=par_d[0];
  par[1]=par_d[1];
  covarpar[0][0]=covarpar_d[0][0];
  covarpar[1][0]=covarpar_d[1][0];
  covarpar[0][1]=covarpar_d[0][1];
  covarpar[1][1]=covarpar_d[1][1];
  chisq=chisq_d;
  printf(" DONE\n");
  
  mfermicut=par_d[0];deltafermi=par_d[1];gammapowerlaw=1;
  errmfermicut=sqrt(covarpar_d[0][0]);errdeltafermi=sqrt(covarpar_d[1][1]);errgammapowerlaw=0;


  printf(" Fermi magnitude  : %g +/- %g\n",mfermicut,errmfermicut);
  printf(" Fermi temperature: %g +/- %g\n",deltafermi,errdeltafermi);
  printf(" Exponential power: %g +/- %g\n\n",gammapowerlaw,errgammapowerlaw);



  cpgsci(6);
  cpgswin(mmin,mmax,0.,1.2);
  cpgmove(mmin,1.);
  for(i=0;i<nbin;i++ ){
    xdum=i*(mmax-mmin)/nbin+mmin;
    cpgdraw(xdum,Fermi(xdum,(double)par_d[0],(double)par_d[1]));
  }
  sumy=sumf=0.0;
  for(i=0;i<nbin;i++) {
    sumf+=Fermi(x_d[i],par_d[0],par_d[1]);
    sumy+=y_d[i]; 
  }  
  fact=sumy/sumf;
  
/*   printf("  fact %g  sumy %g sumf %g \n",fact,sumy,sumf); */

  cpgsci(2);
  cpgswin(mmin,mmax,0.,maxhist*1.2);
  cpgmove(mmin,0.);
  for(i=0;i<nbin;i++ ){
    xdum=i*(mmax-mmin)/nbin+mmin;
    cpgdraw(xdum,fact);
  }
  
  cpgsci(3);
  cpgmove(mmin,0.);
  for(i=0;i<nbin;i++ ){
    xdum=i*(mmax-mmin)/nbin+mmin;
    cpgdraw(xdum,fact*Fermi(xdum,par_d[0],par_d[1]));
  }

  free(x);free(y);free(sigy); 
  free(x_d);free(y_d);free(sigy_d); 
  free(histo);free(cumhisto);
  free(mcal);

  free_matrix_f(covarpar,2,2);
  free_matrix_d(covarpar_d,2,2);

  return 1;

}













int  Seeing(void) {
  fitsfile *fitsimage;
  long naxes[2];
  int nx;
  float *ima;
  int nfound,anynull;
  long fpixel,npixels;
  float nullval;
  int status=0;
  int ipar[4];
  float par[4],sigpar[4];
  float chisq;
  float **covarpar;

/*    float *xmon,*ymon; */
/*    float p0[NMON],p1[NMON],p2[NMON],p3[NMON]; */
/*    float fnul; */

  float *cutx,*cuty,*xpaint,*yfit;
  float *sigcuty;
  int nsum=6.;  /*Half of the pixels to add */
  
  int xi,yi;
  int i,j,k;
  float dmin,dmax;
  int iter;

  float *fitsigma, *fitpeak, *fitxpos, *weight;
  int nseeing=0;

  float minsigma,maxsigma,minpeak,maxpeak,meanpeak;
  float modeseeing,meanseeing,medianseeing;
  float tr[6];

  covarpar=matrix_f(4,4);

  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  
  if( ffopen(&fitsimage, imagefile, READONLY, &status)) fits_report_error(stderr,status);
  if(fits_read_keys_lng(fitsimage, "NAXIS", 1, 2, naxes, &nfound, &status)) fits_report_error(stderr,status);
  nx=naxes[0];

  npixels=naxes[0]*naxes[1];
  if((ima=malloc(npixels*sizeof(float)))==NULL) printf("I cannot dimension ima of %ld elements \n",npixels);
  fpixel=1;
  nullval=0;
  if(fits_read_img(fitsimage, TFLOAT, fpixel, npixels, &nullval, ima, &anynull, &status )) fits_report_error(stderr,status);
  
  nspacut=(int)(25./platescale/2.);
/*   cpgwnad(0.,naxes[0],0.,naxes[1]); */
/*   cpgwnad(834.-nspacut,834.+nspacut,111.-nspacut,111.+nspacut); */
/*   cpggray(ima,naxes[0],naxes[1],834-nspacut,834+nspacut,111-nspacut,111+nspacut,15000,5000,tr); */
/*   cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0); */

    

  if((cutx   =malloc(nspacut*2*sizeof(float)))==NULL) printf("I cannot dimension cutx of %d elements \n",nspacut*2);
  if((cuty   =malloc(nspacut*2*sizeof(float)))==NULL) printf("I cannot dimension cuty of %d elements \n",nspacut*2);
  if((xpaint =malloc(nspacut*2*sizeof(float)))==NULL) printf("I cannot dimension xpaint of %d elements \n",nspacut*2);
  if((yfit   =malloc(nspacut*2*sizeof(float)))==NULL) printf("I cannot dimension yfit of %d elements \n",nspacut*2);
  if((sigcuty=malloc(nspacut*2*sizeof(float)))==NULL) printf("I cannot dimension sigcuty of %d elements \n",nspacut*2);

  
  if((fitsigma=malloc(ndat*sizeof(float)))==NULL) printf("I cannot dimension fitsigma of %d elements \n",ndat);
  if((fitpeak =malloc(ndat*sizeof(float)))==NULL) printf("I cannot dimension fitpeak of %d elements \n",ndat);
  if((fitxpos =malloc(ndat*sizeof(float)))==NULL) printf("I cannot dimension fitxpos of %d elements \n",ndat);
  if((weight  =malloc(ndat*sizeof(float)))==NULL) printf("I cannot dimension weight of %d elements \n",ndat);



/*    if((xmon=malloc(nspacut*2*sizeof(float)))==NULL) printf("I cannot dimension spacut of %d elements \n",nspacut*2); */
/*    if((ymon=malloc(nspacut*2*sizeof(float)))==NULL) printf("I cannot dimension spacut of %d elements \n",nspacut*2); */


  
  ipar[0]=1;ipar[1]=1;ipar[2]=1;ipar[3]=1;
  for(i=0;i<ndat;i++)  {
    if(ilog[i]) {
/*       printf(" %d %d %d %d %d %d\n",xobj[i]>nspacut*2,yobj[i]>nspacut*2,xobj[i]<naxes[0]-1-nspacut*2,yobj[i]<naxes[0]-1-nspacut*2,flux[i]> pow(10.,(mfermicut-a)/b),flux[i]<pow(10.,(mminimum+2-a)/b) );  */

      if(xobj[i]>nspacut*2 && yobj[i]>nspacut*2 && xobj[i]<naxes[0]-1-nspacut*2 && yobj[i]<naxes[0]-1-nspacut*2 && flux[i]> pow(10.,(mfermicut-a)/b) && flux[i]<pow(10.,(mminimum+2-a)/b) ) { 
/*        if(xobj[i]>nspacut*2+20 && yobj[i]>nspacut*2+20 && xobj[i]<naxes[0]-21-nspacut*2 && yobj[i]<naxes[0]-21-nspacut*2 ) { */

/* 	printf(" 666 xobj %f yobj %f\n",xobj[i],yobj[i]); */

	for(j=0;j<nspacut*2;j++) {
	  cutx[j]=0.;
	  cuty[j]=0.;
	  xpaint[j]=(float)j;
	  /* 	printf("Inic\n"); */
	  for(k=0;k<nsum;k++) {
	    xi=(int)(xobj[i]-nsum+k+1);
	    yi=(int)(yobj[i]-nspacut+j+1);
	    if(xi>0 && xi<naxes[0]-1 && yi>0 && yi<naxes[1]-1)  cuty[j]+=ima[xi+yi*nx];
	    xi=(int)(xobj[i]-nspacut+j+1);
	    yi=(int)(yobj[i]-nsum+k+1);
	    if(xi>0 && xi<naxes[0]-1 && yi>0 && yi<naxes[1]-1)   cutx[j]+=ima[xi+yi*nx];
	  }
	  sigcuty[j]=sqrt(fabs(cuty[j]));
	  if(sigcuty[j]==0) sigcuty[j]=1;
/*  	  if(i>143) printf(" xp %f y %f s %f\n",xpaint[j],cuty[j],sigcuty[j]);  */
	}
	MinMax(nspacut*2,cuty,&dmin,&dmax);

	if(dmin==0 && dmax==0) continue;

	xpglob=xpaint;cutyglob=cuty;
/* 	dmin=55000; */
/* 	dmax=68000; */
/* 	cuty[0]=55000;sigcuty[0]=sqrt(cuty[0]); */
/* 	cuty[1]=55000;sigcuty[1]=sqrt(cuty[1]); */
/* 	cuty[2]=55000;sigcuty[2]=sqrt(cuty[2]); */
  	par[0]=dmin;par[1]=dmax-dmin;par[2]=xpaint[nspacut];par[3]=2./platescale;  
	
/* 	pgLimits(nspacut*2,cuty,&dmin,&dmax); */
/* 	cpgpage(); */
/* 	cpgsci(2); */
/* 	cpgswin(0.,nspacut*2+1,dmin,dmax); */
/* 	cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0); */
/* 	cpgline(nspacut*2,xpaint,cuty); */
/* 	cpgsci(3); */


	
/*  	for(l=0;l<NMON;l++) { */
/*  	  for(j=0;j<nspacut*2;j++) { */
/*  	    ymon[j]=cuty[j]+sigcuty[j]*Gasdev(); */
/*  	  } */
/*  	  iter=Mrq(xpaint,ymon,sigcuty,nspacut*2,par,ipar,4,covarpar,&chisq,funseeing);	 */
/*  	  printf("\nS %d Gfit %3d iter MRQ solpar: ",l,iter);   */
/*  	  for(j=0;j<4;j++) printf("  par[%d] %8g",j,par[j]);  */
/*  	  p0[l]=par[0];p1[l]=par[1];p2[l]=par[2];p3[l]=par[3]; */
/*  	} */
/*  	printf("Errores MONTE \n"); */
/*  	StMedia(NMON,p0,&fnul); */
/*    	printf(" errp0 %g ",fnul);  */
/*  	StMedia(NMON,p1,&fnul); */
/*    	printf(" errp0 %g ",fnul);  */
/*  	StMedia(NMON,p2,&fnul); */
/*    	printf(" errp0 %g ",fnul);  */
/*  	StMedia(NMON,p3,&fnul); */
/*    	printf(" errp0 %g ",fnul);  */

/* 	if(i!=339 && i!=340 && !(i>340 && i<700)) { */

 	if(DEBUG) printf("%d  iter INI %f %f %f %f\n",i,par[0],par[1],par[2],par[3]);  
	iter=Mrq(xpaint,cuty,sigcuty,nspacut*2,par,ipar,4,covarpar,&chisq,funseeing);
	if(DEBUG) printf("%d-%d  iter %d END %f %f %f %f chi %f\n",i,ndat,iter,par[0],par[1],par[2],par[3],chisq);  
	sigpar[0]=sqrt(covarpar[0][0])    ;sigpar[1]=sqrt(covarpar[1][1]);
	sigpar[2]=sqrt(covarpar[2][2]);sigpar[3]=sqrt(covarpar[3][3]);

/* 	} */
/* 	else iter=0; */

/*      	printf("\nGfit %3d iter MRQ solpar: ",iter);  */
/*     	for(j=0;j<4;j++) printf("  par[%d] %8g",j,par[j]);   */
/*   	printf(" chi2 %f\n",chisq);  */
/*      	printf("                        : ",iter);    */
/*     	for(j=0;j<4;j++) printf(" epar[%d] %8g",j,sigpar[j]);  */
/*   	printf("\n");  */
 
/*  	for(j=0;j<nspacut*2;j++) yfit[j]=par[0]+par[1]*gaussian(xpaint[j],par[2],par[3]); */

/*  	for(j=0;j<4;j++) { */
/*  	  for(l=0;l<4;l++) { */
/*  	    printf("Covar %d %d teo %f\n",l,j,covarpar[j+l*4]); */
/*  	  } */
/*  	} */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p0,p0,&fnul,&fnul)); */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p0,p1,&fnul,&fnul)); */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p0,p2,&fnul,&fnul)); */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p0,p3,&fnul,&fnul)); */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p1,p0,&fnul,&fnul)); */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p1,p1,&fnul,&fnul)); */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p1,p2,&fnul,&fnul)); */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p1,p3,&fnul,&fnul)); */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p2,p0,&fnul,&fnul)); */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p2,p1,&fnul,&fnul)); */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p2,p2,&fnul,&fnul)); */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p2,p3,&fnul,&fnul)); */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p3,p0,&fnul,&fnul)); */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p3,p1,&fnul,&fnul)); */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p3,p2,&fnul,&fnul)); */
/*  	printf(" Covar MON %f\n",StCovar(NMON,p3,p3,&fnul,&fnul)); */



/*  	cpgsci(1);  */
/*  	cpgswin(0.,nspacut*2.,dmin,dmax); */
/*  	cpgbox("BCTNS",0,0,"BCTNS",0,0);  */
/*  	cpgline(nspacut*2,xpaint,cuty);  */
/*  	cpgsci(4); */
/*  	cpgline(nspacut*2,xpaint,yfit);  */

/*  	amoegauss=1; */
/*  	par[0]=dmin;           sigpar[0]=sqrt(fabs(dmin)); */
/*  	par[1]=dmax-dmin;      sigpar[1]=sqrt(dmax); */
/*  	par[2]=xpaint[nspacut];sigpar[2]=1.; */
/*  	par[3]=2./platescale;  sigpar[3]=par[3]/2.; */
/*  	iter=Amoeba(nspacut*2,xpaint,cuty,4,par,sigpar,(FTOL)*10.,500); */
/*     	printf("\nGfit %3d iter AMO solpar: ",iter);   */
/*  	for(j=0;j<4;j++) printf("  par[%d] %8g",j,par[j]); */
/*  	printf(" chi2 %f\n",Amoe_Funk(nspacut*2,xpaint,cuty,par)); */
/*  	for(j=0;j<nspacut*2;j++) yfit[j]=par[0]+par[1]*gaussian(xpaint[j],par[2],par[3]); */
/*  	cpgsci(2); */
/*  	cpgline(nspacut*2,xpaint,yfit);  */



/*   	 printf("\n Objeto %d x %f y %f f %f\n",i, xobj[i],yobj[i],flux[i]); 	  */
	 
/*  	cpgpage();  */


	if(iter!=0 && (par[3]*platescale)<10. && par[3]>0 && par[1]>0 && (sqrt(par[1]/pow(10.,(mminimum+1-a)/b)))<1.) {
	  fitsigma[nseeing]=par[3]*platescale;
	  fitpeak[nseeing]=par[1];
	  fitxpos[nseeing]=par[2];
	  weight[nseeing]=sqrt(fitpeak[nseeing]/pow(10.,(mminimum+1-a)/b));
/*     	  printf("n %d-%d  we %f ss %f\n",nseeing,ndat,weight[nseeing],fitsigma[nseeing]);    */
	  nseeing++;
	}
/* 	pgLimits(nspacut*2,cuty,&dmin,&dmax); */
/* 	cpgpage(); */
/* 	cpgsci(2); */
/* 	cpgswin(0.,nspacut*2+1,dmin,dmax); */
/* 	cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0); */
/* 	cpgline(nspacut*2,xpaint,cuty); */
/* 	cpgsci(3); */
/* 	cpgmove(0.,0.);  */
/* 	for(k=0;k<nspacut*2;k++) cpgdraw((float)k,par[0]+par[1]*gaussian((float)k,par[2],par[3]));  */

/* 	k=readi(k); */


      }
    }
  }
  if(nseeing<3) {
    printf(" Not possible to obtain mean seeing with less than 3 points\n");
    seeing=0;
    stdseeing=0;
    return(0);
  }

/*    for(i=0;i<nseeing;i++) printf("%d  %f s %f w %f\n",i,fitpeak[i],fitsigma[i],weight[i]); */
  meanseeing=StWeightMedia(nseeing,fitsigma,fitpeak,&stdseeing);
  Quartil(nseeing,fitsigma,&minsigma,&medianseeing,&maxsigma);
  printf(" 1 %f 2 %f 3 %f\n",minsigma,medianseeing,maxsigma);
  minsigma=minsigma-5*(medianseeing-minsigma);
  maxsigma=maxsigma+5*(maxsigma-medianseeing);
  
  modeseeing=StModa(nseeing,fitsigma,25,&minsigma,&maxsigma);
/*   printf(" LA moda es  limite %f %f\n",minsigma,maxsigma); */


  printf(" Mean seeing: %f with std: %f\n",meanseeing,stdseeing);
  printf(" Mode of the seeing: %f\n",modeseeing);
  StMedia(nseeing,fitsigma,&stdseeing);
  printf(" Mean seeing without weight: %f with std: %f\n",StMedia(nseeing,fitsigma,&stdseeing),stdseeing);
  seeing=StWeightMedia(nseeing,fitsigma,fitpeak,&stdseeing);
/*   pgLimits(nseeing,fitsigma,&minsigma,&maxsigma); */
  pgLimits(nseeing,fitpeak,&minpeak,&maxpeak);
  Quartil(nseeing,fitpeak,&minpeak,&meanpeak,&maxpeak);
  minpeak=meanpeak-3*(meanpeak-minpeak);
  maxpeak=meanpeak+3*(maxpeak-meanpeak);
/*   printf(" Limties %f %f %f %f\n",meanpeak,maxpeak,minsigma,maxsigma); */

  
  do  {
    cpgsci(1);
    cpgpage();
    cpgswin(minpeak,maxpeak,minsigma,maxsigma);
    cpgbox("BCTNS",0,0,"BCTNS",0,0);
    cpgpt(nseeing,fitpeak,fitsigma,2);
    cpgsci(2);
    cpgmove(minpeak,meanseeing);
    cpgdraw(maxpeak,meanseeing);
    cpgsci(3);
    cpgmove(minpeak,modeseeing);
    cpgdraw(maxpeak,modeseeing);
/*     cpgcurs(&minpeak,&minsigma,&cnul); */
/*     cpgband(2,1,minpeak,minsigma,&maxpeak,&maxsigma,&cnul); */
  } while(minsigma>maxsigma); 

  seeing=modeseeing*platescale*2.35;
  stdseeing=stdseeing*platescale*2.35;

  free_matrix_f(covarpar,3,3);

  return 1;
}

void PrintLog(void) {

  FILE *fr;
  int filexist=1;
  char wcsar[16],wcsdec[16];

/*   printf(" Esta en el prl\n"); */
  if((fr=fopen(filereg,"r")) ==NULL) {
    filexist=0;
/*     printf(" Mas mejor\n"); */
    if((fr=fopen(filereg,"w")) == NULL) {
      printf(" Couldn't create file %s. Exiting\n",filereg);
      exit(1);
    }
/*     printf("VREADC\n"); */
    fprintf(fr,"#%-51s %-16s  %-16s  %-9s  %-9s  %-8s  %10s  %10s  %10s  %-51s  %-10s  %-9s  %-9s  %-9s  %-9s  %-9s  %-6s  %-7s  %-7s  %-9s  %-9s  %-9s  %-9s  %-9s  %-9s  %-7s  %-7s  %-51s  %-51s\n","Image","Center_RA","Center_DEC","RA_size","DEC_size","Rot(ang)","Epoch","Std_X","Std_Y","Calibration_file","Num_obj","Mag_cal_A","Err_Mag_A","Mag_cal_B","Err_Mag_B","Cov_AB","RMS_MC","max_mag","mode_ma","m_fermi","errm_fermi","del_fer","errdel_fer","Alpha_PL","errAlpha_PL","Seeing","Std_See","Spectra_file","Response_file");
  }
  else {
    fr=fopen(filereg,"a");
    filexist=1;
  }
  ra2str(wcsar,16,alfac,3);
  dec2str(wcsdec,16,deltac,3);
  if(asstatus)
    fprintf(fr,"%-51s  %-16s  %-16s  %9.2g  %9.2g  %8.3f  %10.3f  %10.5f  %10.5f  ",imagefile,wcsar,wcsdec,xdim,ydim,rot,epoch,aststdx,aststdy);
  else 
    fprintf(fr,"%-51s  %-16s  %-16s  %-9s  %-9s  %-8s  %-10s  %10s  %10s  ",imagefile,"INDEF","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF");
  if(mcstatus)
    fprintf(fr,"%-51s  %10d  %9.4f  %9.4g  %9.4f  %9.4g  %9.4g  %6.3f  ",filecross,ndat,a,erra,b,errb,covab,rms);
  else 
    fprintf(fr,"%-51s  %10d  %-9s  %-9s  %-9s  %-9s  %-9s  %-6s  ",filecross,ndat,"INDEF","INDEF","INDEF","INDEF","INDEF","INDEF");
  if(mlstatus)
    fprintf(fr,"%7.3f  %7.3f  %9.4f  %9.4f  %9.4g  %9.4f  %9.4f  %9.4g  ",mmaximum,mmode,mfermicut,errmfermicut,deltafermi,errdeltafermi,gammapowerlaw,errgammapowerlaw);
  else 
    fprintf(fr,"%-7s  %-7s  %-9s  %-9s  %-9s  %-9s  %-9s  %-9s  ","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF","INDEF");
  if(sestatus)
    fprintf(fr,"%7.3f  %7.3f  ",seeing,stdseeing);
  else 
    fprintf(fr,"%-7s  %-7s  ","INDEF","INDEF");
  if(specstatus)
    fprintf(fr,"%-51s  %-51s\n",specfile,respfile);
  else 
    fprintf(fr,"%-51s  %-51s\n","INDEF","INDEF");

  fclose(fr);

}


float logFermi(float x,float mu,float T) {
  return(log(1./(exp((x-mu)/T)+1.)));
}

float PowLaw(float x,float A,float x0,float gam)
{
  return(A*pow((x-x0),gam));
}
float DerPowLaw(float x,float A,float x0,float gam)
{
  return(A*exp(x*gam+x0)); 
/*    return(gamma*A*pow((x-x0),gamma-1.));  */
/*    return(A*pow((x/x0),gamma-1.)); */
}

float logDerPowLaw(float x,float A,float x0,float gam)
{
  return(log(gam)+log(A)+(gam-1.)*log(x-x0));
}
float Fselerf(float x,float A,float B) {
  return(0.5*(1-erf(-pow(10.,-0.4*(x-A))+B)));
}

float Amoe_Funk(int n, float *x, float *y, float *p)
{
  int i;
  double s=0; 
  float fu,logfu;
  float media,sigmax;
  float sumf=0.;
  float sumy=0.;
  float fact;
  media=StMedia(n,x,&sigmax);
/*    for(i=0;i<3;i++) printf(" %f ",p[i]);    */
/*    printf("\n");    */
  
  if(amoegauss) {
    if(p[1]<0) return(1e29);
    if(p[3]<0) return(1e29);
    for(i=0;i<n;i++) sumf+=p[1]*gaussian(x[i],p[2],p[3]);
    for(i=0;i<n;i++) sumy+=y[i]-p[0];
    p[1]=p[1]*sumy/sumf;
    for(i=0; i<n; i++) {
      fu=p[0]+p[1]*gaussian(x[i],p[2],p[3]);
/*       s -=exp(-(fu-y[i])*(fu-y[i])/2./y[i]/y[i]); */
      s +=(fu-y[i])*(fu-y[i])/y[i];
/*       printf(" fu %f y %f s %g\n",fu,y[i],s); */
    }
  }
  else     {
    for(i=0;i<n;i++) {
      sumf+=DerPowLaw(x[i],1.,0.,p[0])*Fermi(x[i],p[1],p[2]);
    }  
    for(i=0;i<n;i++) sumy+=y[i]; 
    fact=sumy/sumf;
/*      printf(" p0 nor %f\n",p[0]);  */

    if(p[2]<=0.) return(1.e29);  
/*      if(p[1]<=1.) return(0.);  */
/*      if(p[2]>=media) return(0.);  */
    /*      if(p[0]<=0.) return(1.e29); */
    /*      if(p[1]<=1.) return(1.e29); */
    /*      if(p[2]>=media) return(1.e29); */
    for(i=0; i<n; i++) {
      if(y[i]>0.) {
	logfu=logDerPowLaw(x[i],p[0],p[2],p[1])+logFermi(x[i],p[3],p[4]);
	fu=DerPowLaw(x[i],fact,0.,p[0])*Fermi(x[i],p[1],p[2]);
	/*    	fu=DerPowLaw(x[i],p[0],p[2],p[1])*Fselerf(x[i],p[3],p[4]);  */
	/*        if(fu<1.) fu=1.;  */
	/*        s -= (log(fu)-log(y[i]))*(log(fu)-log(y[i])); */
	/*        s -= exp((1.-fu/y[i])*(1.-fu/y[i])/(fu/y[i])); */
	/*        s -= exp(-y[i]*(log(fu)-log(y[i]))*(log(fu)-log(y[i])));   */
  	s -= exp(-(fu-y[i])*(fu-y[i])/2./y[i]);    
/*  	s += ((fu-y[i])*(fu-y[i])/y[i]);    */
	/*     	  printf(" f %f logf %f pl %f  fe %f x %f y %f  s  %f\n",fu,logfu,DerPowLaw(x[i],p[0],p[2],p[1]),Fermi(x[i],p[3],p[4]),x[i],y[i],s);   */
	/*       	  printf(" ");   */
      }
/*        else s += ((fu-y[i])*(fu-y[i])); */
    }
  }
/*    printf("   chi      %f   \n",s);   */
  return(s);
}


void funlin(float x,float *y,int n) {
  y[0]=1;
  y[1]=x;
}



void funmaglim_d(double x,double *p,double *y,double *dyda,int n) {
  int i;
  double sumf=0.0,sumy=0.0;
  double fact;
  double xmin,xmax;

  if(DEBUG)  printf(" Antes p1 %f p2 %f delta %f\n",p[0],p[1],deltabin); 

  if(p[0]<xglob_d[0]) {
    p[0]=mmode; 
  }
  if(p[0]>xglob_d[nbin-1]) {
    p[0]=mmode; 
  }
   if(p[1]<deltabin/5.) { 
     p[1]=2*deltabin;  
   } 
  if(p[1]<0) {
    p[1]=0.1; 
  }


  MinMax_d(nbin,xglob_d,&xmin,&xmax);

  sumf=0;
  sumy=0;
  for(i=0;i<nbin;i++) {
    sumf+=Fermi(xglob_d[i],p[0],p[1]);
    sumy+=yglob_d[i]; 
  }
  
  fact=sumy/sumf;

  if(DEBUG)  printf(" fact %f sumy %f sumf %f \n",fact,sumy,sumf);


  *y=fact*Fermi(x,p[0],p[1]); 
  if(DEBUG) printf(" f %g e %g\n",Fermi(x,p[0],p[1]),(exp((x-p[0])/p[1])));
  dyda[0]=*y*Fermi(x,p[0],p[1])*(exp((x-p[0])/p[1]))/p[1];
  dyda[1]=+dyda[1]/p[1]*(x-p[0]);
  if(DEBUG) printf("Despues  y %f p1 %f p2 %f  \n",*y,p[0],p[1]); 


}


void funseeing(float x,float *p,float *y,float *dyda,int n) {

  float fac,ex,arg;
  double sumf=0.0,sumy=0.0;
  int i;
  if(DEBUG2) printf(" antes p0 %f p1 %f p2 %f p3 %f\n",p[0],p[1],p[2],p[3]);    

  if(p[1]<0)  p[1]=-p[1];
  if(p[3]<0)  p[3]=-p[3];
  if(p[3]>2*nspacut) p[3]=2./platescale;   /* No puede ser una gaussiana muy ancha */
  if(p[3]<0.1) p[3]=1;                     /* Ni mucho menos estrecha que un pixel */
  if(p[2]+2*p[3]<0) p[2]=nspacut;
  if(p[2]-2*p[3]>nspacut*2) p[2]=nspacut;

  for(i=0;i<nspacut*2;i++) {
    sumf+=p[1]*gaussian(xpglob[i],p[2],p[3]);
    sumy+=cutyglob[i]-p[0];
    if(DEBUG2) printf(" x %f y %f sumf %f sumy %f\n",xpglob[i],cutyglob[i],sumf,sumy); 
  } 
  if(p[1]==0) p[1]=sumy;
  else        p[1]=p[1]*sumy/sumf;
  
/*   printf(" norma %g %g\n",sumy,sumf); */

/*   cpgmove(0.,0.);  */
/*   for(i=0;i<nspacut*2;i++) cpgdraw((float)i,p[0]+p[1]*gaussian((float)i,p[2],p[3]));  */
  
  
/*   printf(" p0 %f p1 %f p2 %f p3 %f\n",p[0],p[1],p[2],p[3]);    */
  arg=(x-p[2])/p[3];
  ex=gaussian(x,p[2],p[3]);
/*    ex=exp(-arg*arg); */
  fac=p[1]*2.0*arg*ex;
  *y=p[0]+p[1]*ex;
  dyda[0]=1;
  dyda[1]=ex;
  dyda[2]=fac/p[3];
  dyda[3]=fac*arg/p[3];
  
/*   printf(" Sale funsein\n");  */

}

int  SpecInfo(void) {

  return(1);

}
