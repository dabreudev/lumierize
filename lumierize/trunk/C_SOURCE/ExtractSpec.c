#include "modulos.h"
#define NMAX 100000

/* Variables extraccion de espectros */
char rawspec[51];
char wcspec[51];
char fcspec[51];
int dispaxis;
char errrawspec[51];
char errwcspec[51];
char errfcspec[51];
int spapix,specpix;
int isradec;
int waveflag=0;
int fluxflag=0;
int ismultspec=0;

/* char objbox[51]; */

/* Variables fichero BUHO  */
char catfile[51];
int xcol,ycol;

/* Variables para las imagenes de fondo de cielo */

fitsfile *skyfits,*errskyfits;
char fitssky[51];
char fitserrsky[51];

/* int xdircol,ydircol; */
/* FILE *skyfile,*errskyfile; */

/* int skycol,sigcol; */
/* int errskycol; */


/* Variables para la extraccion de catálogo astrométrico */
double xpos,ypos;
int off;


/* Variables para la transformacion astrometrica */
struct WorldCoor *wcsim=0;      /* World coordinate system structure */
char *header;
double alfac, deltac;
float xdim, ydim;
int lhead,nbfits;
fitsfile *tf;


/* Variables para el tratamiento de errores */
int errflag;

char coorfile[51];


char pgdevice[100];

/* Variables para el fichero de salida */
char outfile[51];
FILE *of;


/* Variables para leer FITS */
char imagefile[51];
char errimagefile[51];


/* Variables para la funcion de dispersion */
char dispfile[51];
float *ldodisp;
float *position;
int ndisp;
float A,B,C;
float ldomin,ldomax,deltaldo;


/* Variables para la funcion respuesta */
float r_pix1,r_val1,r_delt1;
char respfile[51];
float *respuesta,*ldoresp;
int nresp;


/* Subrutinas */
float ldo2x(float ldo);
float x2ldo(float x);
float dxdl(float ldo);
float resp(float ldo);
void ReadDisp(char dispfile[]);
void ReadResp(char respfile[], int fluxcal);
void ReadOptions();
void LoadOptions(char filepar[]);
float intimapix_lin(float *ima, int nx, int ny, float xp, float yp, float nullval);


int main(int argc, char **argv)
{




  /* Variables para la extraccion de los espectros */
  float  xi,yi;

  fitsfile *spectrum;
  fitsfile *errspectrum;
  char specfile[51];
  char wavefile[51];
  char fluxfile[51];
  char errspecfile[51];
  char errwavefile[51];
  char errfluxfile[51];

  /* Variables para la extraccion de espectros pero en una sola imagen */

  char key[10];
  fitsfile *wholespec;
  fitsfile *wholeerrspec;
  fitsfile *wholewave;
  fitsfile *wholeerrwave;
  fitsfile *wholeflux;
  fitsfile *wholeerrflux;

  char     wholespecfile[51];
  char     wholewavefile[51];
  char     wholefluxfile[51];
  char     wholeerrspecfile[51];
  char     wholeerrwavefile[51];
  char     wholeerrfluxfile[51];
  long wholenaxes[2],firstpixel[2],lastpixel[2];

 


  float *spec,*x,*wavespec,*ldospec,*fluxspec=NULL,*dumspec;
  float *errspec=NULL,*errwavespec=NULL,*errfluxspec=NULL;
  float ldo;
  int ipix;
  int j,k;
  float crpix=1;
  float crval=0;
  float cdelt=1;
  int fluxcal=0;
  /* Hasta aqui */
  /* Variables para la extraccion de las cajitas */
/*   fitsfile *box; */
/*   long nbox1,nbox2; */
/*   char boxfile[51]; */

  /* Hasta aqui */
    


  /* Variables para el fichero BUHO */
  int  nobj;
  int i,itrueobj;
  float *xp,*yp;
/*   float *xep,*yep,*pa; */
  float *sky,*sig,*errsky;

  float *iobj;
  int *log;
/*   FILE *filebuho; */
  char snul[50],snul1[32],snul2[32];

  /* Variables para el catalogo */
  double *alfa=NULL,*delta=NULL;
/*   double mag[NMAX],magb[NMAX]; */
/*   double num[NMAX]; */
/*   int platenum[NMAX]; */
  float *imasky;
  float *imaerrsky;
  struct headfits hsky,herrsky;
  /* Varaibales para leer el FITS */
  int status=0;
  int nfound, anynull;
  fitsfile *image, *errimage;
  long naxes[2], fpixel,  npixels, ii;
  long errnaxes[2], errfpixel, errnpixels;
  float datamin, datamax, nullval;
  int errnfound;
  float *ima,*errima=NULL;

  /* Variables del dibujo */
/*   float factorarr=1; */
  float tr[6];
  float mean;
  float sigma;
  int nx,ny;
  int errnx,errny;
  int pg1;
  float fnul;
  char cnul;
  FILE *filepar;
/*   FILE *out; */


  char touchchar[500];

  long rand_i;
  int time_i;
  char tempfile[500];

  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  RNDInit();



  if(argc < 2) {
    ReadOptions();
    pg1=cpgopen("?");
  }
  else {
    if((filepar=fopen(argv[1],"r"))==NULL) {
      printf("ERROR: Can't open options file %s\n",argv[1]);
      exit(1);
    }
    LoadOptions( argv[1]);
    pg1=cpgopen(pgdevice);
      
 
  }




  /* Leo la imagen FITS de entrada */
/*   //  printf(" Input FITS file: "); */
/*   //  scanf("%s",imagefile); */
  if( ffopen(&image, imagefile, READONLY, &status)) fits_report_error(stderr,status);
  if(fits_read_keys_lng(image, "NAXIS", 1, 2, naxes, &nfound, &status)) fits_report_error(stderr,status);
  npixels=naxes[0]*naxes[1];
  if((ima=malloc(npixels*sizeof(float)))==NULL) printf("I cannot dimension ima of %ld elements \n",npixels);
  fpixel=1;
  nullval=0;
  datamin=1.0e30;
  datamax=-1.0e30;
  
  printf("...Reading image %s \n",imagefile);
  if(fits_read_img(image, TFLOAT, fpixel, npixels, &nullval, ima, &anynull, &status )) fits_report_error(stderr,status);
  printf("...Computing datamin and datamax \n");
  for (ii=0 ;ii<npixels;ii++) {
/*     //printf(" ima %d %e \n",ii,ima[ii]); */
    if( ima[ii]< datamin) datamin=ima[ii];
    if(ima[ii]> datamax) datamax=ima[ii];
  }
  printf(" Datamin %f Datamax %f \n",datamin,datamax);
  mean=StMedia(npixels,ima,&sigma);
  datamin=mean-sigma*2;
  datamax=mean+sigma*2;
  
  nx=naxes[0];
  ny=naxes[1];
  if(errflag) {
    if(ffopen(&errimage, errimagefile, READONLY, &status)) fits_report_error(stderr,status);
    if(fits_read_keys_lng(errimage, "NAXIS", 1, 2, errnaxes, &errnfound, &status)) fits_report_error(stderr,status);
    errnpixels=errnaxes[0]*errnaxes[1];
    if((errima=malloc(errnpixels*sizeof(float)))==NULL) printf("I cannot dimension errima of %ld elements \n",errnpixels);
    errfpixel=1;
    nullval=0;
/*     printf(" errnpixl %d\n",errnpixels); */
    printf("...Reading error image %s \n",errimagefile);
    if(fits_read_img(errimage, TFLOAT, errfpixel, errnpixels, &nullval, errima, &anynull, &status )) fits_report_error(stderr,status);
    errnx=errnaxes[0];
    errny=errnaxes[1];
    if(naxes[0] != errnaxes[0] || naxes[1] != errnaxes[1] ) {
      fprintf(stderr,"It seems that the input image and the error image are not comptible\n");
      fprintf(stderr,"Image dimensions %ld x %ld  Error image dimersions %ld x %ld\n",naxes[0],naxes[1],errnaxes[0],errnaxes[1]);
      exit(1);
    }
  }
/*   printf(" Salio de aqui\n"); */
  


  
  /* Dibujo la imagen  */
  
  fflush(NULL);
  setvbuf(stdin,"",_IOLBF,0);

  cpgask(0);
/*   //cpglab("pixel","pixel",""); */
  cpgwnad(0.,naxes[0],0.,naxes[1]);

/*   printf(" Es qui\n"); */
  
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  cpggray(ima,naxes[0],naxes[1],1,naxes[0],1,naxes[1],datamax,datamin,tr);
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
  
  
  /* Leo el fichero BUHO */
  
  
  
  nobj=FileNLin(catfile);
/*   printf(" YESSSS\n"); */
  if((xp=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension xp of %d elements \n",nobj);
  if((yp=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension xp of %d elements \n",nobj);
  if((sky=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension sky of %d elements \n",nobj);
  if((errsky=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension errsky of %d elements \n",nobj);
  if((sig=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension sig of %d elements \n",nobj);
  if((log=malloc(nobj*sizeof(int)))==NULL) printf("I cannot dimension xp of %d elements \n",nobj);
  if((iobj=malloc(nobj*sizeof(int)))==NULL) printf("I cannot dimension iobj of %d elements \n",nobj);
  printf("Number of objects %d from catalog %s\n",nobj,catfile);
  
  if(fitsrh(fitssky,&hsky)) {
    fprintf(stderr,"ExtracSpec: ERROR. Error reading image %s",fitssky);
    exit(1);
  }
  if(ffopen(&skyfits, fitssky, READONLY, &status)) fits_report_error(stderr,status);
  if((imasky=malloc(hsky.naxis1*hsky.naxis2*sizeof(float)))==NULL )  { printf("I cannot dimension imaksy of %d bytes",hsky.naxis1*hsky.naxis2);exit(1);} 
  if(fits_read_img(skyfits, TFLOAT, fpixel, hsky.naxis1*hsky.naxis2, &nullval, imasky, &anynull, &status )) fits_report_error(stderr,status);
  
  

  if(fitsrh(fitserrsky,&herrsky)) {
    fprintf(stderr,"ExtracSpec: ERROR. Error reading image %s",fitserrsky);
    exit(1);
  }
  if( ffopen(&errskyfits, fitserrsky, READONLY, &status)) fits_report_error(stderr,status);
  if((imaerrsky=malloc(hsky.naxis1*hsky.naxis2*sizeof(float)))==NULL )  { printf("I cannot dimension imaerrksy of %d bytes",hsky.naxis1*hsky.naxis2);exit(1);} 
  if(fits_read_img(errskyfits, TFLOAT, fpixel, herrsky.naxis1*herrsky.naxis2, &nullval, imaerrsky, &anynull, &status )) fits_report_error(stderr,status);


  ReadNumcol(catfile, 1,iobj,log,&nobj);
/*   printf(" Y deaqui\n"); */
  
  if(!isradec) {
/*     printf(" No tiene que estarqu \n"); */
    ReadNumcol(catfile, xcol,xp,log,&nobj);
    ReadNumcol(catfile, ycol,yp,log,&nobj);
    /*     ReadNumcol(catfile, skycol,sky,log,&nobj); */
    /*     ReadNumcol(catfile, errskycol,errsky,log,&nobj); */
    /* Termino de leerlo */
    for(i=0;i<nobj;i++ ) {
      fits_pv(imasky,&hsky,xp[nobj],yp[nobj],sky+i);
      fits_pv(imaerrsky,&hsky,xp[nobj],yp[nobj],errsky+i);
    }
  }
  else {
/*     printf(" pero si aui1 \n"); */
    if((alfa=malloc(nobj*sizeof(double)))==NULL) printf("I cannot dimension alfa of %d elements \n",nobj);
    if((delta=malloc(nobj*sizeof(double)))==NULL) printf("I cannot dimension delta of %d elements \n",nobj);
/*     printf(" Sabin\n"); */
    ReadWCScol(catfile, xcol,alfa,log,&nobj);
    ReadWCScol(catfile, ycol,delta,log,&nobj);
/*     printf(" Pasado \n"); */

    /* Esto crea una copia de la cabecera de la imagen para luego leerla con fitsrhead
       que no permite leer imagenes comprimidas
       Despues borro el fichero temporal temp_header */
    rand_i=random();
    time_i=(int)time(NULL);
    sprintf(tempfile,"temp%010d_%012ld_header.fits",time_i,rand_i);
    fprintf(stderr," Using %s as temp file \n",tempfile); 
    unlink(tempfile);
    ffinit(&tf,tempfile,&status);
    fits_copy_header(image, tf, &status);
    ffclos(tf,&status);
    if ((header = fitsrhead (tempfile, &lhead, &nbfits)) == NULL) {
      fprintf (stderr, "Cannot read FITS header of file %s\n", imagefile);
      exit(1);
    }
    unlink(tempfile);
    if((wcsim=wcsinit(header))==NULL) {
      printf(" No WCS information found in header\n");  
      exit(1);
    }
/*     printf(" Cafa\n"); */
    printf(" WCS information from header:\n");
    PrintWCS(header,1);
/*     printf(" SAlalad\n"); */
    pix2wcs(wcsim,naxes[0]/2., naxes[1]/2. ,&alfac, &deltac);
    xdim=naxes[0]*fabs((float)((*wcsim).xinc)*3600.);
    ydim=naxes[1]*fabs((float)((*wcsim).yinc)*3600.);
    /*     mag1=mag2=0; */
    /*     printf(" Searching objects round %f %f\n",alfac/15., deltac); */
    /*     printf(" Searching with box %f %f \n",xdim,ydim); */
    /*     nusno=uacread("UAC2",alfac,deltac,xdim/2./3600.,ydim/2./3600.,0.,1,2000.,wcsim->equinox,mag1,mag2,0,NMAX,num,alfa,delta,mag,magb,platenum,0); */
    
    /*     nobj=0; */
    /*     strcpy(usnofile,imagefile); */
    /*     strcat(usnofile,".extract.usno"); */
    /*     if((fusno=fopen(usnofile,"w"))==NULL) { */
    /*       printf("Cannot open file %s\n",usnofile); */
    /*       exit(1); */
    /*     } */
    /*     fprintf(fusno,"#USNO objects with extracted spectra\n"); */
    /*     fprintf(fusno,"#Num     RA(J2000)   DEC(J2000)      Mag      MagB    Plate N.   Xpos    Ypos\n"); */
    for(i=0;i<nobj;i++ ) {
      if(log[i]) {
	alfa[i]=alfa[i]*15.;
	wcs2pix(wcsim,alfa[i],delta[i],&xpos,&ypos,&off);
/* 	printf(" alfa %f delta %f\n",alfa[i],delta[i]); */
/* 	printf(" off %d xpos %f ypos %f\n",off,xpos,ypos); */
	if(!off && xpos>1+specpix && xpos<naxes[0]-specpix && ypos>1+spapix && ypos<naxes[1]-spapix) {
	  log[i]=1;
	  xp[i]         =(float)xpos;
	  yp[i]         =(float)ypos;
	  /* 	alfa[nobj]       =alfa[i]; */
	  /* 	delta[nobj]      =delta[i]; */
	  /* 	mag[nobj]        =mag[i]; */
	  /* 	num[nobj]        =num[i]; */
	  /* 	magb[nobj]       =magb[i]; */
	  /* 	platenum[nobj]   =platenum[i]; */
/* 	  printf(" xp %f yp %f\n",xp[i],yp[i]); */
	  fits_pv(imasky,&hsky,xp[i],yp[i],sky+i);
	  fits_pv(imaerrsky,&hsky,xp[i],yp[i],errsky+i);
	  
/* 	  printf(" Cielo %f \n",sky[i]); */
/* 	  exit(1); */
	}
	else log[i]=0;
	/* 	fprintf(fusno," %14.8f %s %s %6.2f %6.2f  %8d  %7.2f  %7.2f\n",num[i],snul1,snul2,mag[i],magb[i],platenum[i],xp[i],yp[i]); */
      }
    }
    /*     printf(" Number of objects to extract: %d\n",nobj); */
  }


  /* Creo el fichero de salida */

  if((of=fopen(outfile,"w"))==NULL) {
    printf("ERROR: Can't open out file %s\n",outfile);
    exit(1);
  }

  fprintf(of,"#Spectra file for image %s from catalogue file %s\n",imagefile,catfile);
  fprintf(of,"#                                raw_spectra_file                ");
  if(waveflag) fprintf(of,"                wave_calib_spectra_file                     ");
  if(fluxflag) fprintf(of,"              flux_calib_spectra_file            ");
  if(errflag) fprintf(of,"                        raw_err_spectra_file             ");
  if(errflag && waveflag) fprintf(of,"             wave_calib_err_spectra_file                  ");
  if(errflag && fluxflag) fprintf(of,"               flux_calib_err_spectra_file            ");
  fprintf(of,"      X_extraction     Y_extraction  ");
  if(isradec) fprintf(of,"        RA          DEC       ");
  fprintf(of,"     Sky      \n");
  

  /* Abro los ficheros FITS de espectros */

  /* Estas variables las necesito para crear los FITS de espectros */
  ldomin=x2ldo(0);
  ldomax=x2ldo(2*specpix-1);
  deltaldo=(ldomax-ldomin)/(2*specpix-1);
  ReadResp(respfile,   fluxcal);
  naxes[0]=specpix*2;
  naxes[1]=1;

  wholenaxes[0]=specpix*2;
  wholenaxes[1]=nobj;

  if(ismultspec) {
    
    /* Este es el de los espectros crudos */
    strcpy(wholespecfile,rawspec);
    strcat(wholespecfile,".fits");
    if(ffinit(&wholespec,wholespecfile,&status)) fits_report_error(stderr,status);
    if(status) fits_report_error(stderr,status);
    fits_create_img(wholespec, -32,2,wholenaxes,&status);
    if(status) fits_report_error(stderr,status);
    fits_write_key(wholespec,TFLOAT,"CRPIX1",&crpix,"Reference pixel",&status);
    fits_write_key(wholespec,TFLOAT,"CRVAL1",&crval,"Coordinate at reference pixel",&status);
    fits_write_key(wholespec,TFLOAT,"CDELT1",&cdelt,"Coordinate increment per pixel",&status);
    fits_write_key(wholespec,TINT,"NSCANS",&spapix,"Number of scans added",&status);
    fits_write_key(wholespec,TSTRING,"IMAGE",imagefile,"Image where the spectra were extracted",&status);
    fits_write_key(wholespec,TSTRING,"ERRIMAGE",errimagefile,"Error image for IMAGE",&status);
    fits_write_key(wholespec,TSTRING,"CATALOG",catfile,"Catalogue with objects",&status);
    fits_write_key(wholespec,TINT,"NSPEC",&nobj,"Number of objects from catalogue",&status);
    if(status) fits_report_error(stderr,status);
    fits_set_hdrsize(wholespec,nobj*5,&status);
    if(status) fits_report_error(stderr,status);

    if(errflag) {
      
      /* Este es el de los errores de los espectros crudos */
      strcpy(wholeerrspecfile,errrawspec);
      strcat(wholeerrspecfile,".fits");
      if(ffinit(&wholeerrspec,wholeerrspecfile,&status)) fits_report_error(stderr,status);
      if(status) fits_report_error(stderr,status);
      fits_create_img(wholeerrspec, -32,2,wholenaxes,&status);
      if(status) fits_report_error(stderr,status);
      fits_write_key(wholeerrspec,TFLOAT,"CRPIX1",&crpix,"Reference pixel",&status);
      fits_write_key(wholeerrspec,TFLOAT,"CRVAL1",&crval,"Coordinate at reference pixel",&status);
      fits_write_key(wholeerrspec,TFLOAT,"CDELT1",&cdelt,"Coordinate increment per pixel",&status);
      fits_write_key(wholeerrspec,TINT,"NSCANS",&spapix,"Number of scans added",&status);
      fits_write_key(wholeerrspec,TSTRING,"IMAGE",imagefile,"Image where the spectra were extracted",&status);
      fits_write_key(wholeerrspec,TSTRING,"ERRIMAGE",errimagefile,"Error image for IMAGE",&status);
      fits_write_key(wholeerrspec,TSTRING,"CATALOG",catfile,"Catalogue with objects",&status);
      fits_write_key(wholeerrspec,TINT,"NSPEC",&nobj,"Number of objects from catalogue",&status);
      fits_set_hdrsize(wholeerrspec,nobj*5,&status);
      if(status) fits_report_error(stderr,status);
      
    }

    if(waveflag) {
      
      /* Este es el de los espectros calibrados en longitud de onda */
      strcpy(wholewavefile,wcspec);
      strcat(wholewavefile,".fits");
      if(ffinit(&wholewave,wholewavefile,&status)) fits_report_error(stderr,status);
      if(status) fits_report_error(stderr,status);
      fits_create_img(wholewave, -32,2,wholenaxes,&status);
      if(status) fits_report_error(stderr,status);
      fits_write_key(wholewave,TFLOAT,"CRPIX1",&crpix,"Reference pixel",&status);
      fits_write_key(wholewave,TFLOAT,"CRVAL1",&ldomin,"Coordinate at reference pixel",&status);
      fits_write_key(wholewave,TFLOAT,"CDELT1",&deltaldo,"Coordinate increment per pixel",&status);
      fits_write_key(wholewave,TINT,"NSCANS",&spapix,"Number of scans added",&status);
      fits_write_key(wholewave,TSTRING,"IMAGE",imagefile,"Image where the spectra were extracted",&status);
      fits_write_key(wholewave,TSTRING,"ERRIMAGE",errimagefile,"Error image for IMAGE",&status);
      fits_write_key(wholewave,TSTRING,"RESPONSE",respfile,"Fits file used as response function",&status);
      fits_write_key(wholewave,TSTRING,"CATALOG",catfile,"Catalogue with objects",&status);
      fits_write_key(wholewave,TINT,"NSPEC",&nobj,"Number of objects from catalogue",&status);
      fits_set_hdrsize(wholewave,nobj*5,&status);
      if(status) fits_report_error(stderr,status);
      
      if(errflag) {
	/* Este es el de los errores de los espectros calibrados en ldo */
	strcpy(wholeerrwavefile,errwcspec);
	strcat(wholeerrwavefile,".fits");
	if(ffinit(&wholeerrwave,wholeerrwavefile,&status)) fits_report_error(stderr,status);
	if(status) fits_report_error(stderr,status);
	fits_create_img(wholeerrwave, -32,2,wholenaxes,&status);
	if(status) fits_report_error(stderr,status);
	fits_write_key(wholeerrwave,TFLOAT,"CRPIX1",&crpix,"Reference pixel",&status);
	fits_write_key(wholeerrwave,TFLOAT,"CRVAL1",&ldomin,"Coordinate at reference pixel",&status);
	fits_write_key(wholeerrwave,TFLOAT,"CDELT1",&deltaldo,"Coordinate increment per pixel",&status);
	fits_write_key(wholeerrwave,TINT,"NSCANS",&spapix,"Number of scans added",&status);
	fits_write_key(wholeerrwave,TSTRING,"IMAGE",imagefile,"Image where the spectra were extracted",&status);
	fits_write_key(wholeerrwave,TSTRING,"ERRIMAGE",errimagefile,"Error image for IMAGE",&status);
	fits_write_key(wholeerrwave,TSTRING,"RESPONSE",respfile,"Fits file used as response function",&status);
	fits_write_key(wholeerrwave,TSTRING,"CATALOG",catfile,"Catalogue with objects",&status);
	fits_write_key(wholeerrwave,TINT,"NSPEC",&nobj,"Number of objects from catalogue",&status);
	fits_set_hdrsize(wholeerrwave,nobj*5,&status);
	if(status) fits_report_error(stderr,status);
      }
    }
    
    if(fluxflag) {
      /* Este es el de los espectros calibrados en longitud de onda y en flujo*/
      strcpy(wholefluxfile,fcspec);
      strcat(wholefluxfile,".fits");
      if(ffinit(&wholeflux,wholefluxfile,&status)) fits_report_error(stderr,status);
      if(status) fits_report_error(stderr,status);
      fits_create_img(wholeflux, -32,2,wholenaxes,&status);
      if(status) fits_report_error(stderr,status);
      fits_write_key(wholeflux,TFLOAT,"CRPIX1",&crpix,"Reference pixel",&status);
      fits_write_key(wholeflux,TFLOAT,"CRVAL1",&ldomin,"Coordinate at reference pixel",&status);
      fits_write_key(wholeflux,TFLOAT,"CDELT1",&deltaldo,"Coordinate increment per pixel",&status);
      fits_write_key(wholeflux,TINT,"NSCANS",&spapix,"Number of scans added",&status);
      fits_write_key(wholeflux,TSTRING,"IMAGE",imagefile,"Image where the spectra were extracted",&status);
      fits_write_key(wholeflux,TSTRING,"ERRIMAGE",errimagefile,"Error image for IMAGE",&status);
      fits_write_key(wholeflux,TSTRING,"RESPONSE",respfile,"Fits file used as response function",&status);
      fits_write_key(wholeflux,TSTRING,"CATALOG",catfile,"Catalogue with objects",&status);
      fits_write_key(wholeflux,TINT,"NSPEC",&nobj,"Number of objects from catalogue",&status);
      fits_set_hdrsize(wholeflux,nobj*5,&status);
      if(status) fits_report_error(stderr,status);
      
      if(errflag) {
	/* Este es el de los errores de los espectros calibrados en ldo */
	strcpy(wholeerrfluxfile,errfcspec);
	strcat(wholeerrfluxfile,".fits");
	if(ffinit(&wholeerrflux,wholeerrfluxfile,&status)) fits_report_error(stderr,status);
	if(status) fits_report_error(stderr,status);
	fits_create_img(wholeerrflux, -32,2,wholenaxes,&status);
	if(status) fits_report_error(stderr,status);
	fits_write_key(wholeerrflux,TFLOAT,"CRPIX1",&crpix,"Reference pixel",&status);
	fits_write_key(wholeerrflux,TFLOAT,"CRVAL1",&ldomin,"Coordinate at reference pixel",&status);
	fits_write_key(wholeerrflux,TFLOAT,"CDELT1",&deltaldo,"Coordinate increment per pixel",&status);
	fits_write_key(wholeerrflux,TINT,"NSCANS",&spapix,"Number of scans added",&status);
	fits_write_key(wholeerrflux,TSTRING,"IMAGE",imagefile,"Image where the spectra were extracted",&status);
	fits_write_key(wholeerrflux,TSTRING,"ERRIMAGE",errimagefile,"Error image for IMAGE",&status);
	fits_write_key(wholeerrflux,TSTRING,"RESPONSE",respfile,"Fits file used as response function",&status);
	fits_write_key(wholeerrflux,TSTRING,"CATALOG",catfile,"Catalogue with objects",&status);
	fits_write_key(wholeerrflux,TINT,"NSPEC",&nobj,"Number of objects from catalogue",&status);
	fits_set_hdrsize(wholeerrflux,nobj*5,&status);
	if(status) fits_report_error(stderr,status);
      }
    }
    
    
    firstpixel[0]=1;lastpixel[0]=specpix*2;
    
  }


  /* EXTRACCION ESPECTROS 
     --------------------
     --------------------*/
  
  if((spec=malloc(specpix*2*sizeof(float)))==NULL) printf("I cannot dimension spec of %d elements \n",specpix*2);
  if(errflag) if((errspec=malloc(specpix*2*sizeof(float)))==NULL) printf("I cannot dimension errspec of %d elements \n",specpix*2);
  if((x=malloc(specpix*2*sizeof(float)))==NULL) printf("I cannot dimension x of %d elements \n",specpix*2);
  setvbuf(stdin,"",_IOLBF,0);
  /*   pg2=cpgopen(pgdevice); */
  cpgask(0);
  /*   printf("dimensio x %d specpix%d \n",nx,specpix); */
  for(i=0;i<specpix*2;i++) x[i]=(float)i;

  printf(" Number of spectra to extract: %d\n",nobj);

  printf("\n");
  itrueobj=0;
  for(i=0;i<nobj;i++) {
    for(j=0;j<specpix*2;j++) spec[j]=0;
    if(errflag) for(j=0;j<specpix*2;j++) errspec[j]=0;
    if(log[i]) {
      itrueobj++;
      strcpy(specfile,rawspec);
      strcpy(errspecfile,errrawspec);
      /*       if(fromcat) sprintf(snul,"%12.8f",num[i]); */
      sprintf(snul,"%07d",(int)iobj[i]);
      strcat(specfile,snul);
      strcat(specfile,".fits");
      strcat(errspecfile,snul);
      strcat(errspecfile,".fits");
      /*       if(!fromcat) printf(" Creating spec %f with name %s at position [%f,%f] \n",iobj[i],specfile,xp[i],yp[i]); */
      /*       else printf(" Creating spec of object %g with name %s at position [%f,%f] \n",num[i],specfile,xp[i],yp[i]); */
      printf("%07d/%07d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",i+1,nobj); 
/*       printf("%07d/%07d",i+1,nobj); */
/*       printf(" %07d/%07d Image %s Object %d .. ",i+1,nobj,imagefile,iobj[i]); */

      /*       else printf(" %6d/%6d Image %s Object %12.8f .. ",i+1,nobj,imagefile,num[i]); */
/*       printf(" Extracting at %f %f \n",xp[i],yp[i]); */
      for(j=0;j<specpix*2;j++) {
	for(k=0;k<spapix*2;k++) {
	  /* Cuidado!! xp esta en pixels y hay que pasarlo a 
	     dimensiones X Y en X!!. De ahi el -1 */

	  if(dispaxis==1){ 
	    xi=xp[i]-specpix+j+0.5;
	    yi=yp[i]-spapix+k+0.5;
	  }
	  else {
	    xi=(xp[i]-spapix+k+0.5);
	    yi=(yp[i]-specpix+j+0.5);
	  }
/* 	  printf(" Adding at %f %f \n",xi,yi); */
	  spec[j]+=intimapix_lin(ima,nx,ny,xi,yi,sky[i]);
/* 	  printf(" %d spec[j] %f\n",j,spec[j]); */
	  if(errflag) {
	    errspec[j]+=intimapix_lin(errima,nx,ny,xi,yi,errsky[i])*intimapix_lin(errima,nx,ny,xi,yi,errsky[i]);
	  }
	}
	spec[j]=(spec[j]-spapix*2*sky[i]);
	if(errflag) errspec[j]=sqrt(errspec[j]+spapix*2*errsky[i]*spapix*2*errsky[i]);
      }
      datamin=1.0e30;
      datamax=-1.0e30;
      for (ii=0 ;ii<specpix*2;ii++) {
	if(spec[ii]< datamin) datamin=spec[ii];
	if(spec[ii]> datamax) datamax=spec[ii];
      }
      /* Dibujo el espectro */
      
      
      cpgswin(0.,specpix*2,datamin-(datamax-datamin)*.2,datamax+(datamax-datamin)*.2);
      
      cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
      cpgbin(specpix*2,x,spec,1);
      
      /* Salvando el espectro RAW */
      /*       printf(" Saving spec %s\n",specfile); */
/*       printf("R"); */
      status=0;
      
      if(ismultspec) {
	/* Salvando de la nueva manera */
	firstpixel[1]=itrueobj;lastpixel[1]=itrueobj;
	fits_write_subset_flt(wholespec,1,2,wholenaxes,firstpixel,lastpixel,spec,&status);
	if(status) fits_report_error(stderr,status);
	sprintf(key,"N%07d",itrueobj);
	fits_write_key(wholespec,TFLOAT,key,iobj+i,"Number of object in catalogue",&status);
	sprintf(key,"X%07d",itrueobj);
	fits_write_key(wholespec,TFLOAT,key,xp+i,"X position of the object in the image",&status);
	sprintf(key,"Y%07d",itrueobj);
	fits_write_key(wholespec,TFLOAT,key,yp+i,"Y position of the object in the image",&status);
	if(isradec) {
	  ra2str(snul1,32,alfa[i],3);
	  dec2str(snul2,32,delta[i],2);
	  sprintf(key,"R%07d",itrueobj);
	  fits_write_key(wholespec,TSTRING,key,snul1,"R.A. of the extracted object",&status);
	  sprintf(key,"D%07d",itrueobj);
	  fits_write_key(wholespec,TSTRING,key,snul2,"Declination of the extracted object",&status);
	}
	
      }
      else {
	
	/* Salvando de la forma antigua */
	if(ffinit(&spectrum,specfile,&status)) fits_report_error(stderr,status);
	if(status) fits_report_error(stderr,status);
	fits_create_img(spectrum, -32,1,naxes,&status);
	if(status) fits_report_error(stderr,status);
	/*       printf(" ANTESD EET\n"); */
	fits_write_img(spectrum,TFLOAT,1,specpix*2,spec,&status);
	if(status) fits_report_error(stderr,status);
	
	fits_write_key(spectrum,TFLOAT,"XCPIX",xp+i,"X position of the object in the image",&status);
	fits_write_key(spectrum,TFLOAT,"YCPIX",yp+i,"Y position of the object in the image",&status);
	fits_write_key(spectrum,TFLOAT,"SKY",sky+i,"Substracted sky near the object",&status);
	/*       fits_write_key(spectrum,TFLOAT,"SKY_SIG",errsig+i,"Sigma of sky near the object",&status); */
	if(isradec) {
	  ra2str(snul1,32,alfa[i],3);
	  dec2str(snul2,32,delta[i],2);
	  fits_write_key(spectrum,TSTRING,"RA" ,snul1,"R.A. of the extracted object",&status);
	  fits_write_key(spectrum,TSTRING,"DEC",snul2,"Declination of the extracted object",&status);
	}
	fits_write_key(spectrum,TFLOAT,"CRPIX1",&crpix,"Reference pixel",&status);
	fits_write_key(spectrum,TFLOAT,"CRVAL1",&crval,"Coordinate at reference pixel",&status);
	fits_write_key(spectrum,TFLOAT,"CDELT1",&cdelt,"Coordinate increment per pixel",&status);
	fits_write_key(spectrum,TINT,"NSCANS",&spapix,"Number of scans added",&status);
	fits_write_key(spectrum,TSTRING,"IMAGE",imagefile,"Image where the spec was extracted from",&status);
	/*       fits_write_key(spec,TSTRING,"RESPONSE",respfile,"Fits file used as response function",&status); */
	if(status) fits_report_error(stderr,status);
	fits_write_key(spectrum,TFLOAT,"NOBJ",iobj+i,"Number of object in catalogue file",&status);
	if(status) fits_report_error(stderr,status);
	fits_close_file(spectrum,&status);
	if(status) fits_report_error(stderr,status);
	/*            printf(" Making wavelength calibration...\n"); */
	if(itrueobj==-11) {
	  fits_close_file(wholespec,&status);
	  exit(1);
	}
      }
      if(errflag) {
	if(ismultspec) {
	  /* Salvando de la nueva manera */
	  firstpixel[1]=itrueobj;lastpixel[1]=itrueobj;
	  fits_write_subset_flt(wholeerrspec,1,2,wholenaxes,firstpixel,lastpixel,errspec,&status);
	  sprintf(key,"N%07d",itrueobj);
	  fits_write_key(wholeerrspec,TFLOAT,key,iobj+i,"Number of object in catalogue",&status);
	  sprintf(key,"X%07d",itrueobj);
	  fits_write_key(wholeerrspec,TFLOAT,key,xp+i,"X position of the object in the image",&status);
	  sprintf(key,"Y%07d",itrueobj);
	  fits_write_key(wholeerrspec,TFLOAT,key,yp+i,"Y position of the object in the image",&status);
	  if(isradec) {
	    ra2str(snul1,32,alfa[i],3);
	    dec2str(snul2,32,delta[i],2);
	    sprintf(key,"R%07d",itrueobj);
	    fits_write_key(wholeerrspec,TSTRING,key,snul1,"R.A. of the extracted object",&status);
	    sprintf(key,"D%07d",itrueobj);
	    fits_write_key(wholeerrspec,TSTRING,key,snul2,"Declination of the extracted object",&status);
	  }
	}
	else {
	  if(ffinit(&errspectrum,errspecfile,&status)) fits_report_error(stderr,status);
	  if(status) fits_report_error(stderr,status);
	  fits_create_img(errspectrum, -32,1,naxes,&status);
	  if(status) fits_report_error(stderr,status);
	  fits_write_img(errspectrum,TFLOAT,1,specpix*2,errspec,&status);
	  if(status) fits_report_error(stderr,status);
	  fits_write_key(errspectrum,TFLOAT,"XCPIX",xp+i,"X position of the object in the image",&status);
	  fits_write_key(errspectrum,TFLOAT,"YCPIX",yp+i,"Y position of the object in the image",&status);
	  fits_write_key(errspectrum,TFLOAT,"SKY",sky+i,"Substracted sky near the object",&status);
	  /* 	  fits_write_key(errspectrum,TFLOAT,"SKY_SIG",sig+i,"Sigma of sky near the object",&status); */
	  if(isradec) {
	    ra2str(snul1,32,alfa[i],3);
	    dec2str(snul2,32,delta[i],2);
	    fits_write_key(errspectrum,TSTRING,"RA" ,snul1,"R.A. of the extracted object",&status);
	    fits_write_key(errspectrum,TSTRING,"DEC",snul2,"Declination of the extracted object",&status);
	  }
	  fits_write_key(errspectrum,TFLOAT,"CRPIX1",&crpix,"Reference pixel",&status);
	  fits_write_key(errspectrum,TFLOAT,"CRVAL1",&crval,"Coordinate at reference pixel",&status);
	  fits_write_key(errspectrum,TFLOAT,"CDELT1",&cdelt,"Coordinate increment per pixel",&status);
	  fits_write_key(errspectrum,TINT,"NSCANS",&spapix,"Number of scans added",&status);
	  fits_write_key(errspectrum,TSTRING,"IMAGE",imagefile,"Image where the spec was extracted from",&status);
	  if(status) fits_report_error(stderr,status);
	  fits_write_key(errspectrum,TFLOAT,"NOBJ",iobj+i,"Number of object in buho file",&status);
	  /* 	else fits_write_key(errspectrum,TDOUBLE,"USNOID",num+i,"USNO identification name",&status); */
	  if(status) fits_report_error(stderr,status);
	  fits_close_file(errspectrum,&status);
	  if(status) fits_report_error(stderr,status);
	}
      }
      if(waveflag) {
/* 	printf("W"); */
      WAVCAL:
	ldomin=x2ldo(0.);
	ldomax=x2ldo(2*specpix-1.);
	deltaldo=(ldomax-ldomin)/(2*specpix-1.);
      
	/*        Now we create the wavelength calibrated spectra */
	
	if((wavespec=malloc(specpix*2*sizeof(float)))==NULL) printf("I cannot dimension wavespec of %d elements \n",specpix*2);
	if((dumspec=malloc(specpix*2*sizeof(float)))==NULL) printf("I cannot dimension wavespec of %d elements \n",specpix*2);
	if((fluxspec=malloc(specpix*2*sizeof(float)))==NULL) printf("I cannot dimension fluxspec of %d elements \n",specpix*2);
	if((ldospec =malloc(specpix*2*sizeof(float)))==NULL) printf("I cannot dimension ldospec  of %d elements \n",specpix*2);
	if(errflag) {
	  if((errwavespec=malloc(specpix*2*sizeof(float)))==NULL) printf("I cannot dimension errwavespec of %d elements \n",specpix*2);
	  if((errfluxspec=malloc(specpix*2*sizeof(float)))==NULL) printf("I cannot dimension errfluxspec of %d elements \n",specpix*2);
	}
	for(j=0;j<specpix*2;j++) {
	  ldo=ldomin+j*deltaldo;
	  ipix=(int)(ldo2x(ldo)+0.5);
	  /*El 0.5 es para que el redondeo sea mejor (porque (int)0.9 =0 !! */
	  /* 	ipix2=(int)(ldo2x(ldo+4*deltaldo)); */
	  /* 	dxdl: derivada de x con respecto a lambda para transformar. */
	  if(ipix<0) ipix=0;
	  if(ipix>(specpix*2-1)) ipix=specpix*2-1;
	  /* 		printf(" ldo %f pix %d \n",ldo,ipix,dxdl); */
	  /* 		printf(" HA ver pix %d ldo %f pix %f \n",-specpix,ldomin,ldo2x(ldomin)); */
	  wavespec[j]=spec[ipix]*dxdl(ldo);
	  if(errflag) errwavespec[j]=errspec[ipix]*dxdl(ldo);
	  /* 	wavespec[j]=spec[ipix]*1.; */
	  /* 	printf(" Ldo %f wces %f ",ldo,wavespec[j]); */
	  if(resp(ldo)<=0.1) fluxspec[j]=0;
	  else fluxspec[j]=wavespec[j]/resp(ldo);
	  if(errflag) {
	    if(resp(ldo)<=0.1) errfluxspec[j]=0;
	    else errfluxspec[j]=errwavespec[j]/resp(ldo);
	  }
	  /* 	printf(" Raw %f Wave %f Flues %f Resp %f\n",spec[ipix],wsavespec[j],fluxspec[j],resp(ldo)); */
	}
	
	
	if(i==-1) {
	  
	  datamin=1.0e30;
	  datamax=-1.0e30;
	  for (ii=0 ;ii<specpix*2;ii++) {
	    if(wavespec[ii]< datamin) datamin=wavespec[ii];
	    if(wavespec[ii]> datamax) datamax=wavespec[ii];
	    ldo=ldomin+ii*deltaldo;
	    dumspec[ii]=resp(ldo);
	  }
	  
	  cpgpage();
	  printf(" Con A: %f\n",A);
	  printf(" xldo(0) %f xldo(80) %f ldox(6500) %f ldox(7000) %f\n",x2ldo(0.),x2ldo(80.),ldo2x(6500.),ldo2x(7000.));
	  cpgswin(ldomin,ldomax,datamin,datamax);
	  
	  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
	  cpgswin(0.,specpix*2,datamin,datamax);
	  cpgbin(specpix*2,x,wavespec,1);
	  
	  cpgswin(0.,specpix*2,0.,1.);
	  cpgbin(specpix*2,x,dumspec,1);
	  cpgcurs(&fnul,&fnul,&cnul);
	  printf(" cnul %c\n",cnul);
	  if(cnul=='A') A+=10;
	  if(cnul=='X') A-=10;
	  
	  goto WAVCAL;
	}
	/*             printf(" Creating calibrated spec %f with name %s at position [%f,%f] \n",iobj[i],wavefile,xp[i],yp[i]); */
	
	
	/*       Now saving wavelength calibrated spec */
	strcpy(wavefile,wcspec);
	strcpy(errwavefile,errwcspec);
	/*       if(fromcat) sprintf(snul,"%12.8f",num[i]); */
	sprintf(snul,"%07d",(int)iobj[i]);
	strcat(wavefile,snul);
	strcat(wavefile,".fits");
	strcat(errwavefile,snul);
	strcat(errwavefile,".fits");
	crpix=1;
	/*       printf(" saving calibrated spec %s\n",wavefile); */
	printf("F");
	
	if(ismultspec) {
	  
	  firstpixel[1]=itrueobj;lastpixel[1]=itrueobj;
	  fits_write_subset_flt(wholewave,1,2,wholenaxes,firstpixel,lastpixel,wavespec,&status);
	  sprintf(key,"N%07d",itrueobj);
	  fits_write_key(wholewave,TFLOAT,key,iobj+i,"Number of object in catalogue",&status);
	  sprintf(key,"X%07d",itrueobj);
	  fits_write_key(wholewave,TFLOAT,key,xp+i,"X position of the object in the image",&status);
	  sprintf(key,"Y%07d",itrueobj);
	  fits_write_key(wholewave,TFLOAT,key,yp+i,"Y position of the object in the image",&status);
	  if(isradec) {
	    ra2str(snul1,32,alfa[i],3);
	    dec2str(snul2,32,delta[i],2);
	    sprintf(key,"R%07d",itrueobj);
	    fits_write_key(wholewave,TSTRING,key,snul1,"R.A. of the extracted object",&status);
	    sprintf(key,"D%07d",itrueobj);
	    fits_write_key(wholewave,TSTRING,key,snul2,"Declination of the extracted object",&status);
	  }
	  
	}
	else {
	  
	  ffinit(&spectrum,wavefile,&status);
	  /*       printf(" Y AQUI!!\n"); */
	  if(status) fits_report_error(stderr,status);
	  fits_create_img(spectrum, -32,1,naxes,&status);
	  if(status) fits_report_error(stderr,status);
	  fits_write_img(spectrum,TFLOAT,1,specpix*2,wavespec,&status);
	  if(status) fits_report_error(stderr,status);
	  fits_write_key(spectrum,TFLOAT,"XCPIX",xp+i,"X position of the object in the image",&status);
	  fits_write_key(spectrum,TFLOAT,"YCPIX",yp+i,"X position of the object in the image",&status);
	  fits_write_key(spectrum,TFLOAT,"SKY",sky+i,"Substracted sky near the object",&status);
	  /* 	fits_write_key(spectrum,TFLOAT,"SKY_SIG",sig+i,"Sigma of sky near the object",&status); */
	  
	  if(isradec) {
	    ra2str(snul1,32,alfa[i],3);
	    dec2str(snul2,32,delta[i],2);
	    fits_write_key(spectrum,TSTRING,"RA" ,snul1,"R.A. of the extracted object",&status);
	    fits_write_key(spectrum,TSTRING,"DEC",snul2,"Declination of the extracted object",&status);
	  }
	  
	  fits_write_key(spectrum,TFLOAT,"CRPIX1",&crpix,"Reference pixel",&status);
	  fits_write_key(spectrum,TFLOAT,"CRVAL1",&ldomin,"Coordinate at reference pixel",&status);
	  fits_write_key(spectrum,TFLOAT,"CDELT1",&deltaldo,"Coordinate increment per pixel",&status);
	  fits_write_key(spectrum,TINT,"NSCANS",&spapix,"Number of scans added",&status);
	  fits_write_key(spectrum,TSTRING,"IMAGE",imagefile,"Image where the spectrum was extracted from",&status);
	  fits_write_key(spectrum,TSTRING,"RESPONSE",respfile,"Fits file used as response function",&status);
	  fits_write_key(spectrum,TFLOAT,"NOBJ",iobj+i,"Number of object in buho file",&status);
	  /*       else fits_write_key(spectrum,TDOUBLE,"USNOID",num+i,"USNO identification name",&status); */
	  /*       fits_write_key(spec,TSTRING,"DISPFILE",dispfile,"Dispersion law used to wavelength calibration",&status); */
	  fits_close_file(spectrum,&status);
	  if(status) fits_report_error(stderr,status);
	}
	
	if(errflag) {
	  if(ismultspec) {	  
	    fits_write_subset_flt(wholeerrwave,1,2,wholenaxes,firstpixel,lastpixel,errwavespec,&status);
	    sprintf(key,"N%07d",itrueobj);
	    fits_write_key(wholeerrwave,TFLOAT,key,iobj+i,"Number of object in catalogue",&status);
	    sprintf(key,"X%07d",itrueobj);
	    fits_write_key(wholeerrwave,TFLOAT,key,xp+i,"X position of the object in the image",&status);
	    sprintf(key,"Y%07d",itrueobj);
	    fits_write_key(wholeerrwave,TFLOAT,key,yp+i,"Y position of the object in the image",&status);
	    if(isradec) {
	      ra2str(snul1,32,alfa[i],3);
	      dec2str(snul2,32,delta[i],2);
	      sprintf(key,"R%07d",itrueobj);
	      fits_write_key(wholeerrwave,TSTRING,key,snul1,"R.A. of the extracted object",&status);
	      sprintf(key,"D%07d",itrueobj);
	      fits_write_key(wholeerrwave,TSTRING,key,snul2,"Declination of the extracted object",&status);
	    }
	  }
	
	  else {
	    
	    ffinit(&errspectrum,errwavefile,&status);
	    if(status) fits_report_error(stderr,status);
	    fits_create_img(errspectrum, -32,1,naxes,&status);
	    if(status) fits_report_error(stderr,status);
	    fits_write_img(errspectrum,TFLOAT,1,specpix*2,errwavespec,&status);
	    if(status) fits_report_error(stderr,status);
	    fits_write_key(errspectrum,TFLOAT,"XCPIX",xp+i,"X position of the object in the image",&status);
	    fits_write_key(errspectrum,TFLOAT,"YCPIX",yp+i,"X position of the object in the image",&status);
	    fits_write_key(errspectrum,TFLOAT,"SKY",sky+i,"Substracted sky near the object",&status);
	    /* 	  fits_write_key(errspectrum,TFLOAT,"SKY_SIG",sig+i,"Sigma of sky near the object",&status); */
	    
	    if(isradec) {
	      ra2str(snul1,32,alfa[i],3);
	      dec2str(snul2,32,delta[i],2);
	      fits_write_key(errspectrum,TSTRING,"RA" ,snul1,"R.A. of the extracted object",&status);
	      fits_write_key(errspectrum,TSTRING,"DEC",snul2,"Declination of the extracted object",&status);
	    }
	    
	    fits_write_key(errspectrum,TFLOAT,"CRPIX1",&crpix,"Reference pixel",&status);
	    fits_write_key(errspectrum,TFLOAT,"CRVAL1",&ldomin,"Coordinate at reference pixel",&status);
	    fits_write_key(errspectrum,TFLOAT,"CDELT1",&deltaldo,"Coordinate increment per pixel",&status);
	    fits_write_key(errspectrum,TINT,"NSCANS",&spapix,"Number of scans added",&status);
	    fits_write_key(errspectrum,TSTRING,"IMAGE",imagefile,"Image where the spectrum was extracted from",&status);
	    fits_write_key(errspectrum,TSTRING,"RESPONSE",respfile,"Fits file used as response function",&status);
	    fits_write_key(errspectrum,TFLOAT,"NOBJ",iobj+i,"Number of object in buho file",&status);
	    /* 	else fits_write_key(errspectrum,TDOUBLE,"USNOID",num+i,"USNO identification name",&status); */
	    fits_close_file(errspectrum,&status);
	    if(status) fits_report_error(stderr,status);
	    
	  }
	}
      }
      
      if(fluxflag) {
	
	if(ismultspec) {
	  firstpixel[1]=itrueobj;lastpixel[1]=itrueobj;
	  fits_write_subset_flt(wholeflux,1,2,wholenaxes,firstpixel,lastpixel,fluxspec,&status);
	  sprintf(key,"N%07d",itrueobj);
	  fits_write_key(wholeflux,TFLOAT,key,iobj+i,"Number of object in catalogue",&status);
	  sprintf(key,"X%07d",itrueobj);
	  fits_write_key(wholeflux,TFLOAT,key,xp+i,"X position of the object in the image",&status);
	  sprintf(key,"Y%07d",itrueobj);
	  fits_write_key(wholeflux,TFLOAT,key,yp+i,"Y position of the object in the image",&status);
	  if(isradec) {
	    ra2str(snul1,32,alfa[i],3);
	    dec2str(snul2,32,delta[i],2);
	    sprintf(key,"R%07d",itrueobj);
	    fits_write_key(wholeflux,TSTRING,key,snul1,"R.A. of the extracted object",&status);
	    sprintf(key,"D%07d",itrueobj);
	    fits_write_key(wholeflux,TSTRING,key,snul2,"Declination of the extracted object",&status);
	  }
	  
	}
	else {
	  strcpy(fluxfile,fcspec);
	  /* 	if(fromcat) sprintf(snul,"%12.8f",num[i]); */
	  sprintf(snul,"%08d",(int)iobj[i]);
	  strcat(fluxfile,snul);
	  strcat(fluxfile,".fits");
	  strcat(errfluxfile,snul);
	  strcat(errfluxfile,".fits");
	  crpix=1;
	  ffinit(&spectrum,fluxfile,&status);
	  if(status) fits_report_error(stderr,status);
	  fits_create_img(spectrum, -32,1,naxes,&status);
	  if(status) fits_report_error(stderr,status);
	  fits_write_img(spectrum,TFLOAT,1,specpix*2,fluxspec,&status);
	  if(status) fits_report_error(stderr,status);
	  fits_write_key(spectrum,TFLOAT,"XCPIX",xp+i,"X position of the object in the image",&status);
	  fits_write_key(spectrum,TFLOAT,"YCPIX",yp+i,"X position of the object in the image",&status);
	  fits_write_key(spectrum,TFLOAT,"SKY",sky+i,"Substracted sky near the object",&status);
	  /* 	  fits_write_key(spectrum,TFLOAT,"SKY_SIG",sig+i,"Sigma of sky near the object",&status); */
	  if(isradec) {
	    ra2str(snul1,32,alfa[i],3);
	    dec2str(snul2,32,delta[i],2);
	    fits_write_key(spectrum,TSTRING,"RA" ,snul1,"R.A. of the extracted object",&status);
	    fits_write_key(spectrum,TSTRING,"DEC",snul2,"Declination of the extracted object",&status);
	  }
	  fits_write_key(spectrum,TFLOAT,"CRPIX1",&crpix,"Reference pixel",&status);
	  fits_write_key(spectrum,TFLOAT,"CRVAL1",&ldomin,"Coordinate at reference pixel",&status);
	  fits_write_key(spectrum,TFLOAT,"CDELT1",&deltaldo,"Coordinate increment per pixel",&status);
	  fits_write_key(spectrum,TINT,"NSCANS",&spapix,"Number of scans added",&status);
	  fits_write_key(spectrum,TSTRING,"IMAGE",imagefile,"Image where the spectrum was extracted from",&status);
	  fits_write_key(spectrum,TSTRING,"RESPONSE",respfile,"Fits file used as response function",&status);
	  fits_write_key(spectrum,TFLOAT,"NOBJ",iobj+i,"Number of object in buho file",&status);
	  /* 	else fits_write_key(spectrum,TDOUBLE,"USNOID",num+i,"USNO identification name",&status); */
	  /* 	fits_write_key(spectrum,TFLOAT,"NOBJ",iobj+i,"Numbre of object in buho file",&status); */
	  /* 	fits_write_key(spectrum,TSTRING,"DISPFILE",dispfile,"Dispersion law used to wavelength calibration",&status); */
	  fits_close_file(spectrum,&status);
	  if(status) fits_report_error(stderr,status);
	}
	if(errflag ) {
	  if(ismultspec) {
	    firstpixel[1]=itrueobj;lastpixel[1]=itrueobj;
	    fits_write_subset_flt(wholeerrflux,1,2,wholenaxes,firstpixel,lastpixel,errfluxspec,&status);
	    sprintf(key,"N%07d",itrueobj);
	    fits_write_key(wholeerrflux,TFLOAT,key,iobj+i,"Number of object in catalogue",&status);
	    sprintf(key,"X%07d",itrueobj);
	    fits_write_key(wholeerrflux,TFLOAT,key,xp+i,"X position of the object in the image",&status);
	    sprintf(key,"Y%07d",itrueobj);
	    fits_write_key(wholeerrflux,TFLOAT,key,yp+i,"Y position of the object in the image",&status);
	    if(isradec) {
	      ra2str(snul1,32,alfa[i],3);
	      dec2str(snul2,32,delta[i],2);
	      sprintf(key,"R%07d",itrueobj);
	      fits_write_key(wholeerrflux,TSTRING,key,snul1,"R.A. of the extracted object",&status);
	      sprintf(key,"D%07d",itrueobj);
	      fits_write_key(wholeerrflux,TSTRING,key,snul2,"Declination of the extracted object",&status);
	    }
	  }
	  else {
	    
	    strcpy(errfluxfile,errfcspec);
	    /* 	  if(fromcat) sprintf(snul,"%12.8f",num[i]); */
	    sprintf(snul,"%08d",(int)iobj[i]);
	    strcat(errfluxfile,snul);
	    strcat(errfluxfile,".fits");
	    ffinit(&errspectrum,errfluxfile,&status);
	    if(status) fits_report_error(stderr,status);
	    fits_create_img(errspectrum, -32,1,naxes,&status);
	    if(status) fits_report_error(stderr,status);
	    fits_write_img(errspectrum,TFLOAT,1,specpix*2,errfluxspec,&status);
	    if(status) fits_report_error(stderr,status);
	    fits_write_key(errspectrum,TFLOAT,"XCPIX",xp+i,"X position of the object in the image",&status);
	    fits_write_key(errspectrum,TFLOAT,"YCPIX",yp+i,"X position of the object in the image",&status);
	    fits_write_key(errspectrum,TFLOAT,"SKY",sky+i,"Substracted sky near the object",&status);
	    /* 	    fits_write_key(errspectrum,TFLOAT,"SKY_SIG",sig+i,"Sigma of sky near the object",&status); */
	    if(isradec) {
	      ra2str(snul1,32,alfa[i],3);
	      dec2str(snul2,32,delta[i],2);
	      fits_write_key(errspectrum,TSTRING,"RA" ,snul1,"R.A. of the extracted object",&status);
	      fits_write_key(errspectrum,TSTRING,"DEC",snul2,"Declination of the extracted object",&status);
	    }
	    fits_write_key(errspectrum,TFLOAT,"CRPIX1",&crpix,"Reference pixel",&status);
	    fits_write_key(errspectrum,TFLOAT,"CRVAL1",&ldomin,"Coordinate at reference pixel",&status);
	    fits_write_key(errspectrum,TFLOAT,"CDELT1",&deltaldo,"Coordinate increment per pixel",&status);
	    fits_write_key(errspectrum,TINT,"NSCANS",&spapix,"Number of scans added",&status);
	    fits_write_key(errspectrum,TSTRING,"IMAGE",imagefile,"Image where the spectrum was extracted from",&status);
	    fits_write_key(errspectrum,TSTRING,"RESPONSE",respfile,"Fits file used as response function",&status);
	    fits_write_key(errspectrum,TFLOAT,"NOBJ",iobj+i,"Number of object in buho file",&status);
	    /* 	  else fits_write_key(errspectrum,TDOUBLE,"USNOID",num+i,"USNO identification name",&status); */
	    fits_close_file(errspectrum,&status);
	  if(status) fits_report_error(stderr,status);
	  }
	}
      }
      if(ismultspec) {
	fprintf(of," %40s[*,%07d:%07d] ",wholespecfile,itrueobj,itrueobj);
	if(waveflag) fprintf(of," %40s[*,%07d:%07d]  ",wholewavefile,itrueobj,itrueobj);
	if(fluxflag) fprintf(of,"  %40s[*,%07d:%07d] ",wholefluxfile,itrueobj,itrueobj);
	if(errflag) fprintf(of," %40s[*,%07d:%07d] ",wholeerrspecfile,itrueobj,itrueobj);
	if(errflag && waveflag) fprintf(of," %40s[*,%07d:%07d] ",wholeerrwavefile,itrueobj,itrueobj);
	if(errflag && fluxflag) fprintf(of," %51s[*,%07d:%07d] ",wholeerrfluxfile,itrueobj,itrueobj);
      }
      else {
	fprintf(of," %59s ",specfile);
	if(waveflag) fprintf(of," %59s  ",wavefile);
	if(fluxflag) fprintf(of,"  %59s ",fluxfile);
	if(errflag) fprintf(of," %59s ",errspecfile);
	if(errflag && waveflag) fprintf(of," %59s ",errwavefile);
	if(errflag && fluxflag) fprintf(of," %59s ",errfluxfile);
      }
	fprintf(of,"  %12.2f  %12.2f ",xp[i],yp[i]);
	if(isradec)   {  
	  ra2str(snul1,32,alfa[i],3);
	  dec2str(snul2,32,delta[i],2);
	  fprintf(of," %32s %32s ",snul1,snul2);
      }
      fprintf(of,"  %10.2f \n",sky[i]);
    }
  }
  
  printf("\n Finishing extraction\n");

  fclose(of);
  if(ismultspec) {
    wholenaxes[1]=itrueobj;
    fits_modify_key_lng(wholespec,"NSPEC",itrueobj,"&",&status);
    fits_close_file(wholespec,&status);
    ffopen(&wholespec,wholespecfile,READWRITE,&status);
    fits_resize_img(wholespec,-32,2,wholenaxes,&status); 
    fits_close_file(wholespec,&status);

    if(errflag) {
      fits_modify_key_lng(wholeerrspec,"NSPEC",itrueobj,"&",&status);
      fits_close_file(wholeerrspec,&status);
      ffopen(&wholeerrspec,wholeerrspecfile,READWRITE,&status);
      fits_resize_img(wholeerrspec,-32,2,wholenaxes,&status);
      fits_modify_key_lng(wholeerrspec,"NSPEC",itrueobj,"&",&status);
      fits_close_file(wholeerrspec,&status);
    }
      
    if(waveflag) {
      fits_modify_key_lng(wholewave,"NSPEC",itrueobj,"&",&status);
      fits_close_file(wholewave,&status);
      ffopen(&wholewave,wholewavefile,READWRITE,&status);
      fits_resize_img(wholewave,-32,2,wholenaxes,&status);
      fits_modify_key_lng(wholewave,"NSPEC",itrueobj,"&",&status);
      fits_close_file(wholewave,&status);
      
      if(errflag) {
	fits_modify_key_lng(wholeerrwave,"NSPEC",itrueobj,"&",&status);
	fits_close_file(wholeerrwave,&status);
	ffopen(&wholeerrwave,wholeerrwavefile,READWRITE,&status);
	fits_resize_img(wholeerrwave,-32,2,wholenaxes,&status);
	fits_modify_key_lng(wholeerrwave,"NSPEC",itrueobj,"&",&status);
	fits_close_file(wholeerrwave,&status);
      }
    }
    if(fluxflag) {
      fits_modify_key_lng(wholeflux,"NSPEC",itrueobj,"&",&status);
      fits_close_file(wholeflux,&status);
      ffopen(&wholeflux,wholefluxfile,READWRITE,&status);
      fits_resize_img(wholeflux,-32,2,wholenaxes,&status);
      fits_modify_key_lng(wholeflux,"NSPEC",itrueobj,"&",&status);
      fits_close_file(wholeflux,&status);
      
      if(errflag) {
	fits_modify_key_lng(wholeerrflux,"NSPEC",itrueobj,"&",&status);
	fits_close_file(wholeerrflux,&status);
	ffopen(&wholeerrflux,wholeerrfluxfile,READWRITE,&status);
	fits_resize_img(wholeerrflux,-32,2,wholenaxes,&status);
	fits_modify_key_lng(wholeerrflux,"NSPEC",itrueobj,"&",&status);
	fits_close_file(wholeerrflux,&status);
      }
    }
  }
  sprintf(touchchar,"/bin/touch %s.extractdone\n",imagefile);
  system(touchchar);
  
  return(0);
}











void ReadResp(char respfile[],  int fluxcal)
{
  fitsfile *respo;
  int status=0;
  int nfound, anynull;
  long r_naxes[1];
  long r_pixels;
  int i;
  long  fpixel;
  float  nullval;
  char comment[100];

  fluxcal=0;
/*   //printf(" specpix vale %d\n",specpix); */
  if(strcmp(respfile,"NONE")) {
    printf(" HORORORORORO\n");
    fluxcal=1;
    if( ffopen(&respo, respfile, READONLY, &status)) fits_report_error(stderr,status);
    if(fits_read_keys_lng(respo, "NAXIS", 1, 1, r_naxes, &nfound, &status)) fits_report_error(stderr,status);
    if(nfound!=1) {
      printf(" FITS file does not contain 1 dimension\n");
      exit(1);
    }
    
    r_pixels=r_naxes[0];
    if((respuesta=malloc(r_pixels*sizeof(float)))==NULL) printf("I cannot dimension respuesta of %ld elements \n",r_pixels);
    if((ldoresp=malloc(r_pixels*sizeof(float)))==NULL) printf("I cannot dimension ldoresp of %ld elements \n",r_pixels);
    printf(" Dimesiaondan %ld\n",r_pixels);
    fpixel=1;
    nullval=0;
    printf("...Reading file %s \n",respfile);
    if(fits_read_img(respo, TFLOAT, fpixel, r_pixels, &nullval, respuesta, &anynull, &status )) fits_report_error(stderr,status);
    if(fits_read_key(respo, TFLOAT, "CRPIX1", &r_pix1, comment, &status )) fits_report_error(stderr,status);
    if(fits_read_key(respo, TFLOAT, "CRVAL1", &r_val1, comment, &status )) fits_report_error(stderr,status);
    if(fits_read_key(respo, TFLOAT, "CDELT1", &r_delt1, comment, &status )) fits_report_error(stderr,status);
/*     //    printf(" He liefo  %f %f %f\n",r_pix1,r_val1,r_delt1);     */
/*     //exit(1); */
    for(i=0;i<r_pixels;i++) {
      ldoresp[i]=r_val1+r_delt1*(i+1-r_pix1);
    }

    /* Termino de leerlo */
  }
  else {
    r_pixels=specpix*2;
    if((respuesta=malloc(r_pixels*sizeof(float)))==NULL) printf("I cannot dimension respuesta of %ld elements \n",r_pixels);
    if((ldoresp=malloc(r_pixels*sizeof(float)))==NULL) printf("I cannot dimension ldoresp of %ld elements \n",r_pixels);

    for(i=0;i<specpix*2;i++) {
      respuesta[i]=1.;
      ldoresp[i]=ldomin+i*deltaldo;
/*       //printf("NONE %d la respuesta %f ldo %f\n",i,respuesta[i],ldoresp[i]); */
    }


  }
  nresp=r_pixels;

/*   //printf(" 1er ldo %f\n",ldoresp[0]); */
/*   for(i=0;i<r_pixels;i++) { */
/*     printf(" %d ldoresp %f rep %f %f\n",i,ldoresp[i],respuesta[i],resp(ldoresp[i])); */
/*   } */
  
  

}



void ReadDisp(char dispfile[])
{
  FILE *fp;
  int i;    
/*   float dxdl; */
/*   int ipix,ipix2; */
/*   float ldo; */
/*   int j; */
  ndisp=FileNLin(dispfile);
/*   //  printf("Ndisp es %d\n",ndisp); */
  if((fp=fopen(dispfile,"r"))==NULL) {
    printf("Cannot open file %s\n",dispfile);
    exit(1);
  }
  if((ldodisp=malloc(ndisp*sizeof(float)))== NULL) {
    printf("I can't dimension the vector ldodisp of %d elements",ndisp);
    exit(1);
  }
  if((position=malloc(ndisp*sizeof(float)))== NULL) {
    printf("I can't dimension the Vector position of %d elements",ndisp);
    exit(1);
  }
  for (i=0;i<ndisp;i++) {
    fscanf(fp," %f %f",ldodisp+i,position+i);
/*     //printf(" bonito %f %f \n",ldodisp[i],position[i]); */
  }
/*   //cpgopen("?"); */
/*   //cpgswin(3400,6000,0.,0.3); */
/*   //cpgmove(3400.,0.); */
/*   //cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0); */
  /*
    for(j=0;j<specpix*2;j++) {
    ldo=3487.17+j*12.98;
    ipix=(int)(ldo2x(ldo));
    //    ipix2=(int)(ldo2x(ldo+4*12.98));
    //    dxdl= (ipix2-ipix)/4./12.98;
    dxdl=139.e5/((ldo-1606)*(ldo-1606))/15;
    
    //dxdl: derivada de x con respecto a lambda para transformar.
    //    if(ipix<0) ipix=0;
    //    if(ipix>(specpix*2-1)) ipix=specpix*2-1;
    cpgdraw(ldo,dxdl);
    printf(" ldo %f pix %d dxdl %f   A %f\n",ldo,ipix,dxdl,a);
    }
    cpgswin(3400,6000,150.,200);
    cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
    cpgmove(3400.,0.);
    
    for(j=0;j<specpix*2;j++) {
    ldo=3487.17+j*12.98;
    ipix=(int)(ldo2x(ldo));
    cpgdraw(ldo,ipix);
    }
    
    
    
    exit(1);
  */
/*   //cpgend(); */
  

}


float ldo2x(float ldo)
{
/*   //  return(Lagr4(ldodisp,position,ndisp,ldo)); */
  return((A+B/(ldo-C))/15.);

}

float x2ldo(float x)
{
/*   //  return(Lagr4(position,ldodisp,ndisp,x)); */
  return(C-B/(A-x*15.));
}

float dxdl(float ldo)
{
  return(-B/((ldo-C)*(ldo-C))/15.);
}


float resp(float ldo) {
/*   //printf("nresp vale %d\n",nresp); */
  return(Lagr2(ldoresp,respuesta,nresp,ldo));
}




void ReadOptions() {
  dispaxis=0;
  printf(" Input FITS file: ");
  reads(imagefile,imagefile);
  printf(" Input error image for previous file (NONE=no error treatment): ");
  reads(errimagefile,errimagefile);
  if(strcmp(errimagefile,"NONE"))    errflag=1;
  printf(" Input catalogue with objects: ");
  reads(catfile,catfile);
  printf(" Input column in %s with x positions (negative for RA positions): ",catfile);
  xcol=readi(2);
  if(xcol<0) {
    isradec=1;
    xcol=-xcol;
  }
  if (isradec)   printf(" Input column in %s with DEC positions: ",catfile);
  else  printf(" Input column in %s with y positions: ",catfile);
  ycol=readi(3);
  printf(" Root name for espectra: ");
  reads(rawspec,rawspec);
  if(errflag) {
    printf(" Root name for error espectra: ");
    reads(errrawspec,errrawspec);
  }

  printf(" Half number of pixels in spatial direction to extract: ");
  spapix=readi(5);
  printf(" Half of pixels in spectrum direction to extract: ");
  specpix=readi(40);
  while(dispaxis!=1 && dispaxis!=2) {
    printf(" Dispresion axis (1: X , 2: Y): ");
    dispaxis=readi(1);
  }
  printf(" File with dispersion law (Angstroms   Delta(pixel)): ");
  reads(dispfile,dispfile);
  printf(" Constat A: ");
  A=readf(1.);
  printf(" Constat B: ");
  B=readf(1.);
  printf(" Constat C: ");
  C=readf(1.);
/*   //ReadDisp(dispfile); */

  printf(" FITS file with response function of the system (NONE for no response function: ");
  reads(respfile,respfile);
  
  printf(" Root name for wavelength calibrated espectra: ");
  reads(wcspec,wcspec);
  printf(" Root name for flux & wavelength calibrated espectra: ");
  reads(fcspec,fcspec);
  if(errflag) {
    printf(" Root name for wavelength calibrated error espectra: ");
    reads(errwcspec,errwcspec);
    printf(" Root name for flux & wavelength calibrated error espectra: ");
    reads(errfcspec,errfcspec);
  }
/*   printf(" Root name for boxes extracted around spectra: "); */
/*   scanf("%s",objbox); */

}


void LoadOptions(char filepar[50]) {
  fitsfile *parfile;
  char key[9]="";
  int status=0;
/*   //char string[51]; */
  char comment[51];


  printf("Using parameter file <<%s>>\n",filepar);
  if( ffopen2(&parfile,filepar, READONLY, &status)) fits_report_error(stderr,status);

  sprintf(key,"IMAGE");
  ffgky(parfile,TSTRING,key,imagefile,comment,&status);
  fits_report_error(stderr,status);
/*   //  status=0; */
  sprintf(key,"ERRIMAGE");
  ffgky(parfile,TSTRING,key,errimagefile,comment,&status);
  fits_report_error(stderr,status);
  if(strcmp(errimagefile,"NONE"))    errflag=1;
      
  sprintf(key,"CATFILE");
  ffgky(parfile,TSTRING,key,catfile,comment,&status);
  fits_report_error(stderr,status);
/*   //  status=0; */
/*   sprintf(key,"ERRCAT"); */
/*   ffgky(parfile,TSTRING,key,errcatfile,comment,&status); */
/*   fits_report_error(stderr,status); */
/*   //  status=0; */
  sprintf(key,"XCOL");
  ffgky(parfile,TINT,key,&xcol,comment,&status);  
  fits_report_error(stderr,status);
/*   //  status=0; */
  sprintf(key,"YCOL");
  ffgky(parfile,TINT,key,&ycol,comment,&status);
  fits_report_error(stderr,status);
/*   //  status=0; */
/*   sprintf(key,"XDIRCOL"); */
/*   ffgky(parfile,TINT,key,&xdircol,comment,&status);   */
/*   fits_report_error(stderr,status); */
/*   //  status=0; */
/*   sprintf(key,"YDIRCOL"); */
/*   ffgky(parfile,TINT,key,&ydircol,comment,&status);   */
/*   fits_report_error(stderr,status); */
/*   //  status=0; */
  sprintf(key,"SKYFILE");
  ffgky(parfile,TSTRING,key,fitssky,comment,&status);  
  fits_report_error(stderr,status);
/*   //  status=0; */
  sprintf(key,"ERRSKFIL");
  ffgky(parfile,TSTRING,key,fitserrsky,comment,&status);  
  fits_report_error(stderr,status);
/*   //  status=0; */
  sprintf(key,"SPAPIX");
  ffgky(parfile,TINT,key,&spapix,comment,&status);  
  fits_report_error(stderr,status);
/*   //  status=0; */
  sprintf(key,"SPECPIX");
  ffgky(parfile,TINT,key,&specpix,comment,&status);  
  fits_report_error(stderr,status);
/*   //  status=0; */
  sprintf(key,"DISPAXIS");
  ffgky(parfile,TINT,key,&dispaxis,comment,&status);
  fits_report_error(stderr,status);
/*   //  status=0; */
  sprintf(key,"RADEC-XY");
  ffgky(parfile,TLOGICAL,key,&isradec,comment,&status);
  fits_report_error(stderr,status);

  sprintf(key,"MULTSPEC");
  ffgky(parfile,TLOGICAL,key,&ismultspec,comment,&status);
  fits_report_error(stderr,status);
/*   //  status=0; */
  sprintf(key,"RAWSPEC");
  ffgky(parfile,TSTRING,key,rawspec,comment,&status);
  fits_report_error(stderr,status);
/*   //  status=0; */
  sprintf(key,"WCSPEC");
  ffgky(parfile,TSTRING,key,wcspec,comment,&status);
  fits_report_error(stderr,status);
/*   //  status=0; */
  sprintf(key,"FCSPEC");
  ffgky(parfile,TSTRING,key,fcspec,comment,&status);
  fits_report_error(stderr,status);
/*   //  status=0; */
   sprintf(key,"ERRRAWSP");
  ffgky(parfile,TSTRING,key,errrawspec,comment,&status);
  fits_report_error(stderr,status);
/*   //  status=0; */
  sprintf(key,"ERRWCSP");
  ffgky(parfile,TSTRING,key,errwcspec,comment,&status);
  fits_report_error(stderr,status);
/*   //  status=0; */
  sprintf(key,"ERRFCSP");
  ffgky(parfile,TSTRING,key,errfcspec,comment,&status);
  fits_report_error(stderr,status);
/*   sprintf(key,"OBJBOX"); */
/*   ffgky(parfile,TSTRING,key,objbox,comment,&status); */
/*   fits_report_error(stderr,status); */
/*   //  status=0; */
  sprintf(key,"A_DISP");
  ffgky(parfile,TFLOAT,key,&A,comment,&status);
  fits_report_error(stderr,status);
/*   //  status=0; */
  sprintf(key,"B_DISP");
  ffgky(parfile,TFLOAT,key,&B,comment,&status);
  fits_report_error(stderr,status);
/*   //  status=0; */
  sprintf(key,"C_DISP");
  ffgky(parfile,TFLOAT,key,&C,comment,&status);
  fits_report_error(stderr,status);
/*   //  status=0; */
  sprintf(key,"RESPFILE");
  ffgky(parfile,TSTRING,key,respfile,comment,&status);
  fits_report_error(stderr,status);
/*   //  status=0; */
  ffgky(parfile,TSTRING,"OUTFILE",outfile,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"DEVICE",pgdevice,comment,&status);
  fits_report_error(stderr,status);

  if(status!=0) {
    fits_report_error(stderr,status);
    printf("ERROR: Reading parameter file. Not enough keywords.\n");
    exit(1);
  }

  if(strcmp(fcspec,"NONE"))  fluxflag=1;
  if(strcmp(respfile,"NONE")) fluxflag=1;  
  if(strcmp(wcspec,"NONE"))  {waveflag=1;}
  else {waveflag=0;fluxflag=0;}

}


