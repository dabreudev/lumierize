/* #undef __cplusplus */
#include "modulos.h"
#define NMAX     15000
#define NOBJDET  25000
#define NSECMAX  500
#define NCROSSMIN 40
#define WCSVERSION28 
#define DEBUG 0
#define MAXITERMATCH 20
/* #define WCSVERSION26 */


/* //A lo mejor hay que cambiarlo para que ajuste el de grado 2 */


#define  FTOL    0.0001 
struct obj_table {
  float xpos;
  float ypos;
  float flux;
  int log;
};
struct star_table {
  double ar;
  double dec;
  double num;
  double mag; 
  double mag2;
};

char parfilename[100]="";
     
int incompare(const void *X1,const void *X2);
int iscompare(const void *X1,const void *X2); 
int min(int x1,int x2);
void SaveParam(char parfilename[100]);
void LoadParam(char file[100]);
void WriteWCS_keys();
int  ComputeSolution(double *alfa,double *delta,
                     double *xpix,double *ypix);
void ShiftSolution(double *alfas, double *deltas,
                   double *xpixs, double *ypixs);
void InputInitGuess(int *initguess);
void Keyboard(double *xpoint, double *ypoint, double *alfa, double *delta, int *nstar);
void Buho(double *xpoint, double *ypoint,  int *npoint,float flux1,float flux2);
void ClickImage(double *xpoint, double *ypoint, double *alfa, double *delta,  int *nstar);
void File(double *xpoint, double *ypoint, double *alfa, double *delta, int *nstar);
void FindGSC(double *ar,double *dec,double *mag,float alfacen,float deltacen,float xdim,float ydim,int *nstar,double mag1,double mag2,float factx,float facty,int nmax);
void FindAPM(double *ar,double *dec,double *mag,float alfacen,float deltacen,float xdim,float ydim,int *nstar,double mag1,double mag2,float factx,float facty,int nmax);
void FindACT(double *ar,double *dec,double *mag,float alfacen,float deltacen,float xdim,float ydim,int *nstar,double mag1,double mag2,float factx,float facty,int nmax);
void FindUSNO(double *ar,double *dec,double *mag,float alfacen,float deltacen,float xdim,float ydim,int *nstar,double mag1,double mag2,float factx,float facty,int nmax);
void FindTycho2(double *ar,double *dec,double *mag,float alfacen,float deltacen,float xdim,float ydim,int *nstar,double mag1,double mag2,float factx,float facty,int nmax);
void Plot(float *x,float *y,int n,int mark,int color);
void Plotd(double *x,double *y,int n,int mark,int color);
void PlotStars(double *ar,double *dec,double *mag,int nstar,double *alfa,double *delta);
void ClickBoth(double *ar    ,double *dec   ,int nstar,
               double *xpoint,double *ypoint,int npoint,
               double *alfa  ,double *delta ,
               double *xpix,  double *ypix  ,
               int naxes0,   int naxes1);
void ClickImFindStar(double *ar    ,double *dec   ,int nstar,
                     double *alfa  ,double *delta ,
                     double *xpix,  double *ypix  ,
                     int naxes0,   int naxes1);
void DeleteClickIm(double *alfa  ,double *delta ,
		   double *xpix,  double *ypix  );
void DeleteClickErr(double *alfa  ,double *delta ,
		    double *xpix,  double *ypix  ,
                    float *xerror, float *yerror);
void InitWCS(char *inputfile,int *initguess);

void StarMatch_ceg(int *ncross, double *alfa, double *delta, double *xpix, double *ypix, int *off, int ns,double  *sx, double *sy, int ng, double *gra, double *gdec, double *gx, double *gy, double tol, struct WorldCoor *wcs, int debug);
void MatchTol(int *ncross, double *alfa, double *delta, double *xpix, double *ypix, int *off, int ns,double  *sx, double *sy, int ng, double *gra, double *gdec, double *gx, double *gy, double tol, struct WorldCoor *wcs, int debug);



/* void PlotPointRef(float *xpoint,float *ypoint,int npoint); */
/* void PlotCross(float *,float *, ...); */

/* struct WorldCoor *GetWCSFITS(char *filename); */


/* Varaibales para leer el FITS */
char inputfile[51]="";
int status=0;
int nfound, anynull;
fitsfile *fitsimage;
long naxes[2];
long fpixel, nbuffer, npixels, ii;
float datamin, datamax, nullval;
float *buffer;
int nkey;

char *header;

/* //Variables para el fichero BUHO */
char bfile[51];
FILE *buhofile;
int bxcol=14;
int bycol=15;
int bfcol=4;

/* //Variables de la astrometria que deberian ver todas las subrutinas */
double epoch=2000.; 
float rotang=0;
int flip=0; 
struct WorldCoor *wcsim=0;      /* World coordinate system structure */
int errflag=0;
int ncross=0;
int ncrosssec=0;
float alfac=0, deltac=0; 
float xpixsiz=0,ypixsiz=0; /*  En arcsec/pixel */ 
int secflag=0;    /* Flag para saber si hay estrellas secundarias o no */
double shiftra=0.,shiftdec=0.;
float xem,yem,aem,dem;
float xes,yes,aes,des;



/* //Parametros del fichero de parametros */

char astcat[10];
char seccat[10];
int ngrad=3;
int nbrightestobj=50;
int nbrighteststar=50;
char pgdevice[100];
int interact=1;


int main(int argc, char **argv) 
{
  double dnul1,dnul2;
  char snul[1000];
  char option='C';
/*   //Coordenadas de estrellas de referencia: */
  double ar[NMAX],dec[NMAX];
  double xstar[NMAX],ystar[NMAX];     /* //Calculadas. De momento no hacen falta */
  int off[NMAX];
  int nstar=0;
  int nstar_m=0;
/*   //Coordenadas de estrellas secundarias de referencia: */
  double ars[NSECMAX],decs[NSECMAX];
  double xstars[NSECMAX],ystars[NSECMAX];     /* //Calculadas. De momento no hacen falta */
  int offs[NSECMAX];
  double xpixs[NSECMAX],ypixs[NSECMAX];         /* //En la imagen  */
  double alfas[NSECMAX],deltas[NSECMAX];        /* //Coordenadas reales. */
  double mags[NSECMAX];
  int nstarsec=0;
/*   //Coordenadas de puntos de referencia */
  double xpoint[NOBJDET],ypoint[NOBJDET];
  int npoint=0;
  int npoint_m=0;
  int nitermatch=0;
  int ngradtmp;
  int firstwcs=0;
/*   //Correlacion entre ambas: */
/*   //  float xpix[NMAX],ypix[NMAX]; */
/*   //float alfa[NMAX],delta[NMAX]; */
  float tol=20.;
  float dxys,dxy;
  int jgs;
  double xpix[NMAX],ypix[NMAX];         /* //En la imagen  */
  double alfa[NMAX],delta[NMAX];        /* //Coordenadas reales. */
  double mag[NMAX];
/*   //Calculo de las posiciones segun la actual transformacion */
/*   //double xpixcal[NMAX],ypixcal[NMAX];   //De alfa,delta a xpixcal,ypixcal */
/*   //double alfacal[NMAX],deltacal[NMAX];  //De xpix,ypix a alfacal, deltacal */
/*   //float xerror[NMAX],yerror[NMAX],rerror[NMAX]; //Diferencia xpix-xpixcal... */
/*   //float aerror[NMAX],derror[NMAX];              //Diferencia alfa-alfacal... */
/*   float xcal,ycal; */
/*   float psi,eta; */
  int i,j;
  float mean,sigma;
  /* Luego lo cambiare si hace falta */
  float flux1=0,flux2=0;
  double mag1=0,mag2=0;
  float factra=1,factdec=1;

  int initguess=0;
  int fitwcs=0;
/*   struct plt_cte ctes; */

  int verbose;

  char touchchar[500];

  float tr[6];
  FILE *filepar;

  /* Variables del dibujo */


/*   //int daschan;  //Esto es para la WFC */

  if(argc < 2) {
    strcpy(pgdevice,"?");
    /* Entrada de datos */
    printf("\n Input FITS image where to compute astrometry: ");
    reads(inputfile,inputfile);  
/*     scanf("%s",inputfile);   */
  }
  else {
    if((filepar=fopen(argv[1],"r"))==NULL) {
      printf("ERROR: Can't open options file %s\n",argv[1]);
      exit(1);
    }
    strcpy(parfilename,argv[1]);
    LoadParam(parfilename);
    interact=0;
/*     //initguess=1; */
  }

  verbose=interact;

  cpgbeg(0,pgdevice,1,1);

  /* Leo la imagen FITS de entrada */
  if( ffopen(&fitsimage, inputfile, READWRITE, &status)) fits_report_error(stderr,status);
  if(fits_read_keys_lng(fitsimage, "NAXIS", 1, 2, naxes, &nfound, &status)) fits_report_error(stderr,status);
  npixels=naxes[0]*naxes[1];
  /* Ahora compruebo que no hay ninguna 
     solucion ya establecida  antes de cargar la imagen*/
  if(interact) InitWCS(inputfile,&initguess);
  else initguess=0;

  if((buffer=malloc(npixels*sizeof(float)))==NULL) printf("I cannot dimension buffer of %ld elements \n",npixels);
  fpixel=1;
  nullval=0;
  datamin=1.0e30;
  datamax=-1.0e30;
  printf("...Reading image %s \n",inputfile);
  if(fits_read_img(fitsimage, TFLOAT, fpixel, npixels, &nullval, buffer, &anynull, &status )) fits_report_error(stderr,status);
/*   //printf("...Computing datamin and datamax \n"); */
/*   //Cierro el fichero con FITSIO: */
/*   //if(ffclos(fitsimage,&status)) fits_report_error(stderr,status); */



  mean=StMedia(npixels,buffer,&sigma);
  datamin=mean-sigma*3;
  datamax=mean+sigma*3;
  printf(" Datamin cut %f Datamax cut %f \n",datamin,datamax);

/*   for (ii=0 ;ii<npixels;ii++) { */
/*     //printf(" buffer %d %e \n",ii,buffer[ii]); */
/*     if( buffer[ii]< datamin) datamin=buffer[ii]; */
/*     if(buffer[ii]> datamax) datamax=buffer[ii]; */
/*   } */
/*   printf(" Datamin %f Datamax %f \n",datamin,datamax); */
  
  /* Dibujo la imagen  */




  fflush(NULL);
  setvbuf(stdin,"",_IOLBF,0);
  setvbuf(stdout,"",_IOLBF,0);


  if(interact) {
    cpgsvp(0.1,0.5,0.1,0.9);
    cpgwnad(0.,naxes[0],0.,naxes[1]);
    tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
    cpggray(buffer,naxes[0],naxes[1],1,naxes[0],1,naxes[1],datamax,datamin,tr);
    cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
    strcpy(snul,"Image \0");
    strcat(snul,inputfile);
    cpglab("X axis","Y axis",snul);
  }
  

  if(!interact) {
    printf(" Non interactive mode.\n");
    if(!initguess) {
      printf(" Setting initial coordinates AR %f Dec %f ROT %f\n",alfac/15,(float)deltac,rotang);
      wcsim=wcsxinit((double)(alfac),(double)(deltac),(double)(xpixsiz),naxes[0]/2.,naxes[1]/2.,naxes[0],naxes[1],0.,2000,epoch,"TAN");
      if(!flip) wcsdeltset(wcsim, wcsim->cdelt[0],wcsim->cdelt[1],(double)rotang);
      else      wcsdeltset(wcsim,-wcsim->cdelt[0],wcsim->cdelt[1],(double)rotang);
/*       pix2wcs(wcsim,naxes[0]/1., naxes[1]/1. ,&dnul1, &dnul2); */
/*       printf(" 2000 2000 %f %f\n",dnul1,dnul2); */
/*       pix2wcs(wcsim,0., naxes[1]/1. ,&dnul1, &dnul2); */
/*       printf(" 0 2000 %f %f\n",dnul1,dnul2); */
/*       pix2wcs(wcsim,0., 0. ,&dnul1, &dnul2); */
/*       printf(" 0 0 %f %f\n",dnul1,dnul2); */
/*       exit(1); */



      /*       wcsdeltset(wcsim,wcsim->cdelt[0],wcsim->cdelt[1],90.);  */
    }
    /*     SetFITSWCS(header,wcsim); */
    /*     printf(" Pasa de aqui\n"); */
    /*     PrintWCS(header,1); */
    /*     exit(1); */
    
    printf(" Reading sources file\n");

    Buho(xpoint,ypoint,&npoint,0.,0.);
    if(npoint!=0) Plotd(xpoint,ypoint,npoint,4,4); 
    printf(" Using %d sources \n",npoint);
    
    
    if(strsrch(astcat,"ACT") != NULL)  
      FindACT(ar,dec,mag, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstar,0.,0.,1.1,1.1,NMAX);
    else if(strsrch(astcat,"USNO") != NULL)
      FindUSNO(ar,dec,mag, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstar,0.,0.,1.1,1.1,NMAX);
    else if(strsrch(astcat,"GSC") != NULL) 
      FindGSC(ar,dec,mag, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstar,0.,0.,1.1,1.1,NMAX);
    else if(strsrch(astcat,"APM") != NULL) 
      FindAPM(ar,dec,mag, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstar,0.,0.,1.1,1.1,NMAX);
    else if(strsrch(astcat,"Tycho-2") != NULL) 
      FindTycho2(ar,dec,mag, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstar,0.,0.,1.1,1.1,NMAX);


    ncross=0;
    if(nstar!=0) PlotStars(ar,dec,mag,nstar,alfa,delta);
    printf(" Using %d reference stars \n",nstar);
/*     //exit(1); */
    
    tol=20.;
    for(i=0;i<nstar;i++) {
      xstar[i]=0;ystar[i]=0;
      wcs2pix(wcsim,ar[i],dec[i],&xstar[i],&ystar[i],&off[i]);
/*       printf(" %d/%d   x %f  y %f RA %f DEC %f\n",i,nstar,xstar[i],ystar[i],ar[i],dec[i]); */
    }
    for(i=0;i<npoint;i++) {
/*       printf(" %d/%d   x %f  y %f \n",i,npoint,xpoint[i],ypoint[i]); */
    }

/*     npoint=min(nbrightestobj,npoint); */
/*     nstar=min(nbrighteststar,nstar); */

    printf(" Verbose %d\n",verbose);
    npoint_m=85;
    nstar_m=20; 
    nitermatch=0;
    firstwcs=0;
    while(ncross<NCROSSMIN && nitermatch<MAXITERMATCH ) {
      if(nitermatch==1) nstar_m*=2;
      if(nitermatch==2) npoint_m*=1.3;
      if(nitermatch==3) nstar_m*=2;
      if(nitermatch==4) npoint_m*=1.3;
      if(nitermatch==5) nstar_m*=1.6;
      if(nitermatch==5) {
	if(strsrch(astcat,"ACT") != NULL)  
	  FindACT(ar,dec,mag, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstar,0.,0.,1.0,1.0,NMAX);
	else if(strsrch(astcat,"USNO") != NULL)
	  FindUSNO(ar,dec,mag, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstar,0.,0.,1.0,1.0,NMAX);
	else if(strsrch(astcat,"GSC") != NULL) 
	  FindGSC(ar,dec,mag, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstar,0.,0.,1.0,1.0,NMAX);
	else if(strsrch(astcat,"APM") != NULL) 
	  FindAPM(ar,dec,mag, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstar,0.,0.,1.0,1.0,NMAX);
	else if(strsrch(astcat,"Tycho-2") != NULL) 
	  FindTycho2(ar,dec,mag, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstar,0.,0.,1.0,1.0,NMAX);
	nstar_m*=1.2;
	npoint_m*=1.2;
      }
      if(nitermatch>5) {
	nstar_m*=1.2;
	npoint_m*=1.2;
      }

      pix2wcs(wcsim,naxes[0]/2., naxes[1]/2. ,&dnul1, &dnul2);
      alfac=(float)dnul1;deltac=(float)dnul2;

      if(nstar_m>nstar) nstar_m=nstar;
      if(npoint_m>npoint) npoint_m=npoint;



      StarMatch_ceg(&ncross, alfa, delta, xpix, ypix, off, npoint_m,xpoint,ypoint,nstar_m,ar,dec,xstar,ystar,tol,wcsim,0);
    cpgsvp(0.1,0.5,0.1,0.9);
    cpgwnad(0.,naxes[0],0.,naxes[1]);
/*     tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1; */
/*     cpggray(buffer,naxes[0],naxes[1],1,naxes[0],1,naxes[1],datamax,datamin,tr); */
    cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
    strcpy(snul,"Image \0");
    strcat(snul,inputfile);
    cpglab("X axis","Y axis",snul);

      Plotd(xpix,ypix,ncross,2,2);
      Plotd(xpoint,ypoint,npoint_m,3,3);
      Plotd(xstar,ystar,nstar_m,4,4); 
/*       cpgclos(); */
/*       exit(1); */
/* 	printf(" Sale con ncross %d\n",ncross); */
/* 	printf(" WCS after the matching\n"); */
/* 	SetFITSWCS(header,wcsim); */
/* 	PrintWCS(header,1); */

/* 	for(i=0;i<npoint;i++) printf(" %d xpo %f ypo %f\n",i,xpoint[i],ypoint[i]); */
/* 	for(i=0;i<nstar;i++) printf(" %d ra %f dec %f xst %f yst %f\n",i,ar[i],dec[i],xstar[i],ystar[i]); */

/* 	printf(" Todo con ceg\n"); */
/* 	printf(" Number of matchings %d\n",ncross); */

/* 	exit(1); */
      if(ncross>15 && firstwcs==20) { /* De momento esto no lo hago */
	ngradtmp=ngrad;
	ngrad=3;
	printf(" Computing partial solution in the mean time\n");  
	ComputeSolution(alfa,delta,xpix,ypix);
	ngrad=ngradtmp;
	firstwcs=1;
      }
      printf(" Matching stars in iteration %d. %d matches with %d stars and %d objects \n",nitermatch,ncross,nstar_m,npoint_m);
      nitermatch++;
    }
/*     ncross=StarMatch(npoint,xpoint,ypoint,nstar,ar,dec,xstar,ystar,tol,wcsim,0); */
    printf(" At end, using %d matchings\n", ncross);
    
/*     exit(1); */
/*     nm=0; */
/*     for (i = 0; i < npoint; i++) { */
/*       dxys = -1.0; */
/*       igs = -1; */
/*       for (j = 0; j < nstar; j++) { */
/* 	if (!off[j]) { */
/* 	  dxy=(xstar[j]-xpoint[i])*(xstar[j]-xpoint[i])+(ystar[j]-ypoint[i])*(ystar[j]-ypoint[i]); */
/* 	  if (dxy < tol*tol) { */
/* 	    if (dxys < 0.0) { */
/* 	      dxys = dxy; */
/* 	      igs = j; */
/* 	    } */
/* 	    else if (dxy < dxys) { */
/* 	      dxys = dxy; */
/* 	      igs = j; */
/* 	    } */
/* 	  } */
/* 	} */
/*       } */
/*       if (igs > -1) { */
/* 	alfa[nm]=ar[igs]; */
/* 	delta[nm]=dec[igs]; */
/* 	xpix[nm]=xpoint[i]; */
/* 	ypix[nm]=ypoint[i]; */
/* 	nm++; */
/*       } */
/*     } */

/*     printf(" Salio de aqui \n"); */
/*     for(i=0;i<nstar;i++) { */
/*       xstar[i]=0;ystar[i]=0; */
/*       wcs2pix(wcsim,ar[i],dec[i],&xstar[i],&ystar[i],&off[i]); */
/*     } */
/*     ncross=StarMatch(npoint,xpoint,ypoint,nstar,ar,dec,xstar,ystar,tol,wcsim,verbose); */
/*     nm=0; */
/*     for (i = 0; i < npoint; i++) { */
/*       dxys = -1.0; */
/*       igs = -1; */
/*       for (j = 0; j < nstar; j++) { */
/*         if (!off[j]) {     */
/* 	  dxy=(xstar[j]-xpoint[i])*(xstar[j]-xpoint[i])+(ystar[j]-ypoint[i])*(ystar[j]-ypoint[i]); */
/*           if (dxy < tol*tol) { */
/*             if (dxys < 0.0) { */
/*               dxys = dxy; */
/*               igs = j; */
/*             } */
/*             else if (dxy < dxys) { */
/*               dxys = dxy; */
/*               igs = j; */
/*             } */
/*           } */
/*         } */
/*       } */
/*       if (igs > -1) { */
/*         alfa[nm]=ar[igs]; */
/*         delta[nm]=dec[igs]; */
/*         xpix[nm]=xpoint[i]; */
/*         ypix[nm]=ypoint[i]; */
/*         nm++; */
/*       } */
/*     } */
/*     ncross=StarMatch(npoint,xpoint,ypoint,nstar,ar,dec,xstar,ystar,tol,wcsim,verbose); */
/*     nm=0; */
/*     for (i = 0; i < npoint; i++) { */
/*       dxys = -1.0; */
/*       igs = -1; */
/*       for (j = 0; j < nstar; j++) { */
/*         if (!off[j]) {     */
/* 	  dxy=(xstar[j]-xpoint[i])*(xstar[j]-xpoint[i])+(ystar[j]-ypoint[i])*(ystar[j]-ypoint[i]); */
/*           if (dxy < tol*tol) { */
/*             if (dxys < 0.0) { */
/*               dxys = dxy; */
/*               igs = j; */
/*             } */
/*             else if (dxy < dxys) { */
/*               dxys = dxy; */
/*               igs = j; */
/*             } */
/*           } */
/*         } */
/*       } */
/*       if (igs > -1) { */
/*         alfa[nm]=ar[igs]; */
/*         delta[nm]=dec[igs]; */
/*         xpix[nm]=xpoint[i]; */
/*         ypix[nm]=ypoint[i]; */
/*         nm++; */
/*       } */
/*     } */
    /*     Lo hago TRES veces, el matching */
    
/*     ncross=nm; */
    if(nstar!=0) PlotStars(ar,dec,mag,nstar,alfa,delta);
/*     //for(i=0;i<ncross;i++) printf(" %d xpix %f ypix %f\n",i,xpix[i],ypix[i]); */
    cpgsvp(0.1,0.5,0.1,0.9);
    cpgwnad(0.,naxes[0],0.,naxes[1]);

    Plotd(xpix,ypix,ncross,2,2);
/*     Plotd(xstar,ystar,nstar,6,6); */
    printf(" Number of matchings %d\n",ncross);

    if(ncross>ngrad) {
      fitwcs=1;       
      ComputeSolution(alfa,delta,xpix,ypix);
      /*     exit(1); */
      pix2wcs(wcsim,naxes[0]/2., naxes[1]/2. ,&dnul1, &dnul2);
      alfac=(float)dnul1;deltac=(float)dnul2;
      xpixsiz=fabs((float)((*wcsim).xinc)*3600.);
      /*     printf(" HOAL \n"); */
      if(strsrch(seccat,"NONE") == NULL)
	{
	  printf(" Main WCS solution set. Now trying to compute secondary correction \n"); 
	  if(strsrch(seccat,"ACT") != NULL) 
	    FindACT(ars,decs,mags, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstarsec,0.,0.,1.,1.,NSECMAX);
	  else if(strsrch(seccat,"Tycho-2") != NULL)
	    FindTycho2(ars,decs,mags, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstarsec,0.,0.,1.,1.,NSECMAX);
	  printf(" Matching secondary stars with objects\n");
	  if(nstarsec>0) {
	    /* 	  printf(" MEJR\n"); */
	    tol=5./xpixsiz;    /* Permitimos 5 segundos de arco como maximo */
	    /* 	  printf(" xpixsiz %f\n",xpixsiz); */
	    ncrosssec=0.;
	    for(i=0;i<nstarsec;i++) {
	      xstar[i]=0;ystar[i]=0;
	      wcs2pix(wcsim,ars[i],decs[i],&xstars[i],&ystars[i],&offs[i]);
	      dxys=tol*tol*2.;
	      jgs=-1;
	      if (!off[i]) {
		for(j=0;j<npoint;j++) {
		  dxy=(xstars[i]-xpoint[j])*(xstars[i]-xpoint[j])+(ystars[i]-ypoint[j])*(ystars[i]-ypoint[j]);
		  /* 		printf("%d %d dxy %f\n",i,j,dxy); */
		  if(dxy < tol*tol && dxy < dxys) {
		    jgs=j;
		    dxys=dxy;
		  }
		}
	      }
	      if(jgs > -1 ) {
		alfas[ncrosssec]=ars[i];
		deltas[ncrosssec]=decs[i];
		xpixs[ncrosssec]=xpoint[jgs];
		ypixs[ncrosssec]=ypoint[jgs];
		ncrosssec++;
	      }
	    }
	  }
	  if(ncrosssec>0) {
	    
	    secflag=1;
	  }
	  else {
	    printf(" Not possible to use secondary stars\n");
	    secflag=0;
	  }
	  
	}
      else secflag=0;

      printf(" Number of secondary stars used for correction: %d\n",ncrosssec);
      
      printf(" Shifting solution\n");
      
      if (secflag) {
	PlotStars(ar,dec,mag,nstar,alfa,delta);
	ShiftSolution(alfas,deltas,xpixs,ypixs);
      }
      printf(" Changing central AR DEC from %f %f\n",alfac,deltac);
      alfac=wcsim->xref;deltac=wcsim->yref;
      pix2wcs(wcsim,naxes[0]/2., naxes[1]/2. ,&dnul1, &dnul2);
      alfac=(float)dnul1;deltac=(float)dnul2;
      xpixsiz=fabs((float)((*wcsim).xinc)*3600.);
      ypixsiz=fabs((float)((*wcsim).yinc)*3600.);
      printf(" to %f %f\n",alfac,deltac);
    }
    else {
      fitwcs=0;
    }
    if(iswcs(wcsim) && fitwcs) {
      printf(" WCS solution has been computed\n");
      /*     PrintWCS(header,1); */
      /*     printf(" Cabarao\n"); */
/*       SetFITSWCS(header,wcsim); */
      /*     printf(" Pasa de aqui\n"); */
/*       //PrintWCS(header,1); */
      printf(" This solution will be set in the header\n");
      WriteWCS_keys();
      printf(" Header succesfully updated\n");
      /*      printf(" An error ocurred while writing the keywords\n"); */
      sprintf(touchchar,"/bin/touch %s.wcsdone\n",inputfile);
      system(touchchar);
    }
    else printf(" No WCS imformation set. NO WRITING\n");
    if(ffclos(fitsimage,&status)) fits_report_error(stderr,status);
    cpgend();
    exit(0);
  }

  

  /* Menu principal */ 
  while(option!='E') { 
    printf("\n\n Reference stars input (actual number= %d)\n",nstar); 
    printf(" G Find star coordinates in the GSC\n"); 
    printf(" T Find star coordinates in the APM\n"); 
    printf(" U Find star coordinates in the USNO\n"); 
    printf(" J Find star coordinates in the ACT\n"); 
    printf(" O Find star coordinates in the Tycho-2\n"); 
    printf(" V Find secondary reference stars in ACT\n");
    printf(" Z Find secondary reference stars in Tycho-2\n");
    printf(" F Input star coordinates by a file\n"); 
    printf(" K Input star coordinates by keyboard\n"); 
    printf(" Reference points at the image (actual number= %d)\n",npoint); 
    printf(" B Catalogue file with positions of objects\n"); 
    printf(" M Mark points at image\n"); 
    printf(" Y Input points by keyboard\n"); 
    printf(" Cross correlating reference points and starts (actual number = %d, secondary = %d)\n",ncross,ncrosssec); 
    printf(" A Automatic cross correlation\n"); 
    printf(" Q Match stars with current WCS solution given a tolerance\n");
    printf(" C Click, find nearest object in file (B) and input coords. by kbd\n");
    printf(" I Click on image and give coordinates by keyboard\n");
    printf(" P Click on image, find nearest in (B)  and ra-dec plane to find nearest star\n"); 
    printf(" L Click on image and find nearest star in ra-dec plane\n"); 
    printf(" D Click on image and delete reference point\n");
    printf(" N Delete all reference points\n");
    printf(" Others\n"); 
    printf(" H Change cuts\n");
    printf(" R Reset initial plate solution\n");
    printf(" C Compute astrometric solution with current positions\n"); 
    printf(" X Test astrometric solution cliking\n"); 
    printf(" W Write solution in image\n"); 
/*     printf(" S Set coordinates for file in (B)\n");  */
    printf(" S Save parameter file\n"); 
    printf(" E Exit\n");
    option=readc('C'); 
    switch (option) { 
    case 'N':
    case 'n':
      ncross=0;
      break;
    case 'S':
    case 's':
      printf("Name of file  :");
      reads(parfilename,parfilename);
      SaveParam(parfilename);
      break;
    case 'R':
    case 'r':
      initguess=0;
      InputInitGuess(&initguess);
      /* wcsim=wcsxinit((double)(alfac),(double)(deltac),(double)(xpixsiz),naxes[0]/2.,naxes[1]/2.,naxes[0],naxes[1],0.,2000,epoch,"TAN"); */
      /*       No se si es necesario  esto que esta comentado. Creo que no*/
      break;
    case 'H':
    case 'h':
      printf(" Current cuts: %f %f\n",datamin,datamax);
      printf(" New background: ");
      datamin=readf(datamin);
      printf(" New foreground: ");
      datamax=readf(datamax);
      break;
    case 'W':
    case 'w':
      if(iswcs(wcsim)) printf(" WCS information set\n");
      else printf(" No WCS imformation set. NO WRITING\n");
      SetFITSWCS(header,wcsim);
/*       //hputs(header,"HISTORY","WCS solution set by PlateAstrom. Part of OPERA with the help of the WCSTools library"); */
      PrintWCS(header,1);
      printf(" This information will be set in the header\n");
      if(status) fits_report_error(stderr,status);
      /*       return; */
      status=0;
      WriteWCS_keys();

      if(status) printf(" Header succesfully updated\n");
      else printf(" An error ocurred while writing the keywords\n");

      break; 
    case 'P':
    case 'p':
      ClickBoth(ar,dec,nstar,xpoint,ypoint,npoint,alfa,delta,xpix,ypix,naxes[0],naxes[1]);
      printf(" Ha pasado de aqui\n"); 
      break;  
    case 'B' : 
    case 'b' : 
      printf("Input file with X,Y position of objects [%s]: ",bfile);
      reads(bfile,bfile);
/*       scanf("%s",bfile); */
      printf("Input column with x (pixel ) coordinate: ");
      bxcol=readi(bxcol);
      printf("Input column with y (pixel ) coordinate: ");
      bycol=readi(bycol);
      printf("Input column with flux: ");
      bfcol=readi(bfcol);
      printf(" Input fluxes range (none if equal); ");
      flux1=readf(flux1);flux2=readf(flux2);
      Buho(xpoint,ypoint,&npoint,flux1,flux2); 
      break; 
    case '1':
      if (secflag) {
	PlotStars(ar,dec,mag,nstar,alfa,delta);
	ShiftSolution(alfas,deltas,xpixs,ypixs);
      }
      printf(" Changing central AR DEC from %f %f\n",alfac,deltac);
      alfac=wcsim->xref;deltac=wcsim->yref;
      pix2wcs(wcsim,naxes[0]/2., naxes[1]/2. ,&dnul1, &dnul2);
      alfac=(float)dnul1;deltac=(float)dnul2;
      xpixsiz=fabs((float)((*wcsim).xinc)*3600.);
      ypixsiz=fabs((float)((*wcsim).yinc)*3600.);
      printf(" to %f %f\n",alfac,deltac);
      break;
    case 'C': 
    case 'c': 
      fitwcs=1; 
      ComputeSolution(alfa,delta,xpix,ypix);
      if (secflag) {
	PlotStars(ar,dec,mag,nstar,alfa,delta);
	ShiftSolution(alfas,deltas,xpixs,ypixs);
      }
      printf(" Changing central AR DEC from %f %f\n",alfac,deltac);
      alfac=wcsim->xref;deltac=wcsim->yref;
      pix2wcs(wcsim,naxes[0]/2., naxes[1]/2. ,&dnul1, &dnul2);
      alfac=(float)dnul1;deltac=(float)dnul2;
      xpixsiz=fabs((float)((*wcsim).xinc)*3600.);
      ypixsiz=fabs((float)((*wcsim).yinc)*3600.);
      printf(" to %f %f\n",alfac,deltac);
      break;
    case 'I':
    case 'i':
      ClickImage(xpoint,ypoint,alfa,delta,&nstar);
      break;
    case 'L':
    case 'l':
      ClickImFindStar(ar,dec,nstar,alfa,delta,xpix,ypix,naxes[0],naxes[1]);
      break;
    case 'F': 
    case 'f': 
      File(xpoint,ypoint,alfa,delta,&nstar);
      break;
    case 'K': 
    case 'k': 
      Keyboard(xpoint,ypoint,alfa,delta, &nstar);
      break;
    case 'U':
    case 'u':
      if(initguess==0)  {
	InputInitGuess(&initguess);
/* 	wcsim=wcsxinit((double)(alfac),(double)(deltac),(double)(xpixsiz),naxes[0]/2.,naxes[1]/2.,naxes[0],naxes[1],0.,2000,epoch,"TAN"); */
      }
      printf(" Input magnitude range (none if equal):");
      mag1=readd(mag1);mag2=readd(mag2);
      printf(" Rango;: %g %8.4g %g\n",mag1,mag1,mag2);
      printf(" Input factor in RA for searching box: ");
      factra=readf(factra);
      printf(" Input factor in Dec for searching box: ");
      factdec=readf(factdec);
/*       printf(" Para aui\n"); */
      FindUSNO(ar,dec,mag, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstar,mag1,mag2,factra,factdec,NMAX);
/*       printf(" AASS\n"); */
      strcpy(astcat,"USNO");
      break;
    case 'J':
    case 'j':
      if(initguess==0)  {
	InputInitGuess(&initguess);
/* 	wcsim=wcsxinit((double)(alfac),(double)(deltac),(double)(xpixsiz),naxes[0]/2.,naxes[1]/2.,naxes[0],naxes[1],0.,2000,epoch,"TAN"); */
      }
      printf(" Input magnitude range (none if equal):");
      mag1=readd(mag1);mag2=readd(mag2);
      printf(" Input factor in RA for searching box: ");
      factra=readf(factra);
      printf(" Input factor in Dec for searching box: ");
      factdec=readf(factdec);
      FindACT(ar,dec,mag, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstar,mag1,mag2,factra,factdec,NMAX);
      strcpy(astcat,"ACT");
      break;
    case 'O':
    case 'o':
      if(initguess==0)  {
	InputInitGuess(&initguess);
/* 	wcsim=wcsxinit((double)(alfac),(double)(deltac),(double)(xpixsiz),naxes[0]/2.,naxes[1]/2.,naxes[0],naxes[1],0.,2000,epoch,"TAN"); */
      }
      printf(" Input magnitude range (none if equal):");
      mag1=readd(mag1);mag2=readd(mag2);
      printf(" Input factor in RA for searching box: ");
      factra=readf(factra);
      printf(" Input factor in Dec for searching box: ");
      factdec=readf(factdec);
      FindTycho2(ar,dec,mag, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstar,mag1,mag2,factra,factdec,NMAX);
      strcpy(astcat,"Tycho-2");
      break;
    case 'G':
    case 'g':
      if(initguess==0)  {
	InputInitGuess(&initguess);
	wcsim=wcsxinit((double)(alfac),(double)(deltac),(double)(xpixsiz),naxes[0]/2.,naxes[1]/2.,naxes[0],naxes[1],0.,2000,epoch,"TAN");
      }
      printf(" Input magnitude range (none if equal):");
      mag1=readd(mag1);mag2=readd(mag2);
      printf(" Input factor in RA for searching box: ");
      factra=readf(factra);
      printf(" Input factor in Dec for searching box: ");
      factdec=readf(factdec);
      FindGSC(ar,dec,mag, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstar,mag1,mag2,factra,factdec,NMAX);
      strcpy(astcat,"GSC");
      break;
    case 'T':
    case 't':
      if(initguess==0) {
        InputInitGuess(&initguess);
/* 	wcsim=wcsxinit((double)(alfac),(double)(deltac),(double)(xpixsiz),naxes[0]/2.,naxes[1]/2.,naxes[0],naxes[1],0.,2000,epoch,"TAN"); */
      }
      printf(" XPIX %f YPIX %f\n",xpixsiz,ypixsiz);
      printf(" Input magnitude range (none if equal):");
      mag1=readd(mag1);mag2=readd(mag2);
      printf(" Input factor in RA for searching box: ");
      factra=readf(factra);
      printf(" Input factor in Dec for searching box: ");
      factdec=readf(factdec);
      printf(" xpixsiz %f ypixsiz %f\n",xpixsiz,ypixsiz);
      FindAPM(ar,dec,mag, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstar,mag1,mag2,factra,factdec,NMAX);
      strcpy(astcat,"APM");
      break;
    case 'V':
    case 'v':
      if(initguess==0 || nstar==0 || fitwcs==0)  {
	printf(" There are no primary reference stars or solution not computed, no use in looking for secondary stars.\n");
	break;
      }
      printf(" Input magnitude range (none if equal):");
      mag1=readd(mag1);mag2=readd(mag2);
      printf(" Important: setting epoch to %4.2f\n",epoch);
      FindACT(ars,decs,mags, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstarsec,mag1,mag2,1.,1.,NSECMAX);
      secflag=1;

      tol=5./xpixsiz;    /* Permitimos 5 segundos de arco como maximo */
      ncrosssec=0.;
      for(i=0;i<nstarsec;i++) {
	xstar[i]=0;ystar[i]=0;
	wcs2pix(wcsim,ars[i],decs[i],&xstars[i],&ystars[i],&offs[i]);
	dxys=tol*tol*2.;
	jgs=-1;
	if (!off[i]) {
	  for(j=0;j<npoint;j++) {
	    dxy=(xstars[i]-xpoint[j])*(xstars[i]-xpoint[j])+(ystars[i]-ypoint[j])*(ystars[i]-ypoint[j]);
	    if(dxy < tol*tol && dxy < dxys) {
	      jgs=j;
	      dxys=dxy;
	    }
	  }
	}
	if(jgs > -1 ) {
	  alfas[ncrosssec]=ars[i];
	  deltas[ncrosssec]=decs[i];
	  xpixs[ncrosssec]=xpoint[jgs];
	  ypixs[ncrosssec]=ypoint[jgs];
	  ncrosssec++;
	}
      }
      if(ncrosssec>0) {
	strcpy(seccat,"ACT");
	secflag=1;
      }
      else {
	printf(" Not possible to use secondary stars\n");
	secflag=0;
      }

      break;
    case 'Z':
    case 'z':
      if(initguess==0 || nstar==0 || fitwcs==0)  {
	printf(" There are no primary reference stars or solution not computed, no use in looking for secondary stars.\n");
	break;
      }
      printf(" Input magnitude range (none if equal):");
      mag1=readd(mag1);mag2=readd(mag2);
/*       printf(" Input factor for box: "); */
/*       fact=readf(fact); */
      printf(" Important: setting epoch to %4.2f\n",epoch);
      FindTycho2(ars,decs,mags, alfac,deltac,xpixsiz*naxes[0]/cos(deltac/180*3.1415), ypixsiz*naxes[1],&nstarsec,mag1,mag2,1.,1.,NSECMAX);
      secflag=1;

/*       strcpy(astcat,"ACT"); */
      tol=5./xpixsiz;    /* Permitimos 5 segundos de arco como maximo */
      ncrosssec=0.;
      for(i=0;i<nstarsec;i++) {
	xstar[i]=0;ystar[i]=0;
	wcs2pix(wcsim,ars[i],decs[i],&xstars[i],&ystars[i],&offs[i]);
	dxys=tol*tol*2.;
	jgs=-1;
	if (!off[i]) {
	  for(j=0;j<npoint;j++) {
	    dxy=(xstars[i]-xpoint[j])*(xstars[i]-xpoint[j])+(ystars[i]-ypoint[j])*(ystars[i]-ypoint[j]);
/* 	    printf("%d %d dxy %f\n",i,j,dxy); */
	    if(dxy < tol*tol && dxy < dxys) {
	      jgs=j;
	      dxys=dxy;
	    }
	  }
	}
	if(jgs > -1 ) {
	  alfas[ncrosssec]=ars[i];
	  deltas[ncrosssec]=decs[i];
	  xpixs[ncrosssec]=xpoint[jgs];
	  ypixs[ncrosssec]=ypoint[jgs];
	  ncrosssec++;
	}
      }
      if(ncrosssec>0) {
	strcpy(seccat,"Tycho-2");
	secflag=1;
      }
      else {
	printf(" Not possible to use secondary stars\n");
	secflag=0;
      }

      break;
    case 'E':
    case 'e':
      if(ffclos(fitsimage,&status)) fits_report_error(stderr,status);
      sprintf(touchchar,"/bin/touch %s.wcsdone\n",inputfile);
      system(touchchar);
      cpgclos();
      exit(1);
      break;
    case 'A':
    case 'a':
      tol=20.;
      while(tol!=0) {
//	for(i=0;i<npoint;i++) printf("ANTEs %d xpo %f ypo %f\n",i,xpoint[i],ypoint[i]);
//	for(i=0;i<nstar;i++) printf("ANTES %d ra %f dec %f mag %f xst %f yst %f\n",i,ar[i],dec[i],mag[i],xstar[i],ystar[i]);

	StarMatch_ceg(&ncross, alfa,delta,xpix,ypix,off,npoint,xpoint,ypoint,nstar,ar,dec,xstar,ystar,tol,wcsim,0);
	printf(" Sale con ncross %d\n",ncross);
	printf(" WCS after the matching\n");
	SetFITSWCS(header,wcsim);
	PrintWCS(header,1);

	for(i=0;i<npoint;i++) printf(" %d xpo %f ypo %f\n",i,xpoint[i],ypoint[i]);
	for(i=0;i<nstar;i++) printf(" %d ra %f dec %f xst %f yst %f\n",i,ar[i],dec[i],xstar[i],ystar[i]);

	printf(" Todo con ceg\n");

	if(nstar!=0) PlotStars(ar,dec,mag,nstar,alfa,delta);
	Plotd(xstar,ystar,nstar,6,6);
	printf(" Number of matchings %d\n",ncross);
	printf(" New tolerance (0=exit): ");
	tol=readf(0);
      }
      fitwcs=1;
      break;
    case 'Q':
    case 'q':
      if(fitwcs==1 || initguess==1) {
	tol=20.;
	while(tol!=0) {
	  MatchTol(&ncross, alfa,delta,xpix,ypix,off,npoint,xpoint,ypoint,nstar,ar,dec,xstar,ystar,tol,wcsim,0);
	  if(nstar!=0) PlotStars(ar,dec,mag,nstar,alfa,delta);
	  Plotd(xstar,ystar,nstar,6,6);
	  printf(" Number of matchings %d\n",ncross);
	  printf(" New tolerance (0=exit): ");
	  tol=readf(0);
	}
      }
      else printf(" Not solution set yet\n");
      break;
      
    }
    nbrighteststar=nstar;
    nbrightestobj=npoint;
    
    cpgeras();
/*     //Dibujo la imagen */
    cpgsvp(0.1,0.5,0.1,0.9);
    cpgwnad(0.,naxes[0],0.,naxes[1]);
    cpggray(buffer,naxes[0],naxes[1],1,naxes[0],1,naxes[1],datamax,datamin,tr);
    cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
    strcpy(snul,"Image \0");
    strcat(snul,inputfile);
    cpglab("X axis","Y axis",snul);
/*     //printf(" DDibuja\n"); */

/*     //Cosas que hay que hacer cada vez que pasamos por el menu */

/*     //npoint=10; */
    
/*     //printf(" paso 1 \n"); */
    Plotd(xpoint,ypoint,npoint,4,4);
/*     //printf(" paso 2 \n"); */
/*     //    Plotd(xstar,ystar,nstar,3,3); */
/*     //printf(" paso 3 \n"); */
    Plotd(xpix,ypix,ncross,2,2);
/*     //printf(" paso 4 \n"); */
    if(nstar!=0) PlotStars(ar,dec,mag,nstar,alfa,delta);
/*     //if(nstar!=0) PlotStars(ar,dec,mag,nstar,alfacal,deltacal); */
/*     //printf(" paso 5 \n"); */
    

    Plotd(xstar,ystar,nstar,6,6);

/*     //    PlotStartsRef(); */
/*     //PlotPointRef(); */
/*     //PlotCroos(); */

  }
  
     


   /* Una vez que tengo ya todas las estrellas */
/*   ctes.ap=23; */
/*   printf("ANTESSSSSS\n"); */

/*   //  PlateCtes_3(nstar, xpoint, ypoint, alfa, delta,alfac, deltac, &ctes); */
/*   printf("DESPIESSSS\n"); */
/*   printf(" la A %f\n",ctes.ap); */
/*   printf(" Plate solution:\n    psi= %e + %e x+ %e y \n    eta= %e + %e x+ %e y \n",ctes.ap,ctes.bp,ctes.cp,ctes.ae,ctes.be,ctes.ce); */
/*   fflush(NULL); */
/*   fflush(stdin); */
/*   factorarr=300; */
/*   while(factorarr!=0) { */
/*     for(i=0;i<nstar;i++) { */
/*       //Ecu2Plac(alfa[i],delta[i],alfac,deltac,&psi,&eta); */
/*       xcal=(ctes.ce*psi-ctes.cp*eta-ctes.ce*ctes.ap+ctes.cp*ctes.ae)/(-ctes.be*ctes.cp+ctes.bp*ctes.ce); */
      
/*       ycal=(ctes.be*psi-ctes.bp*eta-ctes.be*ctes.ap+ctes.bp*ctes.ae)/(ctes.be*ctes.cp-ctes.bp*ctes.ce); */
/*       printf("xcal %f ycal %f\n",xcal,ycal); */
/*       printf("psi    %f eta    %f\n",psi,eta); */
/*       cpgarro(xpoint[i],ypoint[i],xpoint[i]-(xpoint[i]-xcal)*factorarr,ypoint[i]-(ypoint[i]-ycal)*factorarr); */
      
/*       psi=ctes.bp*xpoint[i]+ctes.cp*ypoint[i]+ctes.ap; */
/*       eta=ctes.be*xpoint[i]+ctes.ce*ypoint[i]+ctes.ae; */
/*       printf("psical %f etacal %f\n",psi,eta); */
/*       //Plac2Ecu(psi,eta,alfac,deltac,alfacal+i,deltacal+i); */
/*       //error[i]=sqrt(((alfa[i]-alfacal[i])*(alfa[i]-alfacal[i])+(delta[i]-deltacal[i])*(delta[i]-deltacal[i])))/pi*180*3600; */
/*       //printf("error in arcseconds %f\n",error[i]); */
      
/*     } */
/*     fflush(NULL); */
/*     fflush(stdin); */
/*     printf("Enter another factor for arrows [%f] :",factorarr); */
/*     scanf("%f",&factorarr); */
/*     factorarr=0; */
/*   } */
/*   //mean=StMedia(nstar,error,&sigma); */
/*   //printf(" Error is distributed with mean %e arcseconds and deviation %f arcseconds\n",mean,sigma); */
  



  cpgend();
  return(0);
   
}



void Keyboard(double *xpoint, double *ypoint, double *alfa, double *delta,int *nstar)
{
  float pi=4*atan(1);

  int i;
  double ymax,xmax;
  double ymin,xmin;
  float xcen,ycen;
  float a,b;
  /* Pregunto por la estrella a utilizar para la astrometria */
  
  i=0;
/*   //  printf("Input coordinates for the center of the plate\n"); */
/*   //  alfac=Readhms(); */
/*   //  deltac=Readgms(); */
 do {
  printf("Input coordinates x,y in pixels (0 0 = exit): ");
  xpoint[i]=(double)readf(0);  ypoint[i]=(double)readf(0);
/*   scanf(" %f %f",xpoint+i,ypoint+i); */
/*   //printf("x, y %f %f\n",xpoint[i],ypoint[i]); */
    
   alfa[i]=Readhms();
   delta[i]=Readgms();
/*    //printf("es   %d ala fre %e %e\n",i,alfa[i],delta[i]); */
    
   i++;
    
  } while(!(xpoint[i] == 0 && ypoint[i]== 0));
  
  
  *nstar=i;
  
  MinMax_d(*nstar,xpoint,&xmin,&xmax);
  MinMax_d(*nstar,ypoint,&ymin,&ymax);
  xcen=xmax-xmin;
  ycen=ymax-ymin;
  b=(xpoint[1]*alfa[0]-xpoint[0]*alfa[1])/(xpoint[1]-xpoint[0]);
  a=(alfa[0]-alfa[1])/(xpoint[0]-xpoint[1]);
  alfac=a*xcen+b;
  b=(xpoint[1]*delta[0]-xpoint[0]*delta[1])/(xpoint[1]-xpoint[0]);
  a=(delta[0]-delta[1])/(xpoint[0]-xpoint[1]);
  deltac=a*ycen+b;
  printf("alfac %f deltac %f\n",alfac*180/pi/15,deltac*180/pi);



}
void Buho(double *xpoint, double *ypoint, int *npoint,float flux1,float flux2)
{

/*   float pi=4*atan(1); */
/*   //  static int intcompare; */

/*   //char bfile[51]; */
  int nlin,i;
  float *xp,*yp;
  float *flux;
  int *log;
  struct obj_table *todo;
/*   //static int xcol=14,ycol=15,fcol=4; */
/*   //static float flux1=0,flux2=0; */
  float fluxmin=1e35,fluxmax=-1e35;

  nlin=FileNLin(bfile);

  if(nlin> NOBJDET) {
    printf(" More than %d objects. Recompile to allow more objects\n",NOBJDET);
    exit(1);
  }

  xp=vector_f(nlin);
  yp=vector_f(nlin);
  flux=vector_f(nlin);
  log=vector_i(nlin);
  if((todo=malloc(nlin*sizeof(struct obj_table)))==NULL) {
    printf("I cannot dimension todo of %d bytes \n",nlin*sizeof(struct obj_table));
    exit(1);
  }

  for(i=0;i<nlin;i++)     log[i]=0;
  ReadNumcol(bfile,bxcol,xp,log,&nlin);
  ReadNumcol(bfile,bycol,yp,log,&nlin);
  ReadNumcol(bfile,bfcol,flux,log,&nlin);

  for(i=0;i<nlin;i++) {
    todo[i].xpos=xp[i];
    todo[i].ypos=yp[i];
    todo[i].flux=flux[i];
    todo[i].log=log[i];
  }
  qsort(todo,nlin,sizeof(struct obj_table),incompare);
  for(i=0;i<nlin;i++) {
    xp[i]=  todo[i].xpos;
    yp[i]=  todo[i].ypos;
    flux[i]=todo[i].flux;
    log[i]= todo[i].log;
  }
  *npoint=-1;
  for(i=0;i<nlin;i++)  {
    if(log[i])
      {
	if((flux[i]>flux1 && flux[i]<flux2) || flux1==flux2) {
	  if(flux[i]<fluxmin) fluxmin=flux[i];
	  if(flux[i]>fluxmax) fluxmax=flux[i];
	  (*npoint)++;
	  xpoint[*npoint]=xp[i];
	  ypoint[*npoint]=yp[i];
/* 	  //printf(" npo %d x %f y %f\n",*npoint,xpoint[*npoint],ypoint[*npoint]); */
/* 	  //	printf(" npo %d\n",*npoint); */
	}
      }
  }
  
  (*npoint)++;
  printf("\n Using %d Objects from file %s\n",(*npoint),bfile);
  printf(" Flux range: %f %f \n",fluxmin,fluxmax);
  for(i=0;i<(*npoint);i++)
    {
/*       //printf(" %d %f %f flux %f\n",i,xpoint[i],ypoint[i],flux[i]); */
    }
  free(xp);
  free(yp);
  free(flux);
  free(log);
  free(todo);
}

void ClickImage(double *xpoint, double *ypoint, double *alfa, double *delta,int *nstar)
{
}
void File(double *xpoint, double *ypoint, double *alfa, double *delta,int *nstar)
{
/*   float pi=4*atan(1); */

  char file[51];
  int n,i,ii;
/*   int lxpoint[NMAX],lypoint[NMAX]; */
  int ladat[NMAX],lddat[NMAX];
/*   float arl[NMAX],decl[NMAX]; */
/*   float xpl[NMAX],ypl[NMAX]; */
/*   int xcol, ycol; */
  static int  acol=1, dcol=2;
  printf("Input file with coordinates: ");
  reads(file,file);
/*   scanf("%s",file); */

/* //  printf("Input column with x (pixel ) coordinate: "); */
/* //  scanf("%d",&xcol); */
/* //  printf("Input column with y (pixel ) coordinate: "); */
/* //  scanf("%d",&ycol); */
  printf("Input column with RA (format hh:mm:ss): ");
  acol=readi(acol);
/*   scanf("%d",&acol); */
  printf("Input column with declination (format dd.mmss): ");
  dcol=readi(dcol);
/*   scanf("%d",&dcol); */

/* //  ReadNumcol(file,xcol,xpl,lxpoint,&n); */
/* //  ReadNumcol(file,ycol,ypl,lypoint,&n); */
/* //  ReadNumcol(file,acol,arl,ladat,&n); */
/* //  ReadNumcol(file,dcol,decl,lddat,&n); */
  ReadWCScol(file,acol,alfa,ladat,&n);
  ReadWCScol(file,dcol,delta,lddat,&n);

  ii=-1;
  for(i=0;i<n;i++)
    {
/* //      if(lxpoint[i] && lypoint[i] && ladat[i] && lddat[i]) */
      if( ladat[i] && lddat[i])
	{
	  i++;
/* //	  xpoint[ii]=xpl[i]; */
/* //	  ypoint[ii]=ypl[i]; */
/* 	  alfa[ii]=alfa[i]*15; */
/* //	  delta[ii]=decl[i]; */
	  
	}
    }
  *nstar=ii;
  printf("\n Using %d stars from file %s\n",*nstar,file);
  for(i=0;i<*nstar;i++)
    {
      printf(" %f %f %f %f\n",xpoint[i],ypoint[i],alfa[i],delta[i]);
    }







}







void FindGSC(double *ar,double *dec,double *mag,float alfacen,float deltacen,float xdim,float ydim,  int *nstar,double mag1,double mag2,float factra,float factdec,int nmax)
{

  char fileref[51];
  FILE *fref;
  char snul1[32],snul2[32];

  struct star_table todo[NMAX];
  double gnum[NMAX],gmag[NMAX];
  double arl[NMAX],decl[NMAX];
  int gtype[NMAX];
  int i;  
  int verbose=interact;
/*   static double mag1=0,mag2=0; */
/*   static float fact=1.0;     //Para buscar en una caja un poco mas grande */
  *nstar=0;

  printf(" Searching GSC stars round %f hours %f degrees\n",alfacen/15,deltacen);
  printf("  Box %f x %f arcmin\n",xdim/60.*factra,ydim/60.*factdec);
#ifdef WCSVERSION26
  *nstar=gscread(alfacen,deltacen,xdim/3600./2.*factra,ydim/3600./2.*factdec,0.,WCS_J2000,2000.,epoch,mag1,mag2,-1,nmax,gnum,arl,decl,gmag,gtype,verbose);
#endif
#ifdef WCSVERSION28
  *nstar=gscread(alfacen,deltacen,xdim/3600./2.*factra,ydim/3600./2.*factdec,0.,0,WCS_J2000,2000.,epoch,mag1,mag2,nmax,gnum,arl,decl,gmag,gtype,verbose);
#endif

  printf(" Antes \n");
  printf(" Found stars: %d\n",*nstar);

  for(i=0;i<*nstar;i++) {
    todo[i].num= gnum[i];
    todo[i].ar =  arl[i];
    todo[i].dec= decl[i];
    todo[i].mag=gmag[i];
  }
  qsort(todo,*nstar,sizeof(struct star_table),iscompare);
  for(i=0;i<*nstar;i++) {
    arl[i] =  todo[i].ar;
    decl[i]=  todo[i].dec; 
    gmag[i]=  todo[i].mag;
    gnum[i]=  todo[i].num;
/*     //printf(" Star %8.0f AR %g Dec %g  Mag %f\n",todo[i].num,todo[i].ar,todo[i].dec,todo[i].mag); */

  }

  strcpy(fileref,inputfile);
  strcat(fileref,".gsc");

  if((fref=fopen(fileref,"w"))==NULL) {
    printf("\nERROR: Can't open file %s\n",fileref);
    exit(1);
  }
  printf(" Saving reference stars used for astrometry in file %s\n",fileref);
  fprintf(fref,"#Reference stars in the searched file from the GSC catalogue\n");
  fprintf(fref,"#round %f hours %f degrees, box %f x %f arcmin for epoch %f\n",alfacen/15,deltacen,xdim/60.*factra,ydim/60.*factdec,epoch);
  fprintf(fref,"#Num     RA(J2000)   DEC(J2000)      Mag  \n");
  for(i=0;i<*nstar;i++) {
    ar[i]=arl[i];dec[i]=decl[i];mag[i]=gmag[i];
/*     printf(" Star %8.0f AR %g Dec %g  Mag %f\n",gnum[i],ar[i],dec[i],gmag[i]); */
    ra2str(snul1,32,ar[i],2);
    dec2str(snul2,32,dec[i],1);
    fprintf(fref," %8.0f %s %s %6.2f\n",gnum[i],snul1,snul2,gmag[i]);
  }
  fclose(fref);
}



void FindUSNO(double *ar,double *dec,double *mag,float alfacen,float deltacen,float xdim,float ydim,  int *nstar,double mag1,double mag2,float factra,float factdec,int nmax)
{

  char fileref[51];
  FILE *fref;
  char snul1[32],snul2[32];

  struct star_table todo[NMAX];
  double gnum[NMAX],gmag[NMAX],gmagb[NMAX];
  double arl[NMAX],decl[NMAX];
  int platenum[NMAX];
  int i;  
  int verbose=interact;
  static int savedflag=0;
  *nstar=0;

  printf(" Searching USNO stars round %f hours %f degrees\n",alfacen/15,deltacen);
  printf("  Box %f x %f arcmin\n",xdim/60.*factra,ydim/60.*factdec);
#ifdef WCSVERSION26
  *nstar=uacread("UAC2",alfacen,deltacen,xdim/3600./2.*factra,ydim/3600./2.*factdec,0.,WCS_J2000,2000.,epoch,mag1,mag2,0,nmax,gnum,arl,decl,gmag,gmagb,platenum,verbose);
#endif
#ifdef WCSVERSION28
  *nstar=uacread("UAC2",0,alfacen,deltacen,xdim/3600./2.*factra,ydim/3600./2.*factdec,0.,WCS_J2000,2000.,epoch,mag1,mag2,nmax,gnum,arl,decl,gmag,gmagb,platenum,verbose);
#endif


  /*      Esta puesto igual para la magnitud R y la B. */ 

  printf(" Found stars: %d\n",*nstar);
  for(i=0;i<*nstar;i++) {
    todo[i].num= gnum[i];
    todo[i].ar =  arl[i];
    todo[i].dec= decl[i];
    todo[i].mag=gmag[i];
    todo[i].mag2=gmagb[i];


     /*     printf(" Star %14.8g AR %g Dec %g  Mag %f Margb %f\n",gnum[i],arl[i],decl[i],gmag[i],gmagb[i]); */ 
    
  }
  qsort(todo,*nstar,sizeof(struct star_table),iscompare);

  for(i=0;i<*nstar;i++) {
    arl[i] =  todo[i].ar;
    decl[i]=  todo[i].dec; 
    gmag[i]=  todo[i].mag;
    gmagb[i]=  todo[i].mag2;
    gnum[i]=  todo[i].num;


     /*     printf(" Star %8.0f AR %g Dec %g  Mag %f\n",todo[i].num,todo[i].ar,todo[i].dec,todo[i].mag); */ 

  }
  if(!savedflag) {

    strcpy(fileref,inputfile);
    strcat(fileref,".usno");

    if((fref=fopen(fileref,"w"))==NULL) {
      printf("\nERROR: Can't open file %s\n",fileref);
      exit(1);
    }
    printf(" Saving reference stars used for astrometry in file %s\n",fileref);
    fprintf(fref,"#Reference stars in the searched file from the USNO catalogue\n");
    fprintf(fref,"#round %f hours %f degrees, box %f x %f arcmin for epoch %f\n",alfacen/15,deltacen,xdim/60.*factra,ydim/60.*factdec,epoch);
    fprintf(fref,"#Num           RA(J2000)   DEC(J2000)      MagR     MagB    Plate N.\n");
    for(i=0;i<*nstar;i++) {
      ar[i]=arl[i];dec[i]=decl[i];mag[i]=gmag[i];
      ra2str(snul1,32,ar[i],2);
      dec2str(snul2,32,dec[i],1);
      fprintf(fref," %14.8f %s %s %6.2f %6.2f  %8d\n",gnum[i],snul1,snul2,gmag[i],gmagb[i],platenum[i]);
    }
    fclose(fref);
    savedflag=1;
  }
}

void FindACT(double *ar,double *dec,double *mag,float alfacen,float deltacen,float xdim,float ydim,  int *nstar,double mag1,double mag2,float factra,float factdec,int nmax)
{

  char fileref[51];
  FILE *fref;
  char snul1[32],snul2[32];

  struct star_table todo[NMAX];
  double gnum[NMAX],gmag[NMAX],gmagb[NMAX];
  double gpra[NMAX],gpdec[NMAX];
  double arl[NMAX],decl[NMAX];
  int gtype[NMAX];
  int i;  
  int verbose=interact;
/*   //  static double mag1=0,mag2=0; */
/*   //  static float fact=1.0;     //Para buscar en una caja un poco mas grande */
  *nstar=0;

  printf(" Searching ACT stars round %f hours %f degrees\n",alfacen/15,deltacen);
  printf("  Box %f x %f arcmin\n",xdim/60.*factra,ydim/60.*factdec);
/*   //  printf(" Input magnitude range (none if equal):"); */
/*   //  mag1=(double)readf((float)mag1);mag2=(double)readf((float)mag2); */
/*   //  printf(" Input factor for box: "); */
/*   //  fact=readf(fact); */

#ifdef WCSVERSION26
  *nstar=actread(alfacen,deltacen,xdim/3600./2.*factra,ydim/3600./2.*factdec,0.,WCS_J2000,2000.,epoch,mag1,mag2,nmax,gnum,arl,decl,gmag,gmagb,gtype,verbose);
#endif

#ifdef WCSVERSION28
  *nstar=actread(alfacen,deltacen,xdim/3600./2.*factra,ydim/3600./2.*factdec,0.,0,WCS_J2000,2000.,epoch,mag1,mag2,nmax,gnum,arl,decl,gpra,gpdec,gmag,gmagb,gtype,verbose);
#endif


/*   //Esta puesto igual para la magnitud V y la B.  */
  printf(" Found stars: %d\n",*nstar);
  for(i=0;i<*nstar;i++) {
/*     //printf(" Sin %8.0f AR %g Dec %g  Mag %f\n",gnum[i],arl[i],decl[i],gmag[i]); */
    todo[i].num =gnum[i];
    todo[i].ar=arl[i];
    todo[i].dec=decl[i];
    if(gmag[i]!=0)    todo[i].mag=gmag[i];
    else              todo[i].mag=40.;
  }
  qsort(todo,*nstar,sizeof(struct star_table),iscompare);

  for(i=0;i<*nstar;i++) {

    arl[i] =  todo[i].ar;
    decl[i]=  todo[i].dec; 
    gmag[i]=todo[i].mag;
    gnum[i]= todo[i].num;
/*     //printf(" Con %8.0f AR %g Dec %g  Mag %f\n",gnum[i],arl[i],decl[i],gmag[i]); */
  }
/*   //exit(1); */

  strcpy(fileref,inputfile);
  strcat(fileref,".act");

  if((fref=fopen(fileref,"w"))==NULL) {
    printf("\nERROR: Can't open file %s\n",fileref);
    exit(1);
  }
  printf(" Saving reference stars used for astrometry in file %s\n",fileref);
  fprintf(fref,"#Reference stars in the searched file from the ACT catalogue\n");
  fprintf(fref,"#round %f hours %f degrees, box %f x %f arcmin for epoch %f\n",alfacen/15,deltacen,xdim/60.*factra,ydim/60.*factdec,epoch);
  fprintf(fref,"#Num     RA(J2000)   DEC(J2000)      Mag  \n");
  for(i=0;i<*nstar;i++) {
    ar[i]=arl[i];dec[i]=decl[i];mag[i]=gmag[i];
/*     printf(" Star %8.0f AR %g Dec %g  Mag %f\n",gnum[i],ar[i],dec[i],gmag[i]); */
    ra2str(snul1,32,ar[i],2);
    dec2str(snul2,32,dec[i],1);
    fprintf(fref," %8.0f %s %s %6.2f\n",gnum[i],snul1,snul2,gmag[i]);
  }
  fclose(fref);
}

void FindTycho2(double *ar,double *dec,double *mag,float alfacen,float deltacen,float xdim,float ydim,  int *nstar,double mag1,double mag2,float factra,float factdec,int nmax)
{

  char fileref[51];
  FILE *fref;
  char snul1[32],snul2[32];

  struct star_table todo[NMAX];
  double gnum[NMAX],gmag[NMAX],gmagb[NMAX];
  double gpra[NMAX],gpdec[NMAX];
  double arl[NMAX],decl[NMAX];
  int gtype[NMAX];
  int i;  
  int verbose=interact;
  *nstar=0;

  printf(" Searching Tycho-2 stars round %f hours %f degrees\n",alfacen/15,deltacen);
  printf("  Box %f x %f arcmin\n",xdim/60.*factra,ydim/60.*factdec);
  printf(" na %d %d\n",nmax,NMAX);
#ifdef WCSVERSION28
  *nstar=ty2read(alfacen,deltacen,xdim/3600./2.*factra,ydim/3600./2.*factdec,0.,0,WCS_J2000,2000.,epoch,mag1,mag2,nmax,gnum,arl,decl,gpra,gpdec,gmag,gmagb,gtype,verbose);
#endif
#ifdef WCSVERSION26
  printf(" Tycho-2 catalogue is not available with this version of WCSTOOLS (2.6)\n");
#endif

/*   //Esta puesto igual para la magnitud V y la B.  */
  printf(" Found stars: %d\n",*nstar);
  for(i=0;i<*nstar;i++) {
    todo[i].num =gnum[i];
    todo[i].ar=arl[i];
    todo[i].dec=decl[i];
    if(gmag[i]!=0)    todo[i].mag=gmag[i];
    else              todo[i].mag=40.;
  }
  qsort(todo,*nstar,sizeof(struct star_table),iscompare);

  for(i=0;i<*nstar;i++) {

    arl[i] =  todo[i].ar;
    decl[i]=  todo[i].dec; 
    gmag[i]=todo[i].mag;
    gnum[i]= todo[i].num;
/*     //printf(" Con %8.0f AR %g Dec %g  Mag %f\n",gnum[i],arl[i],decl[i],gmag[i]); */
  }
/*   //exit(1); */

  strcpy(fileref,inputfile);
  strcat(fileref,".ty2");

  if((fref=fopen(fileref,"w"))==NULL) {
    printf("\nERROR: Can't open file %s\n",fileref);
    exit(1);
  }
  printf(" Saving reference stars used for astrometry in file %s\n",fileref);
  fprintf(fref,"#Reference stars in the searched file from the Tycho-2 catalogue\n");
  fprintf(fref,"#round %f hours %f degrees, box %f x %f arcmin for epoch %f\n",alfacen/15,deltacen,xdim/60.*factra,ydim/60.*factdec,epoch);
  fprintf(fref,"#Num     RA(J2000)   DEC(J2000)      Mag  \n");
  for(i=0;i<*nstar;i++) {
    ar[i]=arl[i];dec[i]=decl[i];mag[i]=gmag[i];
/*     printf(" Star %8.0f AR %g Dec %g  Mag %f\n",gnum[i],ar[i],dec[i],gmag[i]); */
    ra2str(snul1,32,ar[i],2);
    dec2str(snul2,32,dec[i],1);
    fprintf(fref," %8.0f %s %s %6.2f\n",gnum[i],snul1,snul2,gmag[i]);
  }
  fclose(fref);
}


void FindAPM(double *ar,double *dec,double *mag,float alfacen,float deltacen,float xdim,float ydim,  int *nstar,double mag1,double mag2,float factra,float factdec,int nmax)
{
/*   char fileapm[51]; */

  char fileref[51];
  FILE *fref;
  char snul1[32],snul2[32];
  struct star_table todo[NMAX];

  double gnum[NMAX],gmag[NMAX];
  double arl[NMAX],decl[NMAX];
  int gtype[NMAX];
  int i;  
  int verbose=interact;
/*   //  static double mag1=0,mag2=0; */
/*   //  static float fact=1.0;     //Para buscar en una caja un poco mas grande */
  *nstar=0;
  printf(" Searching APM stars round %f  hours %f degrees\n",alfacen/15,deltacen);
  printf(" Searching APM stars round %d %d %4.2f   %d %d %4.2f degrees\n"
    ,(int)(alfacen/15),(int)((alfacen/15-(int)(alfacen/15))*60),(((alfacen/15-(int)(alfacen/15))*60)-(int)((alfacen/15-(int)(alfacen/15))*60))*60,
    (int)(deltacen),(int)((deltacen-(int)(deltacen))*60),(((deltacen-(int)(deltacen))*60)-(int)((deltacen-(int)(deltacen))*60))*60);
  printf(" WARNING: not valid for negative\n");

  printf("  Box %f x %f arcmin\n",xdim/60.*factra,ydim/60.*factdec);
/*   //  printf(" Input magnitude range (none if equal):"); */
/*   //  mag1=readf(mag1);mag2=readf(mag2); */
/*   //  printf(" Input factor for box: "); */
/*   //  fact=readf(fact); */
  *nstar=apmread(alfacen,deltacen,xdim/3600./2.*factra,ydim/3600./2.*factdec,0.,WCS_J2000,2000.,epoch,mag1,mag2,-1,nmax,gnum,arl,decl,gmag,gtype,verbose);
  printf("paso apmread\n");
/*   printf(" File with APM file: "); */
/*   scanf("%s",fileapm); */
/*   *nstar=apmread_file(fileapm,alfac,deltac,xdim/3600./2.*fact,ydim/3600./2.*fact,0.,1,2000.,2000.,mag1,mag2,-1,NMAX,gnum,arl,decl,gmag,gtype,0); */
  printf(" Found stars: %d\n",*nstar);

  for(i=0;i<*nstar;i++) {
    todo[i].num= gnum[i];
    todo[i].ar =  arl[i];
    todo[i].dec= decl[i];
    todo[i].mag= gmag[i];
    todo[i].mag2=gtype[i];
    
  }
  qsort(todo,*nstar,sizeof(struct star_table),iscompare);

  for(i=0;i<*nstar;i++) {
    arl[i] =  todo[i].ar;
    decl[i]=  todo[i].dec; 
    gmag[i]=  todo[i].mag;
    gtype[i]= todo[i].mag2;
    gnum[i]=  todo[i].num;
/*     //printf(" Star %8.0f AR %g Dec %g  Mag %f\n",todo[i].num,todo[i].ar,todo[i].dec,todo[i].mag); */

  }



  strcpy(fileref,inputfile);
  strcat(fileref,".apm");
      
  if((fref=fopen(fileref,"w"))==NULL) {
    printf("\nERROR: Can't open file %s\n",fileref);
    exit(1);
  }
  printf(" Saving reference stars used for astrometry in file %s\n",fileref);
  fprintf(fref,"#Reference stars in the searched file from the APM catalogue\n");
  fprintf(fref,"#round %f hours %f degrees, box %f x %f arcmin for epoch %f\n",alfacen/15,deltacen,xdim/60.*factra,ydim/60.*factdec,epoch);
  fprintf(fref,"#Num     RA(J2000)   DEC(J2000)    Red_mag  Obj_type\n");
  for(i=0;i<*nstar;i++) {
    ar[i]=arl[i];dec[i]=decl[i];mag[i]=gmag[i];
/*     printf(" Star %8.0f AR %g Dec %g\n",gnum[i],ar[i],dec[i]); */
    ra2str(snul1,32,ar[i],2);
    dec2str(snul2,32,dec[i],1);
    fprintf(fref," %8.0f %s %s %6.2f  %8d\n",gnum[i],snul1,snul2,gmag[i],gtype[i]);
  }
  fclose(fref);


}




void Plot(float *x,float *y,int n,int mark,int color)
{
  cpgsci(color);
  cpgpt(n,x,y,mark);
  cpgsci(1);
}


void Plotd(double *x,double *y,int n,int mark,int color)
{
  int i;
  cpgsci(color);
  for(i=0;i<n;i++)  cpgpt1(x[i],y[i],mark);
  cpgsci(1);
}



void PlotStars(double *ar, double *dec, double *mag, int nstar,double *alfa,double *delta)
{
  int i;
  float armin,armax,decmin,decmax;
  double magmin,magmax;
  float ch;
  float x1,x2,y1,y2;
  float xtemp;

  cpgqch(&ch);
  
  cpgsvp(0.6,0.9,0.6,0.9);
  MinMax_d(nstar,mag,&magmin,&magmax);

  pgLimits_d(nstar,ar,&armin,&armax);
  pgLimits_d(nstar,dec,&decmin,&decmax);
  if(rotang==0 || rotang==180) {
    if(rotang ==0) {
      x1=armax;x2=armin;y1=decmin;y2=decmax;
    }
    else if(rotang==180) {
      x1=armin;x2=armax;y1=decmax;y2=decmin;
    }
    else {
      x1=armax;x2=armin;y1=decmin;y2=decmax;
    }
    if(flip) {
      xtemp=x1;
      x1=x2;
      x2=xtemp;
    }



/*     //Esto es para que los ejes esten a la misma escala */
    cpgwnad(x1*3600.,x2*3600.,y1*3600.,y2*3600.); 
/*     //cpgeras(); */
/*     //Ar va en grados, luego para pasarlo a segundos de tiempo, multiplico  */
/*     // por 3600 y divido por 15. La declinacion la pasao a arcsec */
    cpgswin(x1*3600/15.,x2*3600/15.,y1*3600,y2*3600);
 
    cpgtbox("ZXBCTNSH",0.0,0,"ZBCTMSDV",0.0,0);
    cpglab("RA","Dec","Reference stars");
/*     //Ahora paso a las coordenadas que  tengo que pintar otra vez */
    cpgswin(x1,x2,y1,y2);
    for( i=0;i< nstar;i++) { 
      cpgsch((float)(0.5+(magmax-mag[i])*4./(magmax-magmin)));
/*       //printf(" SIE %f\n",(float)(0.5+mag[i]*4./(magmax-magmin))); */
      cpgpt1(ar[i],dec[i],2);
    }

    cpgsci(4);
    for( i=0;i< ncross;i++) { 
      cpgpt1(alfa[i],delta[i],4);
    }
    cpgsci(1);
  }
  else {
    if(rotang ==90) {
      x1=decmax;x2=decmin;y1=armax;y2=armin;
    }
    else if(rotang==-90 || rotang == 270) {
      x1=decmin;x2=decmax;y1=armin;y2=armax;
/*       printf(" ESTOY AQUI!!!! %f %f %f %f\n",x1,x2,y1,y2); */
    }
    else {
      x1=decmin;x2=decmax;y1=armin;y2=armax;
    }
/*     printf(" x1 %f x2 %f\n",x1,x2); */
    if(flip) {
      xtemp=x1;
      x1=x2;
      x2=xtemp;
    }
/*     printf(" x1 %f x2 %f\n",x1,x2); */
    cpgwnad(x1*3600,x2*3600,y1*3600,y2*3600);
    cpgswin(x1*3600,x2*3600,y1*3600/15,y2*3600/15);
    cpgtbox("ZBCTNSD",0.0,0,"ZXBCTNSHV",0.0,0);
    cpglab("Dec","RA","Reference stars");
    cpgswin(x1,x2,y1,y2);
    for( i=0;i< nstar;i++) {
      cpgsch((float)(0.5+(magmax-mag[i])*4./(magmax-magmin)));
/*       //printf(" SIE %f mag %f\n",(float)(0.5+(magmax-mag[i])*4./(magmax-magmin)),mag[i]); */
/*       //printf(" delta %f alfa %f\n",dec[i],ar[i]); */
      cpgpt1(dec[i],ar[i],2);
    }
    cpgsci(4);
    for( i=0;i< ncross;i++) {
      cpgpt1(delta[i],alfa[i],4);
    }
    
    cpgsci(1);
  }
  cpgsch(ch);

  
}





void ClickBoth(double *ar    ,double *dec   ,int nstar,
               double *xpoint,double *ypoint,int npoint,
               double *alfa  ,double *delta ,
               double *xpix,  double *ypix  ,
               int naxes0,int naxes1)
{
  char option='A';
  float xcur,ycur;
  int istar,ipoint=0,j;
  char cnul;
  float armin,armax,decmin,decmax;
  float mindist=1.e15,dist;
  float x1,x2,y1,y2;
  float xtemp;
  pgLimits_d(nstar,ar,&armin,&armax);
  pgLimits_d(nstar,dec,&decmin,&decmax);
  



  while(option!='E' && option!='W' && option!='e' && option!= 'w') {
    
    
    printf(" Click on image...\n");
    cpgsvp(0.1,0.5,0.1,0.9);
/*     //cpglab("pixel","pixel",""); */
    cpgwnad(0.,naxes0,0.,naxes1);
    cpgcurs(&xcur,&ycur,&cnul);
    mindist=1.e15;
    for(j=0;j<npoint;j++) {
      dist=((xpoint[j]-xcur)*(xpoint[j]-xcur)+(ypoint[j]-ycur)*(ypoint[j]-ycur));
/*       printf(" Is %d, di %f min %f\n",ipoint,dist,mindist); */
      if(mindist*mindist>dist*dist) {
	mindist=dist;
	ipoint=j;
      }
    }
    printf(" Cursor position %f %f\n",xcur,ycur);
    printf(" Nearest object %d at position %f %f\n",ipoint,xpoint[ipoint],ypoint[ipoint]);
    cpgsci(5);
    cpgsch(3.);
    cpgpt1((float)xpoint[ipoint],(float)ypoint[ipoint],3);
    cpgsci(1);
    cpgsch(1.);
    
    printf(" Click on ar-dec plane...\n");
    cpgsvp(0.6,0.9,0.6,0.9);
    /*     Esto es para que los ejes esten a la misma escala */
    if(rotang==0 || rotang==180) {
      if(rotang ==0) {
	x1=armax;x2=armin;y1=decmin;y2=decmax;
      }
      else if(rotang==180) {
	x1=armin;x2=armax;y1=decmax;y2=decmin;
      }
      else {
	x1=armax;x2=armin;y1=decmin;y2=decmax;
      }
      if(flip) {
	xtemp=x1;
	x1=x2;
	x2=xtemp;
      }
      
      cpgwnad(x1*3600,x2*3600,y1*3600,y2*3600); 
      cpgswin(x1,x2,y1,y2);
/*       cpgswin(x1*3600/15,x2*3600/15,y1*3600,y2*3600); */
    }
    else {
      if(rotang ==90) {
	x1=decmax;x2=decmin;y1=armax;y2=armin;
      }
      else if(rotang==-90 || rotang == 270) {
	x1=decmin;x2=decmax;y1=armin;y2=armax;
/* 	printf(" ESTOY DENTROOOOO!!! %f %f %f %f\n",x1,x2,y1,y2); */
      }
      else {
	x1=decmin;x2=decmax;y1=armin;y2=armax;
      }
      if(flip) {
	xtemp=x1;
	x1=x2;
	x2=xtemp;
      }
      cpgwnad(x1*3600,x2*3600,y1*3600,y2*3600);
      cpgswin(x1,x2,y1,y2);
/*     cpgswin(x1*3600,x2*3600,y1*3600/15,y2*3600/15); */
    }
    mindist=1.e15;
    
    cpgcurs(&xcur,&ycur,&cnul);
/*     printf(" Cursor position %f %f \n",xcur,ycur); */


    istar=-1;
    for(j=0;j<nstar;j++) {
      if(rotang==0 || rotang==180)  dist=(((float)ar[j]-xcur)*((float)ar[j]-xcur)+((float)dec[j]-ycur)*((float)dec[j]-ycur));
      else  dist=(((float)ar[j]-ycur)*((float)ar[j]-ycur)+((float)dec[j]-xcur)*((float)dec[j]-xcur));
  
/*       printf(" Is %d,n %d di %f min %f\n",istar,j,dist,mindist); */
      if(mindist*mindist>dist*dist) {
	mindist=dist;
	istar=j;
      }
    }
    printf(" Cursor position %f %f \n",xcur,ycur);
    printf(" Nearest object %d at position %f %f\n",istar,ar[istar],dec[istar]);
    cpgsci(5);
    cpgsch(3.);
    if(rotang==0 || rotang==180)  cpgpt1(ar[istar],dec[istar],3);
    else  cpgpt1(dec[istar],ar[istar],3);
    cpgsci(1);
    cpgsch(1.);
    
    printf("  A Accept and continue\n");
    printf("  N Do not accept and continue\n");
    printf("  W Accept and exit\n");
    printf("  E Do not accept and exit\n");
    
    option=readc('A');
    switch (option) {
    case 'A':
    case 'a':
      alfa[ncross]=ar[istar];delta[ncross]=dec[istar];
      xpix[ncross]=xpoint[ipoint];ypix[ncross]=ypoint[ipoint];
      ncross++;
      break;
    case 'W':
    case 'w':
      alfa[ncross]=ar[istar];delta[ncross]=dec[istar];
      xpix[ncross]=xpoint[ipoint];ypix[ncross]=ypoint[ipoint];
      ncross++;
      break;
    } 
    
    
  }
  
}

void ClickImFindStar(double *ar    ,double *dec   ,int nstar,
                     double *alfa  ,double *delta ,
                     double *xpix,  double *ypix  ,
                     int naxes0,int naxes1)
{
  char option='A';
  float xcur,ycur;
  int istar,j;
  char cnul;
  float armin,armax,decmin,decmax;
  float mindist=1.e15,dist;
  double xpos,ypos; 

  pgLimits_d(nstar,ar,&armin,&armax);
  pgLimits_d(nstar,dec,&decmin,&decmax);
  



  while(option!='E' && option!='W' && option!='e' && option!= 'w') {
    
    
    printf(" Click on image...\n");
    cpgsvp(0.1,0.5,0.1,0.9);
/*     //cpglab("pixel","pixel",""); */
    cpgwnad(0.,naxes0,0.,naxes1);
    cpgcurs(&xcur,&ycur,&cnul);
    xpos=(double)xcur;ypos=(double)ycur; 
    printf(" Cursor position %f %f\n",xcur,ycur);
/*     printf(" Nearest object position %f %f\n",xpoint[ipoint],ypoint[ipoint]); */
    cpgsci(5);
    cpgsch(3.);
    cpgpt1(xcur,ycur,3);
    cpgsci(1);
    cpgsch(1.);
    
    printf(" Click on ar-dec plane...\n");
    cpgsvp(0.6,0.9,0.6,0.9);
/*     //Esto es para que los ejes esten a la misma escala */
    if(rotang==0 || rotang==180) {
      cpgwnad(armax*3600,armin*3600,decmin*3600,decmax*3600); 
      cpgswin(armax,armin,decmin,decmax);
    }
    else {
/*       //cpgwnad(decmax*3600,decmin*3600,armin*3600,armax*3600); */
/*       //cpgswin(decmax,decmin,armin,armax); */
      cpgwnad(decmin*3600,decmax*3600,armin*3600,armax*3600);
      cpgswin(decmin,decmax,armin,armax);
    }
    mindist=1.e15;
    
    cpgcurs(&xcur,&ycur,&cnul);
    printf(" Cursor position %f %f \n",xcur,ycur);


    istar=-1;
    for(j=0;j<nstar;j++) {
      if(rotang==0 || rotang==180)  dist=(((float)ar[j]-xcur)*((float)ar[j]-xcur)+((float)dec[j]-ycur)*((float)dec[j]-ycur));
      else  dist=(((float)ar[j]-ycur)*((float)ar[j]-ycur)+((float)dec[j]-xcur)*((float)dec[j]-xcur));
/*       printf(" Is %d,n %d di %f min %f\n",istar,j,dist,mindist); */
      if(mindist*mindist>dist*dist) {
	mindist=dist;
	istar=j;
      }
    }
    printf(" Cursor position %f %f la i %d\n",xcur,ycur,istar);
    printf(" Nearest object position %f %f\n",ar[istar],dec[istar]);
    cpgsci(5);
    cpgsch(3.);
    if(rotang==0 || rotang==180)  cpgpt1(ar[istar],dec[istar],3);
    else  cpgpt1(dec[istar],ar[istar],3);
    cpgsci(1);
    cpgsch(1.);
    
    printf("  A Accept and continue\n");
    printf("  N Do not accept and continue\n");
    printf("  W Accept and exit\n");
    printf("  E Do not accept and exit\n");
    
    option=readc('A');
    switch (option) {
    case 'A':
    case 'a':
      alfa[ncross]=ar[istar];delta[ncross]=dec[istar];
      xpix[ncross]=xpos;ypix[ncross]=ypos;
      ncross++;
      break;
    case 'W':
    case 'w':
      alfa[ncross]=ar[istar];delta[ncross]=dec[istar];
      xpix[ncross]=xpos;ypix[ncross]=ypos;
      ncross++;
      break;
    }
    
    
  }
  
}


int ComputeSolution(double *alfa,double *delta,
                     double *xpix,double *ypix)
{
  float xmin,xmax,ymin,ymax;
  int i,j;
  static float factorarr=300;
  char snul[1000];
  double dnul1,dnul2;
  double xpixcal[NMAX],ypixcal[NMAX];   /* //De alfa,delta a xpixcal,ypixcal */
  double alfacal[NMAX],deltacal[NMAX];  /* //De xpix,ypix a alfacal, deltacal */
  float xerror[NMAX],yerror[NMAX],rerror[NMAX]; /* //Diferencia xpix-xpixcal... */
  float aerror[NMAX],derror[NMAX];              /* //Diferencia alfa-alfacal... */
  float tr[6];
  char option='A';
  
/*   //Variables for error computation */
  struct WorldCoor *wcserr[1000];      /* World coordinate system structure */
  int nbootstrap;
  double alfaboot[NMAX],deltaboot[NMAX];
  double xpixboot[NMAX], ypixboot[NMAX];
  double alfaerr,deltaerr,xpixerr,ypixerr;
/*   //char *errheader; */
/*   int lhead,nbfits; */
  float alfacalNE[1000],deltacalNE[1000];
  float alfacalNW[1000],deltacalNW[1000];
  float alfacalSE[1000],deltacalSE[1000];
  float alfacalSW[1000],deltacalSW[1000];
  float alfacalC[1000], deltacalC[1000];
  float first,third,median;

  int verbose;

  int coe1,coe2;
  double cofi[30];
  int pas=1;

  int off;

  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  printf(" Computing calibration with %d stars \n",ncross); 
  setdefwcs(0);
  verbose=interact;


  if(interact) {
    printf(" Degree of polynomial to fit: ");
    ngrad=readf(ngrad);
  }
  else printf(" Using polynomial degree %d\n",ngrad);

  if(ncross>ngrad) FitPlate(wcsim,xpix,ypix,alfa,delta,ncross,ngrad,verbose); 
  printf(" Solution computed. Plotting\n");
  while(ncross>0 && option!='E') { 
    if(interact || (pas==2 && option=='F') ) {
      cpgeras();
      cpgsvp(0.1,0.5,0.1,0.9);
      cpgwnad(0.,naxes[0],0.,naxes[1]);
      if(interact)  
	cpggray(buffer,naxes[0],naxes[1],1,naxes[0],1,naxes[1],datamax,datamin,tr);
      cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
      strcpy(snul,"Image \0");
      strcat(snul,inputfile);
      cpglab("X axis","Y axis",snul);
      if(DEBUG) printf(" pinto xpix, ypix\n");
      if(interact)  
	Plotd(xpix,ypix,ncross,2,2);
      if(DEBUG) printf(" ya he pinto xpix, ypix\n");
    }

/*     printf(" tipo de proyeccion: %d\n",(*wcsim).prjcode);  */
/*     printf(" crval1 %f crde1 %f crpix %f\n",(*wcsim).xref,(*wcsim).xinc,(*wcsim).xrefpix);  */
/*     printf(" crval2 %f crde2 %f crpix %f\n",(*wcsim).yref,(*wcsim).yinc,(*wcsim).yrefpix); */
/*     Esto es una prueba lñinela xpix[0]= 90.2749;ypix[0]= 525.3647; */
    
    if(DEBUG) printf(" Pasa por aqui\n"); 
    GetPlate (wcsim, &coe1, &coe2, cofi);

/*     printf(" ceo1 %d coe2 %d\n",coe1,coe2); */
/*     for(i=0;i<30;i++) printf(" coe %d:  %g\n",i,cofi[i]); */
/*     //i=readi(1); */
/*     FitPlate(wcsim,xpix,ypix,alfa,delta,ncross,ngrad,0); */
/*     //TNXFit(wcsim,xpix,ypix,alfa,delta,ncross,1,1); */
/*     printf(" Y por aqui\n"); */
    GetPlate (wcsim, &coe1, &coe2, cofi);
/*     //printf(" ceo1 %d coe2 %d\n",coe1,coe2); */
/*     //for(i=0;i<30;i++) printf(" coe %d:  %f\n",i,cofi[i]); */

 
    
    for(i=0;i<ncross;i++) {
      
      if(DEBUG) printf(" Em 1 alfa[i] %g delta[i] %g\n",alfa[i],delta[i]); 
      if(DEBUG) printf(" wcs->radecin %s\n",wcsim->radecin);
      wcs2pix(wcsim,alfa[i], delta[i], &dnul1, &dnul2,&off); 
      if(DEBUG) printf(" Cal 1\n");
      xpixcal[i]=dnul1;ypixcal[i]=dnul2;
      if(DEBUG) printf(" Em 2\n"); 
      pix2wcs(wcsim,xpix[i], ypix[i] , &dnul1, &dnul2);
      if(DEBUG) printf(" Cal 2\n");
      alfacal[i]=dnul1;deltacal[i]=dnul2;
      if(DEBUG) printf(" Cal 3\n");
/*       printf(" xpix %f ypix %f alfa %f delta %f alfacal %f deltacal %f\n",xpix[i],ypix[i],alfa[i],delta[i],alfacal[i],deltacal[i]); */
      
      xerror[i]=(xpix[i]-xpixcal[i]);
      if(DEBUG) printf(" Bien 1\n"); 
      yerror[i]=(ypix[i]-ypixcal[i]);
      if(DEBUG) printf(" Bien 2\n"); 
      rerror[i]=sqrt((xpix[i]-xpixcal[i])*(xpix[i]-xpixcal[i])+(ypix[i]-ypixcal[i])*(ypix[i]-ypixcal[i]));
      if(DEBUG) printf(" Bien 3\n"); 
      aerror[i]=alfa[i]-alfacal[i];
      if(DEBUG) printf(" Bien 4\n"); 
      derror[i]=delta[i]-deltacal[i];
      if(DEBUG) printf(" Bien 5 \n"); 
      if(factorarr!=0 && (interact || (pas==2 && option=='F') )) {
	cpgsci(3);
	cpgsvp(0.1,0.5,0.1,0.9);
	/* 	cpglab("pixel","pixel",""); */
/* 	//printf(" ein ka\n"); */
	cpgwnad(0.,naxes[0],0.,naxes[1]);
	/* 	printf("ENTRO\n"); */
/* 	printf(" BIEn e  %f %f %f %f  %f\n",xpix[i],ypix[i],xpixcal[i],ypixcal[i],factorarr); */
	cpgarro((float)(xpix[i]),(float)(ypix[i]),(float)(xpix[i]-(xpix[i]-xpixcal[i])*factorarr),(float)(ypix[i]-(ypix[i]-ypixcal[i])*factorarr));
/* 	//	  printf("error in arcseconds %f\n",error[i]); */
      }
      
    }
    if(DEBUG) printf(" YA calculo los errores\n"); 
    aem=StMedia(ncross,aerror,&aes);
    dem=StMedia(ncross,derror,&des);
    xem=StMedia(ncross,xerror,&xes);
    yem=StMedia(ncross,yerror,&yes);
    printf(" Mean AR_real - AR_calculated: %g with standard deviation %g\n",aem,aes);
    printf(" Mean DEC_real - DEC_calculated: %g with standard deviation %g\n",dem,des);
    printf(" Mean X_real - X_calculated: %g with standard deviation %g\n",xem,xes);
    printf(" Mean Y_real - Y_calculated: %g with standard deviation %g\n",yem,yes);
    
    if(interact|| (pas==2 && option=='F') ) {
      cpgsci(1);
      pgLimits(ncross,xerror,&xmin,&xmax);
      pgLimits(ncross,yerror,&ymin,&ymax);
      cpgsvp(0.6,0.9,0.1,0.5);
      cpgwnad(xmin,xmax,ymin,ymax);
      cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
      cpgpt(ncross,xerror,yerror,1);
      cpglab("X pixel","Y pixel","");
    }
    
    printf(" A Change arrows factor\n");
    printf(" F Fit plate solution again\n");
/*     printf(" S Shift solution using secondary stars\n"); */
    printf(" P Change polynomial degree\n");
    printf(" I Delete reference point in image\n");
    printf(" D Delete reference point in error plot\n"); 
    printf(" C Automatic sigma-clipping algorithm to delete points\n"); 
    printf(" R Compute error in astrometric coeficients\n");
    printf(" E Astrometry is good! and exit\n");
    if(!interact) {
/*       printf(" Entra opciones!!!\n"); */
      if(option=='A') option='C';
      else if(option=='C' && pas==1) option='F';
      else if(option=='F' && pas==1) {
	pas=2;
	option='C';
      }
      else if(option=='C' && pas==2) option='F';
      else if(option=='F' && pas==2) option='E';
/*       printf(" Sale con %c\n",option); */
    }
    else option=readc('F');
/*     //option=readc('E');   */
    switch (option) { 
    case 'F':
    case 'f':
      FitPlate(wcsim,xpix,ypix,alfa,delta,ncross,ngrad,0);
      break;
/*     case 'S': */
/*     case 's': */
/*       PlotStars(ar,dec,mag,nstar,alfa,delta); */
/*       ShiftSolution(alfas,deltas,xpixs,ypixs); */
/*       printf(" Changing central AR DEC from %f %f\n",alfac,deltac); */
/*       alfac=wcsim->xref;deltac=wcsim->yref; */
/*       pix2wcs(wcsim,naxes[0]/2., naxes[1]/2. ,&dnul1, &dnul2); */
/*       alfac=(float)dnul1;deltac=(float)dnul2; */
/*       xpixsiz=fabs((float)((*wcsim).xinc)*3600.); */
/*       ypixsiz=fabs((float)((*wcsim).yinc)*3600.); */
/*       printf(" to %f %f\n",alfac,deltac); */


/*       break; */
    case 'R':
    case 'r':
      errflag=1;
      printf(" Input number of bootstrapping solutions (max=1000): ");
      nbootstrap=readi(100);
      printf(" Input error of reference stars in alfa  coordinate (arcsec) : ");
      alfaerr=(double)readf(0.);
      printf(" Input error of reference stars in delta coordinate (arcsec) : ");
      deltaerr=(double)readf(0.);
      printf(" Input error of  reference points in X coordinate (pixel) : ");
      xpixerr=(double)readf(0.);
      printf(" Input error of  reference points in Y coordinate (pixel) : ");
      ypixerr=(double)readf(0.);

      for(i=0;i<nbootstrap;i++) {
	for(j=0;j<ncross;j++) {
	  alfaboot[j]=alfa[j]+Gasdev()*alfaerr/3600.;
	  deltaboot[j]=delta[j]+Gasdev()*deltaerr/3600.;
	  xpixboot[j]=xpix[j]+Gasdev()*xpixerr;
	  ypixboot[j]=ypix[j]+Gasdev()*ypixerr;
	} 
/* 	//printf(" Hasta aqui bien\n"); */
	wcserr[i]=(struct WorldCoor *) calloc (1, sizeof(struct WorldCoor));

	memcpy(wcserr[i],wcsim,sizeof(struct WorldCoor));
	FitPlate(wcserr[i],xpixboot,ypixboot,alfaboot,deltaboot,ncross,ngrad,0);
/* 	//printf(" Hasta aqui bien mejor\n"); */
/* 	//if(iswcs(wcserr[i])) printf(" WCS information set\n"); */
/* 	//else printf(" No WCS imformation set. NO WRITING\n"); */
/* 	//memcpy(errheader,header,lhead); */
/* 	//free(errheader); */
/* 	//if ((errheader = fitsrhead (inputfile, &lhead, &nbfits)) == NULL) { */
/* 	//  fprintf (stderr, "Cannot read FITS header of file %s\n", inputfile); */
/* 	//  exit(1); */
/* 	//} */
/* 	//SetFITSWCS(errheader,wcserr[i]); */
	printf(" Iteration number %d \n",i);
/* 	//PrintWCS(errheader,1); */
	pix2wcs(wcserr[i],1., (float)naxes[1] , &dnul1, &dnul2);
	alfacalNW[i]=(float)dnul1;deltacalNW[i]=(float)dnul2;
	pix2wcs(wcserr[i], (float)naxes[0] ,1., &dnul1, &dnul2);
	alfacalSE[i]=(float)dnul1;deltacalSE[i]=(float)dnul2;
	pix2wcs(wcserr[i],1., 1. , &dnul1, &dnul2);
	alfacalSW[i]=(float)dnul1;deltacalSW[i]=(float)dnul2;
	pix2wcs(wcserr[i],(float)naxes[0], (float)naxes[1] , &dnul1, &dnul2);
	alfacalNE[i]=(float)dnul1;deltacalNE[i]=(float)dnul2;
	pix2wcs(wcserr[i],(float)naxes[0]/2., (float)naxes[1]/2. , &dnul1, &dnul2);
	alfacalC[i]=(float)dnul1;deltacalC[i]=(float)dnul2;

	

      }
      aem=StMedia(nbootstrap, alfacalNW,&aes);
      dem=StMedia(nbootstrap,deltacalNW,&des);
      Quartil(nbootstrap,alfacalNW,&first,&median,&third);
      aes=(third-first)/1.35;
      Quartil(nbootstrap,deltacalNW,&first,&median,&third);
      des=(third-first)/1.35;
      printf(" Estimated error for upleft corner: %f in alfa %f in delta\n",aes*3600., des*3600.);
      aem=StMedia(nbootstrap, alfacalNE,&aes);
      dem=StMedia(nbootstrap,deltacalNE,&des);
      Quartil(nbootstrap,alfacalNE,&first,&median,&third);
      aes=(third-first)/1.35;
      Quartil(nbootstrap,deltacalNE,&first,&median,&third);
      des=(third-first)/1.35   ;
      printf(" Estimated error for upright corner: %f in alfa %f in delta\n",aes*3600., des*3600.);
      aem=StMedia(nbootstrap, alfacalSW,&aes);
      dem=StMedia(nbootstrap,deltacalSW,&des);
      Quartil(nbootstrap,alfacalSW,&first,&median,&third);
      aes=(third-first)/1.35;
      Quartil(nbootstrap,deltacalSW,&first,&median,&third);
      des=(third-first)/1.35;
      printf(" Estimated error for downleft corner: %f in alfa %f in delta\n",aes*3600., des*3600.);
      aem=StMedia(nbootstrap, alfacalSE,&aes);
      dem=StMedia(nbootstrap,deltacalSE,&des);
      Quartil(nbootstrap,alfacalSE,&first,&median,&third);
      aes=(third-first)/1.35;
      Quartil(nbootstrap,deltacalSE,&first,&median,&third);
      des=(third-first)/1.35;
      printf(" Estimated error for downright corner: %f in alfa %f in delta\n",aes*3600., des*3600.);
      aem=StMedia(nbootstrap, alfacalC,&aes);
      dem=StMedia(nbootstrap,deltacalC,&des);
/*       //printf(" Estimated error for center of image: %f in alfa %f in delta\n",aes*3600., des*3600.); */
      Quartil(nbootstrap,alfacalC,&first,&median,&third);
      aes=(third-first)/1.35;
      aem=median;
      Quartil(nbootstrap,deltacalC,&first,&median,&third);
      des=(third-first)/1.35;
      dem=median;
      printf(" Estimated error for center of image: %f in alfa %f in delta\n",aes*3600., des*3600.);

      cpgsvp(0.6,0.9,0.6,0.9);
      cpgswin((float)(-aes*3.*3600.),(float)(aes*3.*3600.),0.,nbootstrap/5.);
      cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
      cpgswin((float)(aem-aes*3.),(float)(aem+aes*3.),0.,nbootstrap/5.);
      cpghist(nbootstrap,alfacalC,(float)(aem-aes*3.),(float)(aem+aes*3.),20,1);
      cpgsci(2);
      cpgswin((float)(-des*3.*3600.),(float)(des*3.*3600.),0.,nbootstrap/5.);
      cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
      cpgswin((float)(dem-des*3.),(float)(dem+des*3.),0.,nbootstrap/5.);
      cpghist(nbootstrap,deltacalC,(float)(dem-des*3.),(float)(dem+des*3.),20,1);
      cpgsci(1);
      printf(" Everithing in arcseconds\n");

/*       i=readi(1); */
      break;
    case 'P':
    case 'p':
      printf(" Degree of polynomial to fit: ");
      ngrad=readf(ngrad);
      break;
    case 'A':
    case 'a':
      printf("Enter another factor for arrows :");
      factorarr=readf(0);
      break;
    case 'I':
    case 'i':
      DeleteClickIm(alfa,delta,xpix,ypix);
      
      break;
    case 'D':
    case 'd':
      DeleteClickErr(alfa,delta,xpix,ypix,xerror,yerror);
      break;
    case 'C':
    case 'c':
/*       printf(" Entra C\n"); */
      fflush(stdout);
      fflush(stderr);
      setvbuf(stdin,"",_IOLBF,0);
      setvbuf(stdout,"",_IOLBF,0);

/*       //printf(" ncross %d\n",ncross); */
      for(j=0;j<ncross;j++) {
/* 	//printf("%d  ERROR  %f %f   %f %f\n",j,(xerror[j]-xem)/xes,(yerror[j]-yem)/yes,xerror[j],yerror[j]); */
	if(abs(xerror[j]-xem)/xes> 1.8  || abs(yerror[j]-yem)/yes> 1.8) {
	  printf(" Deleting object number %d out of %d\n",j,ncross);	     
	  printf(" xerror[j]-xem %f yerror[j]-yem %f \n",xerror[j]-xem,yerror[j]-yem);
	  memmove(  xpix+j,xpix+j+1  ,((ncross)-j-1)*sizeof(double));
	  memmove(  ypix+j,ypix+j+1  ,((ncross)-j-1)*sizeof(double));
	  memmove(  alfa+j,alfa+j+1  ,((ncross)-j-1)*sizeof(double));
	  memmove( delta+j,delta+j+1 ,((ncross)-j-1)*sizeof(double));
	  memmove(xerror+j,xerror+j+1,((ncross)-j-1)*sizeof(float));
	  memmove(yerror+j,yerror+j+1,((ncross)-j-1)*sizeof(float));
	  ncross--;
	  j--;
	}
      }
/*       j=readi(j); */
/*       //exit(1); */
      break;
    case 'E':
    case 'e':
      printf(" %d cross-corelated objects have been used for astrometric solution\n",ncross);
      break;
    }
  }
  
  /*   printf(" 222 Ncorss vale  %d \n",ncross); */
  
  if((ncross)==0)  { 
    printf(" Not enough points to fit solution\n");
    return(0);
  }



  return(1);
}





void WriteWCS_keys() {
  
  double  ra, dec; 
  char wcstemp[16];
  char comment[50];
  static char wcsproj[8]="TAN";           /* WCS projection name */
  char snul[FLEN_VALUE];
  int ok1,ok2;
  double dnul;
  
  if(iswcs(wcsim)) ;
  else {
    printf(" No WCS imformation set\n");
    return;
  }
  

  
  
    /* Rename old center coordinates */
/*     if (!ksearch (header,"WRA") && ksearch (header,"RA")) */
/*         hchange (header,"RA","WRA"); */
/*     if (!ksearch (header,"WDEC") && ksearch (header,"DEC")) */
/*         hchange (header,"DEC","WDEC"); */

/*     if (!ksearch (header,"WEQUINOX") && ksearch (header,"EQUINOX")) */
/*         hchange (header, "EQUINOX", "WEQUINOX"); */


    ok1=ffgcrd(fitsimage,"WRA",snul,&status);status=0; 
    ok2=ffgcrd(fitsimage,"RA",snul,&status); 
    if(ok1==202 &&  ok2==0)        fits_modify_name(fitsimage,"RA","WRA",&status); 
    ok1=ffgcrd(fitsimage,"WDEC",snul,&status);status=0;
    ok2=ffgcrd(fitsimage,"DEC",snul,&status);
    if(ok1==202 && ok2==0)  fits_modify_name(fitsimage,"DEC","WDEC",&status);
    ok1=ffgcrd(fitsimage,"WEQUINOX",snul,&status);status=0;
    ok2=ffgcrd(fitsimage,"EQUINOX",snul,&status);
    if(ok1==202 && ok2==0) fits_modify_name(fitsimage,"EQUINOX","WEQUINOX",&status); 
    


    /* Only change EPOCH if it is used instead of EQUINOX */
/*     else if (!ksearch (header,"WEPOCH") && ksearch (header,"EPOCH")) */
/*         hchange (header, "EPOCH", "WEPOCH"); */

    else {
      ok1=ffgcrd(fitsimage,"WEPOCH",snul,&status);status=0;
      ok2=ffgcrd(fitsimage,"EPOCH",snul,&status);
      if(ok1==202 && ok2==0)   fits_modify_name(fitsimage,"EPOCH","WEPOCH",&status);
    }

    strcpy(comment,"Comparision catalogue");
    strcpy(wcstemp,astcat);
    status=0;
    fits_delete_key(fitsimage, "ASTCAT", &status); status=0;
    fits_write_key(fitsimage,TSTRING,"ASTCAT",wcstemp,comment,&status);
    strcpy(comment,"Number of primary stars used");
    fits_delete_key(fitsimage, "NCROSS", &status); status=0;
    fits_write_key(fitsimage,TINT,"NCROSS",&ncross,comment,&status);
    if(secflag) {
      strcpy(comment,"Secondary   catalogue");
      strcpy(wcstemp,seccat);
      status=0;
      fits_delete_key(fitsimage, "SECCAT", &status); status=0;
      fits_write_key(fitsimage,TSTRING,"SECCAT",wcstemp,comment,&status);
      strcpy(comment,"Number of secondary stars used");
      fits_delete_key(fitsimage, "NCROSSEC", &status); status=0;
      fits_write_key(fitsimage,TINT,"NCROSSEC",&ncrosssec,comment,&status);
      strcpy(comment,"Shift in RA applied from secondary stars ('')");
      dnul=shiftra*3600.;
      fits_delete_key(fitsimage, "SHIFTRA", &status); status=0;
      fits_write_key(fitsimage,TDOUBLE,"SHIFTRA",&(dnul),comment,&status);
      strcpy(comment,"Shift in DEC applied from secondary stars ('')");
      dnul=shiftdec*3600.;
      fits_delete_key(fitsimage, "SHIFTDEC", &status); status=0;
      fits_write_key(fitsimage,TDOUBLE,"SHIFTDEC",&(dnul),comment,&status);

    }


    status=0;
    strcpy(comment,"Standard dev of solution in X");
    fits_delete_key(fitsimage, "STDX", &status); status=0;
    fits_write_key(fitsimage,TFLOAT,"STDX",&xes,comment,&status);
    strcpy(comment,"Standard dev of solution in Y");
    fits_delete_key(fitsimage, "STDY", &status); status=0;
    fits_write_key(fitsimage,TFLOAT,"STDY",&yes,comment,&status);


/*     printf(" Rtrneia \n"); */

    /* Set new center coordinates */
/*     hputra (header,"RA",wcsim->xref); */
/*     hputdec (header,"DEC",wcsim->yref); */
/*     hputr8 (header, "EQUINOX", wcsim->equinox); */
/*     if (hgetr8 (header, "WEPOCH", &ep)) */
/*         hputr8 (header, "EPOCH", wcsim->equinox); */
/*     else if (!hgetr8 (header, "EPOCH", &ep)) */
/*         hputr8 (header, "EPOCH", wcsim->equinox); */
/*     hputs (header, "RADECSYS", wcsim->radecsys); */


/*     //QUITARLEAS */
/*     //if(status) fits_report_error(stderr,status); */
/*     //if(status) fits_report_error(stderr,status); */
/*     //return; */
    pix2wcs(wcsim,naxes[0]/2., naxes[1]/2. ,&ra, &dec);
    ra2str(wcstemp,16,ra,3);
    strcpy(comment,"Center RA");
    status=0;
    fits_write_key(fitsimage,TSTRING,"RA",wcstemp,comment,&status);
    dec2str(wcstemp,16,dec,3);
    strcpy(comment,"Center DEC");
    fits_write_key(fitsimage,TSTRING,"DEC",wcstemp,comment,&status);
    strcpy(comment,"Equinox of coordinates");
    fits_write_key(fitsimage,TDOUBLE,"EQUINOX",&((*wcsim).equinox),comment,&status);
/*     //printf(" Cerca if \n"); */
    ok1=ffgcrd(fitsimage,"WEPOCH",snul,&status);status=0;
    strcpy(comment,"Epoch of coordinates");
    fits_write_key(fitsimage,TDOUBLE,"EPOCH",&((*wcsim).epoch),comment,&status);
    strcpy(comment,"Coordinate system");
    fits_delete_key(fitsimage, "RADECSYS", &status); status=0;
    if(fits_write_key(fitsimage,TSTRING,"RADECSYS",(*wcsim).radecsys,comment,&status)) fits_report_error(stderr,status);
/*     //printf(" status %d",status); */

    
    
    
    /* Set standard FITS WCS keywords */
/*     //printf(" La proyeccion es %s = %d\n",(*wcsim).ptype,(*wcsim).prjcode); */
    strcpy (wcstemp, "RA---");
    strcat (wcstemp, wcsproj);
    strcpy(comment,"Coordinate type in X");
    fits_delete_key(fitsimage, "CTYPE1", &status); status=0;
    if(fits_write_key(fitsimage,TSTRING,"CTYPE1",wcstemp,comment,&status))  fits_report_error(stderr,status) ;
/*     //hputs  (header, "CTYPE1", wcstemp); */
    strcpy (wcstemp, "DEC--");
    strcat (wcstemp, wcsproj);
/*     //hputs  (header, "CTYPE2", wcstemp); */
    strcpy(comment,"Coordinate type in Y");
    fits_delete_key(fitsimage, "CTYPE2", &status); status=0;
    if(fits_write_key(fitsimage,TSTRING,"CTYPE2",wcstemp,comment,&status))  fits_report_error(stderr,status);
    strcpy(comment,"crval 1");
    fits_delete_key(fitsimage, "CRVAL1", &status); status=0;
    if(fits_write_key(fitsimage,TDOUBLE,"CRVAL1",&((*wcsim).xref),comment,&status))  fits_report_error(stderr,status);
    fits_report_error(stderr,status) ;

    strcpy(comment,"crval 2");
    fits_delete_key(fitsimage, "CRVAL2", &status); status=0;
    if(fits_write_key(fitsimage,TDOUBLE,"CRVAL2",&((*wcsim).yref),comment,&status))  fits_report_error(stderr,status);
    fits_report_error(stderr,status) ;
    
    strcpy(comment,"crpix 1");
    fits_delete_key(fitsimage, "CRPIX1", &status); status=0;
    if(fits_write_key(fitsimage,TDOUBLE,"CRPIX1",&((*wcsim).xrefpix),comment,&status))  fits_report_error(stderr,status);
    strcpy(comment,"crpix 2");
    fits_delete_key(fitsimage, "CRPIX2", &status); status=0;
    if(fits_write_key(fitsimage,TDOUBLE,"CRPIX2",&((*wcsim).yrefpix),comment,&status))  fits_report_error(stderr,status);
/*     //hputnr8 (header, "CRVAL1", 9, (*wcsim).xref); */
/*     //hputnr8 (header, "CRVAL2", 9, (*wcsim).yref); */
/*     //hputnr8 (header, "CRPIX1", 4, (*wcsim).xrefpix); */
/*     //hputnr8 (header, "CRPIX2", 4, (*wcsim).yrefpix); */
    if ((*wcsim).rotmat) {
      printf(" Using rotation matrix\n");
      strcpy(comment,"Rotation matrix");
      fits_delete_key(fitsimage, "CD1_1", &status); status=0;
      if(fits_write_key(fitsimage,TDOUBLE,"CD1_1",&((*wcsim).cd[0]),comment,&status))  fits_report_error(stderr,status);
      strcpy(comment,"Rotation matrix");
      fits_delete_key(fitsimage, "CD1_2", &status); status=0;
      if(fits_write_key(fitsimage,TDOUBLE,"CD1_2",&((*wcsim).cd[1]),comment,&status))  fits_report_error(stderr,status);
      strcpy(comment,"Rotation matrix");
      fits_delete_key(fitsimage, "CD2_1", &status); status=0;
      if(fits_write_key(fitsimage,TDOUBLE,"CD2_1",&((*wcsim).cd[2]),comment,&status))  fits_report_error(stderr,status);
      strcpy(comment,"Rotation matrix");
      fits_delete_key(fitsimage, "CD2_2", &status); status=0;
      if(fits_write_key(fitsimage,TDOUBLE,"CD2_2",&((*wcsim).cd[3]),comment,&status))  fits_report_error(stderr,status);
/*       //hputnr8 (header, "CD1_1", 9, (*wcsim).cd[0]); */
/*       //hputnr8 (header, "CD1_2", 9, (*wcsim).cd[1]); */
/*       //hputnr8 (header, "CD2_1", 9, (*wcims).cd[2]); */
/*       //hputnr8 (header, "CD2_2", 9, (*wcsim).cd[3]); */
      fits_delete_key(fitsimage,"CDELT1",&status);status=0;
      fits_delete_key(fitsimage,"CDELT2",&status);status=0;
      fits_delete_key(fitsimage,"CROTA1",&status);status=0;
      fits_delete_key(fitsimage,"CROTA2",&status);status=0;
/*       hdel (header, "CDELT1"); */
/*       hdel (header, "CDELT2"); */
/*       hdel (header, "CROTA1"); */
/*       hdel (header, "CROTA2"); */
    }
    else {
      printf(" Using old rotation convention\n");
      strcpy(comment,"cdelt 1");
      fits_delete_key(fitsimage, "CDELT1", &status); status=0;
      if(fits_write_key(fitsimage,TDOUBLE,"CDELT1",&((*wcsim).xinc),comment,&status))  fits_report_error(stderr,status);
      strcpy(comment,"cdelt 2 ");
      fits_delete_key(fitsimage, "CDELT2", &status); status=0;
      if(fits_write_key(fitsimage,TDOUBLE,"CDELT2",&((*wcsim).yinc),comment,&status))  fits_report_error(stderr,status);
      strcpy(comment,"crot 1 ");
      fits_delete_key(fitsimage, "CROTA1", &status); status=0;
      if(fits_write_key(fitsimage,TDOUBLE,"CROTA1",&((*wcsim).rot),comment,&status))  fits_report_error(stderr,status);
      strcpy(comment,"crot 2");
      fits_delete_key(fitsimage, "CROTA2", &status); status=0;
      if(fits_write_key(fitsimage,TDOUBLE,"CROTA2",&((*wcsim).rot),comment,&status))  fits_report_error(stderr,status);
      /*       hputnr8 (header, "CDELT1", 9, (*wcsim).xinc); */
/*       hputnr8 (header, "CDELT2", 9, (*wcsim).yinc); */
/*       hputnr8 (header, "CROTA1", 3, (*wcsim).rot); */
/*       hputnr8 (header, "CROTA2", 3, (*wcsim).rot); */
      fits_delete_key(fitsimage,"CD1_1",&status);status=0;
      fits_delete_key(fitsimage,"CD1_2",&status);status=0;
      fits_delete_key(fitsimage,"CD2_1",&status);status=0;
      fits_delete_key(fitsimage,"CD2_2",&status);status=0;
/*       hdel (header, "CD1_1"); */
/*       hdel (header, "CD1_2"); */
/*       hdel (header, "CD2_1"); */
/*       hdel (header, "CD2_2"); */
    }
    fits_report_error(stderr,status);
    
    /* Set plate fit, if present */
    if ((*wcsim).ncoeff1 > 0) {
      char keyword[16];
      int i;
/*       //printf(" escribo polinm,o\n"); */
/*       //printf(" Status %d\n",status); */
      for (i = 0; i < (*wcsim).ncoeff1; i++) {
	sprintf (keyword, "CO1_%d",i+1);
	strcpy(comment,"Astrometric coefficient");
/* 	//printf(" pol %d  %f  key %s\n",i,(*wcsim).x_coeff[i],keyword); */
	fits_delete_key(fitsimage, keyword, &status); status=0;
	if(fits_write_key(fitsimage,TDOUBLE,keyword,&((*wcsim).x_coeff[i]),comment,&status))  fits_report_error(stderr,status);
/* 	//printf(" sta %d\n",status); */
/* 	//hputr8 (header, keyword, (*wcsim).x_coeff[i]); */
      }
    }
    if ((*wcsim).ncoeff2 > 0) {
      char keyword[16];
      int i;
/*       //printf(" escribo polinm,o 222\n"); */
      for (i = 0; i < (*wcsim).ncoeff2; i++) {
	sprintf (keyword, "CO2_%d",i+1);
	strcpy(comment,"Astrometric coefficient");
	fits_delete_key(fitsimage, keyword, &status); status=0;
	if(fits_write_key(fitsimage,TDOUBLE,keyword,&((*wcsim).y_coeff[i]),comment,&status))  fits_report_error(stderr,status);
/* 	//hputr8 (header, keyword, (*wcsim).y_coeff[i]); */
      }
    }
/*     //    printf(" SALGo\n"); */
    return;
}



void DeleteClickIm(double *alfa  ,double *delta ,
		   double *xpix,  double *ypix)
{

/*   char option; */
  float xcur,ycur;
  int idel=-1,j;
  char cnul;
/*   float armin,armax,decmin,decmax; */
  float mindist=1.e30,dist;
  
  
  printf(" Click on image...\n");
  cpgsvp(0.1,0.5,0.1,0.9);
  cpgwnad(0.,naxes[0],0.,naxes[1]);
  cpgcurs(&xcur,&ycur,&cnul);
  for(j=0;j<ncross;j++) {
    dist=((xpix[j]-xcur)*(xpix[j]-xcur)+(ypix[j]-ycur)*(ypix[j]-ycur));
    printf(" Is %d, di %f min %f\n",idel,dist,mindist);
    if(mindist>dist) {
      mindist=dist;
      idel=j;
    }
  }
  printf(" Cursor position %f %f\n",xcur,ycur);
  printf(" Nearest object position %f %f\n",xpix[idel],ypix[idel]);
  cpgsci(5);
  cpgsch(3.);
  cpgpt1((float)xpix[idel],(float)ypix[idel],3);
  cpgsci(1);
  cpgsch(1.);
  

  printf(" Deleting object number %d out of %d\n",idel,ncross);


  memmove(xpix+idel,xpix+idel+1,((ncross)-idel-1)*sizeof(double));
  memmove(ypix+idel,ypix+idel+1,((ncross)-idel-1)*sizeof(double));
  memmove(alfa+idel,alfa+idel+1,((ncross)-idel-1)*sizeof(double));
  memmove(delta+idel,delta+idel+1,((ncross)-idel-1)*sizeof(double));
  ncross--;
}

void DeleteClickErr(double *alfa  ,double *delta ,
		   double *xpix,  double *ypix  , 
                   float *xerror, float *yerror)
{

/*   char option; */
  float xcur,ycur;
  int idel=-1,j;
  char cnul;
  float xmin,xmax,ymin,ymax;
  float mindist=1.e30,dist;
    
    
  printf(" Click on error plot...\n");
  pgLimits(ncross,xerror,&xmin,&xmax);
  pgLimits(ncross,yerror,&ymin,&ymax);
  cpgsvp(0.6,0.9,0.1,0.5);
  cpgwnad(xmin,xmax,ymin,ymax);

  cpgcurs(&xcur,&ycur,&cnul);
  for(j=0;j<ncross;j++) {
    dist=((xerror[j]-xcur)*(xerror[j]-xcur)+(yerror[j]-ycur)*(yerror[j]-ycur));
    printf(" Is %d, di %f min %f\n",idel,dist,mindist);
    if(mindist>dist) {
      mindist=dist;
      idel=j;
    }
  }
  printf(" Cursor position %f %f\n",xcur,ycur);
  printf(" Nearest object position %f %f\n",xerror[idel],yerror[idel]);
  cpgsci(5);
  cpgsch(3.);
/*   //  cpgpt1((float)xerror[idel],(float)yerror[idel],3); */
  cpgsci(1);
  cpgsch(1.);
  

  printf(" Deleting object number %d out of %d\n",idel,ncross);


  memmove(xpix+idel,xpix+idel+1,((ncross)-idel-1)*sizeof(double));
  memmove(ypix+idel,ypix+idel+1,((ncross)-idel-1)*sizeof(double));
  memmove(alfa+idel,alfa+idel+1,((ncross)-idel-1)*sizeof(double));
  memmove(delta+idel,delta+idel+1,((ncross)-idel-1)*sizeof(double));
  (ncross)--;
}

void InputInitGuess( int *initguess)
{
/* //  double epoch; */
  char flipping;
  printf(" Input aproximate RA central coordinate in hours.\n");
  alfac=readf(alfac/15.);
  alfac=alfac*15;
  printf(" Input aproximate DEC central coordinate in degrees.\n");
  deltac=readf(deltac);
  printf(" Input aproximate plate scale in arcsec/pix: ");
  xpixsiz=readf(xpixsiz);
  ypixsiz=xpixsiz;
  printf(" Input observation epoch (ej: 1998.5): ");
  epoch=(double)readf(epoch);
  printf(" Rotation of image with respect to sky: ");
  rotang=readf(rotang);
  printf(" Flipping?: ");
  flipping=readc('n');
  if(flipping=='y')   flip=1;
  else flip=0;
  wcsim=wcsxinit((double)(alfac),(double)(deltac),(double)(xpixsiz),naxes[0]/2.,naxes[1]/2.,naxes[0],naxes[1],0.,2000,epoch,"TAN");
  if(!flip) wcsdeltset(wcsim, wcsim->cdelt[0],wcsim->cdelt[1],(double)rotang);
  else      wcsdeltset(wcsim,-wcsim->cdelt[0],wcsim->cdelt[1],(double)rotang);
  if(iswcs(wcsim)) printf(" Initial WCS solution has been set\n");

  *initguess=1;
}


void InitWCS(char *inputfile,int *initguess)
{
  int lhead,nbfits;
  char snul[1000];
  double dnul1,dnul2;
  char really;
  char flipping;

  if ((header = fitsrhead (inputfile, &lhead, &nbfits)) == NULL) {
    fprintf (stderr, "Cannot read FITS header of file %s\n", inputfile);
    exit(1);
  }
/*   printf(" CAbaezero <<%s>>\n",header); */
/*   //PrintWCS(header,1); */
/*   //printf("caca\n"); */

/*   //  exit(1); */
  
  
  if((wcsim=wcsinit(header))==NULL) {
    printf(" No WCS solution found in header\n");  
/*     //hgetra (header, "RA", alfac);     Estos son los que hay que buscar. */
/*     //hgetdec (header, "DEC", deltac); */
/*     //hgetr8 (header, "SECPIX", xpixsiz); */
/*     //    hgets (header, "CAT-RA",alfac);     // Estos son para INT WFC */
/*     //printf(" alfa %f \n",alfac); */
    
    if(hgets (header, "RA", 80,snul)!=0) {
      alfac=str2ra(snul);
/*       printf(" s <<%s>> alfa %f \n",snul,alfac); */
      hgets (header, "DEC", 80,snul);
      deltac=str2dec(snul);
      /* printf(" s <<%s>> delta %f \n",snul,deltac); */
      printf(" But... Central coordinates found in FITS header: AR= %f DEC= %f\n",alfac,deltac);
      if(hgetr4(header, "SECPPIX",&xpixsiz) == 0) {
	if(interact) {
	  printf(" Input aproximate plate scale in arcsec/pix: ");
	  xpixsiz=readf(1.);
	}
      }	
      ypixsiz=xpixsiz;
/*       alfac=55.5255;deltac=0.02507;xpixsiz=0.37;ypixsiz=0.37;  */
/*       printf(" Epocha %f\n",epoch); */
/*       exit(1); */
      if(hgetr8(header, "EPOCH", &epoch) == 0) {
	if(interact) {
	  printf(" Input observation epoch (ej: 1998.5): ");
	  epoch=(double)readf(epoch);
	}
      }
      wcsim=wcsxinit((double)(alfac),(double)(deltac),(double)(xpixsiz),naxes[0]/2.,naxes[1]/2.,naxes[0],naxes[1],0.,2000,epoch,"TAN");
      if(interact) {
	printf(" Rotation of image with respect to sky: ");
	rotang=readf(rotang);
	printf(" Flipping?: ");
	flipping=readc('n');
	if(flipping=='y')   flip=1;
	else flip=0;
      }
      if(!flip) wcsdeltset(wcsim, wcsim->cdelt[0],wcsim->cdelt[1],(double)rotang);
      else      wcsdeltset(wcsim,-wcsim->cdelt[0],wcsim->cdelt[1],(double)rotang);
      printf(" Useful information found in FITS header:\n");
      printf(" Central coordinates: AR= %f DEC= %f\n",alfac,deltac);
      printf(" Plate scale: %f\n",xpixsiz);
      printf(" Epoch: %f\n",epoch);
      *initguess=1;
    }
    
  }
  else {
/*     wcseq(header,wcsim); */
    printf(" WCS information from FITS header:\n");
    PrintWCS(header,1);
    pix2wcs(wcsim,naxes[0]/2., naxes[1]/2. ,&dnul1, &dnul2);
    alfac=(float)dnul1;deltac=(float)dnul2;
    xpixsiz=fabs((float)((*wcsim).xinc)*3600.);
    ypixsiz=fabs((float)((*wcsim).yinc)*3600.);
/*     //pix2wcs(wcsim,naxes[0]/2., naxes[1]/2. ,&alfac, &deltac); */
    *initguess=1;
    pix2wcs(wcsim,naxes[0]/2., naxes[1]/2. ,&dnul1, &dnul2);
    printf(" Using central coordinates RA = %f Dec = %f\n",dnul1,dnul2);
    printf(" Central pixels %f %f %ld %ld\n",naxes[0]/2., naxes[1]/2.,naxes[0],naxes[1]);
    if(hgetr8(header, "EPOCH", &epoch) == 0)     (*wcsim).epoch=epoch;
    printf(" Using epoch %f\n",epoch);
/*     pix2wcs(wcsim,10000.,10000. ,&dnul1, &dnul2); */
/*     printf(" alfa %f delta %f\n",dnul1,dnul2); */

    if(interact) {
      
      printf(" There is already a WCS solution set.\n Continuing will delete this solution, but it will be used as a starting point.\n Do you really want to continue (y/n)?: ");
      /* if(!interact) exit(1); */
      really=readc('n');
      if(really=='y') {
      }
      else {
	printf(" There is already a WCS solution set. Exiting\n");
	exit(1);
      }
    }
  }
}






void LoadParam(char file[100])
{
  fitsfile *parfile;
/*   char keyf[9]="",key[9]=""; */
  int status=0;
/*   //char string[51]; */
  char comment[51];

  printf("Using parameter file <<%s>>\n",file);
  if( ffopen2(&parfile,file, READONLY, &status)) fits_report_error(stderr,status
);
  
  
  ffgky(parfile,TSTRING,"IMAGE",inputfile,comment,&status);
  fits_report_error(stderr,status);
/*   //  printf(" IMAGE\n"); */
  ffgky(parfile,TSTRING,"ASTCAT",astcat,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"SECCAT",seccat,comment,&status);
  fits_report_error(stderr,status);
/*   //  printf(" ASTCAT\n"); */
  ffgky(parfile,TINT,"POLDEG",&ngrad,comment,&status);
  fits_report_error(stderr,status);
/*   //  printf(" POLDEG\n"); */
  ffgky(parfile,TSTRING,"OBJCAT",bfile,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"XCOL",&bxcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"YCOL",&bycol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"FLUXCOL",&bfcol,comment,&status);
  fits_report_error(stderr,status);
/*   //  printf(" OBJCAT\n"); */
  ffgky(parfile,TINT,"NOBJ",&nbrightestobj,comment,&status);
  fits_report_error(stderr,status);
/*   //  printf(" NOBJ\n"); */
  ffgky(parfile,TINT,"NSTAR",&nbrighteststar,comment,&status);
  fits_report_error(stderr,status);
/*   //  printf(" NSTAR\n"); */
  ffgky(parfile,TSTRING,"DEVICE",pgdevice,comment,&status);
  fits_report_error(stderr,status);
/*   // printf(" ROTANG\n"); */
  ffgky(parfile,TFLOAT,"ROTANG",&rotang,comment,&status);
  fits_report_error(stderr,status);
/*   // printf(" FLIP\n"); */
  ffgky(parfile,TLOGICAL,"FLIP",&flip,comment,&status);
  fits_report_error(stderr,status);
/*   //  exit(1); */
/*   //  printf(" DEVICE\n"); */
  ffgky(parfile,TFLOAT,"RA",&alfac,comment,&status);
  fits_report_error(stderr,status);
  alfac=alfac*15;
/*   //  printf(" RA\n"); */
  ffgky(parfile,TFLOAT,"DEC",&deltac,comment,&status);
  fits_report_error(stderr,status);
/*   //  printf(" DEC\n"); */
/*     printf(" ecpovh %f antes\n",epoch); */

  ffgky(parfile,TDOUBLE,"EPOCH",&epoch,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TFLOAT,"SECPIX",&xpixsiz,comment,&status);
  fits_report_error(stderr,status);
/*   //  printf(" SECPIX\n"); */
  ypixsiz=xpixsiz;

  if(status!=0) {
    fits_report_error(stderr,status);
    printf("ERROR: Reading parameter file. Not enough keywords.\n");
    exit(1);
  }

  printf(" Parameter file read\n");
/*   printf(" ecpovh %f ra %f dec %f sep %f rot %f pol %d\n",epoch,alfac,deltac,xpixsiz,rotang,ngrad); */

}


int min(int x1,int x2) {
  if(x1>x2) return(x2);
  else return(x1);
}


/* Que solucion astrometrica esta usando?
   El tipo de proyeccion es 2: WCS_TAN
   El "pix2wcs" llama a "wcspos" para ese tipo de proyeccion.
   El "wcspos" llama al "wcsrev"
   El "wcsrev" llama al "celrev"
   El "celrev" hace una cosa rara. Llama a la subrutina
   cel->prjrev(), que es un puntero a una subrutina.
   Para el caso TAN, la subrutina se llama "tanrev"


   if (linrev(pixcrd, lin, imgcrd)) {
   world[j] = imgcrd[j] + crval[j];
   

   r = sqrt(x*x + y*y);
   if (r == 0.0) {
   *phi = 0.0;
   } else {
   *phi = atan2deg(x, -y);
   }
   *theta = atan2deg(prj->r0, r);
   


   wcscon (wcs->syswcs,wcs->sysout,eqin,eqout,&xp,&yp,wcs->epoch);
   *xpos = xp;
   *ypos = yp;
   
*/


void ShiftSolution(double *alfas, double *deltas,
                   double *xpixs, double *ypixs)
{

  double mura,mudec,mod;   /* Dirección a la que apunta la corrección */
/*   double xpixscal[NSECMAX],ypixscal[NSECMAX];  */  /* De alfa,delta a xpixcal,ypixcal */
  double alfascal[NSECMAX],deltascal[NSECMAX];  /* De xpix,ypix a alfacal, deltacal */
  double alfashift[NSECMAX],deltashift[NSECMAX];  /* De xpix,ypix a alfacal, deltacal */
  double dnul1,dnul2;
  double oldalfac,olddeltac;
  double newalfac,newdeltac;
  int i;
  char *radecsys;
  char option;

  static float factorarr=3000;

  cpgsclp(1);
/*   printf(" NUEo.\n"); */
  
  for(i=0;i< ncrosssec; i++ ) {
    pix2wcs(wcsim,xpixs[i], ypixs[i] , &dnul1, &dnul2);
    alfascal[i]=dnul1;deltascal[i]=dnul2;
    alfashift[i] = alfas[i]-alfascal[i];
    deltashift[i]=deltas[i]-deltascal[i];
    if(rotang==0 || rotang==180) 
      cpgarro((float)(alfas[i]),(float)(deltas[i]),(float)(alfas[i]-(alfashift[i])*factorarr),(float)(deltas[i]-(deltashift[i])*factorarr));
    else 
      cpgarro((float)(deltas[i]),(float)(alfas[i]),(float)(deltas[i]-(deltashift[i])*factorarr),(float)(alfas[i]-(alfashift[i])*factorarr));
    printf(" RA %f DEC %f RACAL %f DECCAL %f  XPIX %f YPIX %f\n",alfas[i],deltas[i],alfascal[i],deltascal[i],xpixs[i],ypixs[i]);
    printf("      Shift    %f %f\n",alfashift[i]*3600., deltashift[i]*3600.);
  }

  wcsfull(wcsim,&oldalfac, &olddeltac,&dnul1,&dnul2);
  mura =+16.88*sin(olddeltac/180.*PI)*sin(oldalfac/180.*PI)+9.75*cos(olddeltac/180.*PI);
  mudec=+16.88*cos(olddeltac/180.*PI);
  mod=(sqrt(mura*mura+mudec*mudec));
  mura=mura/mod;
  mudec=mudec/mod;

  printf(" Correction must be in the direction %f %f\n",mura,mudec);
  cpgsci(3);
  if(rotang==0 || rotang==180) 
    cpgarro((float)oldalfac,(float)olddeltac,(float)oldalfac+mura*factorarr,(float)olddeltac+mudec*factorarr);
  else 
    cpgarro((float)olddeltac,(float)oldalfac,(float)olddeltac+mudec*factorarr,(float)oldalfac+mura*factorarr);
  cpgsci(1);

  if((radecsys=malloc(32*sizeof(char)))==NULL ) printf("I cannot dimension radecsys\n");
  radecsys=getradecsys(wcsim);
  printf(" In principle, central coordinates:  RA = %f Dec = %f\n",oldalfac,olddeltac);
  /* De la otra manera cambio las coordenadas centrales, pero wcsshift 
 lo que cambia realmente es wcs->xref, wcs->yref, que no tienen por que ser 
 las coordenadas reales del centro */
  oldalfac =(*wcsim).xref;
  olddeltac=(*wcsim).yref;
  
  printf(" Shifting central coordinates from RA = %f Dec = %f\n",oldalfac,olddeltac);

  shiftra=StMedia_d(ncrosssec,alfashift,&dnul1);
  shiftdec=StMedia_d(ncrosssec,deltashift,&dnul1);

  cpgsci(2);
  if(rotang==0 || rotang==180) 
    cpgarro((float)oldalfac,(float)olddeltac,(float)oldalfac-shiftra*factorarr,(float)olddeltac-shiftdec*factorarr);
  else 
    cpgarro((float)olddeltac,(float)oldalfac,(float)olddeltac-shiftdec*factorarr,(float)oldalfac-shiftra*factorarr);
  cpgsci(1);

  
  printf(" Shifting %f %f\n",shiftra*3600.,shiftdec*3600.);

  newalfac=oldalfac+shiftra;
  newdeltac=olddeltac+shiftdec;
  printf(" to RA = %f Dec = %f\n",newalfac,newdeltac);

  wcsshift(wcsim,newalfac,newdeltac,radecsys);
  
  oldalfac =(*wcsim).xref;
  olddeltac=(*wcsim).yref;
  
  printf(" Comprobando RA = %f Dec = %f\n",oldalfac,olddeltac);

  wcsfull(wcsim,&oldalfac, &olddeltac,&dnul1,&dnul2);
  printf(" At the end, central coordinates:  RA = %f Dec = %f dentro de radecsys %s\n",oldalfac,olddeltac,radecsys);
  
  if(interact) {
    printf(" Make the shift [y/n]?: ");
    option=readc('y'); 
    if(option=='y');
    else {
      wcsshift(wcsim,oldalfac,olddeltac,radecsys);
      wcsfull(wcsim,&oldalfac, &olddeltac,&dnul1,&dnul2);
      printf(" Going back to RA = %f Dec = %f\n",oldalfac,olddeltac);
    }   
  }

}



void SaveParam(char parfilename[100])
{

/*   int i,j; */
  int nc=0,nt;
  char ch51[51],ch21[21];
  float secpix;
  FILE *parfile;
  secpix=xpixsiz;
/*   printf("Name of file (%s) :",parfilename); */
/*   getline(parfilename,100); */
  if((parfile=fopen(parfilename,"w")) ==NULL) {
    printf("ERROR: Can't open file\n");
    return;
  }
  fprintf(parfile,"COMMENT  Parameter file for PlateAstrom                                        \n");
  nc++;

  sprintf(ch51,"'%s'",inputfile);
  fprintf(parfile,"IMAGE   = %-51.51s / Input FITS file\n",ch51);
  sprintf(ch21,"'%s'",astcat   );
  fprintf(parfile,"ASTCAT  = %-21s/ Primary star cat (GSC,APM,USNO,ACT,Tycho-2)   \n",ch21);
  if(secflag)  sprintf(ch21,"'%s'",seccat   );
  else  sprintf(ch21,"'%s'","NONE");
  fprintf(parfile,"SECCAT  = %-21s/ Secondary star cat (ACT,Tycho-2)              \n",ch21);
  sprintf(ch51,"'%s'",bfile);
  fprintf(parfile,"OBJCAT  = %-51.51s / Object file    \n",ch51     );
  fprintf(parfile,"XCOL    =%21d / Column in OBJCAT with object X position       \n",bxcol                           );
  fprintf(parfile,"YCOL    =%21d / Column in OBJCAT with object Y position       \n",bycol                           );
  fprintf(parfile,"FLUXCOL =%21d / Column in OBJCAT with object flux position    \n",bfcol                           );
  fprintf(parfile,"NOBJ    =%21d / Number of objects to match                    \n",nbrightestobj                   );
  fprintf(parfile,"NSTAR   =%21d / Number of objects to match                    \n",nbrighteststar                  );
  fprintf(parfile,"RA      =%21f / Approximate RA of the center of the image     \n",alfac/15.                       );
  fprintf(parfile,"DEC     =%21f / Approximate DEC of the center of the image    \n",deltac                          );
  fprintf(parfile,"EPOCH   =%21f / Epoch of the observation (Date of the iamge)  \n",epoch                          );
  fprintf(parfile,"SECPIX  =%21f / Approximate plate scale                       \n",secpix                          );
  fprintf(parfile,"ROTANG  =%21f / Rotation of the image in RA-DEC plane         \n",rotang                          );
  fprintf(parfile,"FLIP    =%21d / Is image flipped? (0=NO, 1=YES)               \n",flip                            );
  fprintf(parfile,"POLDEG  =%21d / Polinomial degree for plate solution          \n",ngrad                           );
  sprintf(ch51,"'%s'","/xserve");
  fprintf(parfile,"DEVICE  = %-51.51s / PGPLOT Device  \n",ch51        );
  nc+=17;
  
  fprintf(parfile,"COMMENT  Blank lines to assure the last block of 2880 caracters is complete    \n");
  fprintf(parfile,"COMMENT                                                                        \n");
  nc+=2;
  
  for(nt=nc;nt<(1+(int)(nc/36))*36-1;nt++) {
/*     //printf("NC %d %d\n",nt,(1+(int)(nc/36.))*36); */
    fprintf(parfile,"COMMENT                                                                        \n");
  }
  fprintf(parfile,"END                                                                            \n");
  fclose(parfile);
}

int incompare(const void *X1,const void *X2) 
{  
  struct obj_table *x1,*x2;
  x1=(struct obj_table*) X1;
  x2=(struct obj_table*) X2;
  if((*x1).flux  < (*x2).flux) return(+1);
  if((*x1).flux  > (*x2).flux) return(-1);
  if((*x1).flux == (*x2).flux) return( 0);
  return(0);
}
int iscompare(const void *X1,const void *X2)
{  
  struct star_table *x1,*x2;
  x1=(struct star_table*)X1;
  x2=(struct star_table*)X2;
  if((*x1).mag  < (*x2).mag) return(-1);
  if((*x1).mag  > (*x2).mag) return(+1);
  if((*x1).mag == (*x2).mag) return( 0);
  return(0);
}

void StarMatch_ceg(int *ncross, double *alfa, double *delta, double *xpix, double *ypix, int *off, int ns,double  *sx, double *sy, int ng, double *gra, double *gdec, double *gx, double *gy, double tol, struct WorldCoor *wcs, int debug)
{
  int nmatch;
  int i,j;
  float dxys,dxy;
  int igs;
  int nrealmatch;
  for(i=0;i<ng;i++) {
    gx[i]=0;gy[i]=0;
    wcs2pix(wcsim,gra[i],gdec[i],&gx[i],&gy[i],&off[i]);
  }

  nmatch=StarMatch(ns,sx,sy,ng,gra,gdec,gx,gy,tol,wcsim,debug);

  printf(" De Star;Atc %d\n",nmatch);

  for(i=0;i<ng;i++) {
    gx[i]=0;gy[i]=0;
    wcs2pix(wcsim,gra[i],gdec[i],&gx[i],&gy[i],&off[i]);
  }
 

  nrealmatch=0;
  for (i = 0; i < ns; i++) {
    dxys = -1.0;
    igs = -1;
    for (j = 0; j < ng; j++) {
      if (!off[j]) {
	dxy=(gx[j]-sx[i])*(gx[j]-sx[i])+(gy[j]-sy[i])*(gy[j]-sy[i]);
/* 	printf(" i %d j %d dxy %f tol %f dxys %f igs %d\n",i,j,dxy,tol,dxys,igs); */
	if (dxy < tol*tol) {
	  if (dxys < 0.0) {
	    dxys = dxy;
	    igs = j;
	  }
	  else if (dxy < dxys) {
	    dxys = dxy;
	    igs = j;
	  }
	}
      }
    }
    if (igs > -1) {
      alfa[nrealmatch]=gra[igs];
      delta[nrealmatch]=gdec[igs];
      xpix[nrealmatch]=sx[i];
      ypix[nrealmatch]=sy[i];
      nrealmatch++;
    }
  }

  *ncross=nrealmatch;
}



void MatchTol(int *ncross, double *alfa, double *delta, double *xpix, double *ypix, int *off, int ns,double  *sx, double *sy, int ng, double *gra, double *gdec, double *gx, double *gy, double tol, struct WorldCoor *wcs, int debug)
{
  int i,j;
  float dxys,dxy;
  int igs;
  int nrealmatch;
  for(i=0;i<ng;i++) {
    gx[i]=0;gy[i]=0;
    wcs2pix(wcsim,gra[i],gdec[i],&gx[i],&gy[i],&off[i]);
  }

  for(i=0;i<ng;i++) {
    gx[i]=0;gy[i]=0;
    wcs2pix(wcsim,gra[i],gdec[i],&gx[i],&gy[i],&off[i]);
  }
 

  nrealmatch=0;
  for (i = 0; i < ns; i++) {
    dxys = -1.0;
    igs = -1;
    for (j = 0; j < ng; j++) {
      if (!off[j]) {
	dxy=(gx[j]-sx[i])*(gx[j]-sx[i])+(gy[j]-sy[i])*(gy[j]-sy[i]);
/* 	printf(" i %d j %d dxy %f tol %f dxys %f igs %d\n",i,j,dxy,tol,dxys,igs); */
	if (dxy < tol*tol) {
	  if (dxys < 0.0) {
	    dxys = dxy;
	    igs = j;
	  }
	  else if (dxy < dxys) {
	    dxys = dxy;
	    igs = j;
	  }
	}
      }
    }
    if (igs > -1) {
      alfa[nrealmatch]=gra[igs];
      delta[nrealmatch]=gdec[igs];
      xpix[nrealmatch]=sx[i];
      ypix[nrealmatch]=sy[i];
      nrealmatch++;
    }
  }

  *ncross=nrealmatch;
}
