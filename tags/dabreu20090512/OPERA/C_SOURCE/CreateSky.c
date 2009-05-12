/*  */
/*   Obtiene una imagen del fondo de cielo y de la desviacion del fondo */
/*   ajustando una gaussiana a la moda de regiones XSCELLxYSCELL de */
/*   la imagen de entrada */
/*   La imagen de salida tendra, por tanto, NAXIS1/XSCELL x NAXIS2/YSCELL */
/*   pixels */
/*   Las imagenes de salida, de fondo y sigma, se escriben en formato FITS */
/*   BITPIX=16, y se actualizan los descriptores CRVAL y CDELT de modo */
/*   que las coordenadas de los pixels de estas imagenes correspondan */
/*   a los centros de las cajas de la imagen original */
/*  */
/*   ccpgp CreateSky '-lFITS -lUtil' */
/*   ccl   CreateSky '-lFITS -lUtil' */
/*  ----------------------------------------------------------------------- */

/* #include <stdlib.h> */
/* #include <stdio.h> */
/* #include <malloc.h> */
/* #include <math.h> */
#include "modulos.h"

#define MAXITER 10000
#define DEBUG 0

static int icompare(x1,x2)
     float *x1,*x2; 
{  
/*   //printf(" x1 %f x2 %f\n",*x1,*x2); */
  if(*x1  < *x2) return(-1);
  if(*x1 == *x2) return(0);
  if(*x1  > *x2) return(+1);

  return(0);
}

float Sky(int n, float *ima, int nbin, float *min, float *max, float *sig,
          float ftol, float *p);
float Amoe_Funk(int n, float *x, float *y, float *p);
void GraphHisto(int n, float *caja, int nbin, float min, float max,
                float moda, float sigma, float *p);
void Decimil(int n,float *x,float moda, float *decmin, float *decmax);

int main(int argc, char **argv)
{

  FILE *fp;
  char fitsin[19],fitssky[19],fitssig[19],fitserrsky[19],dev[19];
  int xscell,yscell,xmin,xmax,ymin,ymax,nbin;
/*   float xbin,ybin; */
  int dimsx,dimsy;
  float imamin,imamax,ftol, fondo,sigma;
  float imin,imax;
  float *imasky,*imasig,*imaerrsky,*caja,*oldimasky,*oldcaja;
  struct headfits hin,hout;
  int x1,x2,y1,y2,c=0,i,j,pgask;
  long blc[2],trc[2],incr[2];
  int pgbegx,pgbegy;
  int nact,ntot;
  char chc[71],chl[80],graphmod;
  float deltpixx,deltpixy;
  float tr[6]={0.,1.,0.,0.,0.,1.},p[5];
  float first,third;
  float median;

  
  /* Variables para el FITS  */
  fitsfile *infits,*skyfits,*sigfits,*errskyfits;
  int status=0;
  int nfound, anynull;

  long naxes[2],  npixels;
/*   float datamin, datamax, nullval; */



/* // ------------------- */
/* // Inicio del programa */
/* // ------------------- */
/*   //system("clear"); */
/*   //printf(" This is CreateSky \n\n"); */
/*   //system("banner CreateSky"); */
/*   //  SlpTime(stdout); */
/*   //printf("\n              @OAL(1996)\n\n"); */
 
  if(argc != 2) {
    printf("Use: CreateSky parameter_file\n\n");
    exit(1);
    }

/* // ------------------------- */
/* // Lectura de los parametros */
/* // ------------------------- */
  if((fp=fopen(argv[1],"r"))==NULL) {
    fprintf(stderr,"CreateSky: ERROR. I cannot open parameter file %s\n",argv[1]);
    exit(1);
  }
  c += f_kfs(fp,"IMAGE   =",fitsin);
  c += f_kfs(fp,"SKY     =",fitssky);
  c += f_kfs(fp,"SIG     =",fitssig);
  c += f_kfs(fp,"ERRSKY  =",fitserrsky);
  c += f_kfi(fp,"XSCELL  =",&xscell);
  c += f_kfi(fp,"YSCELL  =",&yscell);
  c += f_kfi(fp,"XMIN    =",&xmin);
  c += f_kfi(fp,"XMAX    =",&xmax);
  c += f_kfi(fp,"YMIN    =",&ymin);
  c += f_kfi(fp,"YMAX    =",&ymax);
  c += f_kff(fp,"IMAMIN  =",&imamin);
  c += f_kff(fp,"IMAMAX  =",&imamax);
  c += f_kfi(fp,"NBIN    =",&nbin);
/*   c += f_kff(fp,"FTOL    =",&ftol); */
  ftol=5e-6;
  c += f_kfs(fp,"DEVICE  =",dev);
  c += f_kfi(fp,"PGASK   =",&pgask);
  c += f_kfc(fp,"GRAPHMOD=",&graphmod);
  c += f_kfi(fp,"PGBEGX  ",&pgbegx);
  c += f_kfi(fp,"PGBEGY  ",&pgbegy);
  fclose(fp);
  if(c != 18) {
    printf("Not enough parameters in parameter file %s (%d)\n",argv[1],c);
    exit(1);
    }

  if(nbin>xscell*yscell) {
    printf(" Cannot use more bins (%d) than pixels in cells (%d x %d)\n",nbin,xscell,yscell);
    exit(1);
  }

/* // ---------------------------------------- */
/* // Impresion de los parametros del programa */
/* // ---------------------------------------- */
  printf("\n\n");
  printf(" Parameters of the file  %s\n",argv[1]);
  printf("Input image ............................. %19s\n",fitsin);
  printf("Sky bakground output image............... %19s\n",fitssky);
  printf("Sky sigma output image................... %19s\n",fitssig);
  printf("Sky error output image................... %19s\n",fitserrsky);
  printf("Binning for computing sky................ %d x %d pixels\n",
		       xscell,yscell);
  printf("Image subsection......................... [%d:%d,%d:%d]\n",
		 xmin,xmax,ymin,ymax);
  printf("Maximum and minimum values for histogram. %f %f\n",imamin,imamax);
  printf("            No. de bins ................. %d\n",nbin);
  printf("---------------------------------------------------\n\n");

 
/* // -------------------- */
/* // Abro el modo grafico */
/* // -------------------- */
  if(graphmod == 'F') {     
    pgbegx=1;
    pgbegy=1;
    }
  if(graphmod == 'C' || graphmod == 'F') {
    cpgbeg(0,dev,pgbegx,pgbegy);
    cpgask(pgask);
    }

/* // ------------------------------------- */
/* // Leo la cabecera de la imagen original */
/* // ------------------------------------- */
  if(fitsrh(fitsin,&hin)) exit(1);  
  printf("Input imagen:\n");
  fitsinfo(fitsin,stdout);
  printf("\n--------------------------------------------\n");

  /* Ahora leo la informacion a traves de fitsio */
  if( ffopen(&infits, fitsin, READONLY, &status)) fits_report_error(stderr,status);
  if(fits_read_keys_lng(infits, "NAXIS", 1, 2, naxes, &nfound, &status)) fits_report_error(stderr,status);
  if(nfound!=2) {
    printf(" FITS file does not contain 2 dimensions\n");
    exit(1);
  }
  npixels=naxes[0]*naxes[1];




/* // --------------------------------- */
/* // Test de los valores de la subcaja */
/* // --------------------------------- */
  if(xmin == 0 && ymin == 0 && xmax == 0 && ymax == 0) {
    xmin=1;
    xmax=hin.naxis1;
    ymin=1;
    ymax=hin.naxis2;
    printf("Ajustando la subcaja al tamano maximo de la imagen\n");
    printf("Subcaja: [%d:%d;%d:%d]\n",xmin,xmax,ymin,ymax);
    }

  if(xmin < 1 || ymin < 1 || xmax > hin.naxis1 || ymax > hin.naxis2
     || xmax <= xmin  || ymax <= ymin) {
    fprintf(stderr,"CreateSky: ERROR. Subcaja mal definida\n");
    fprintf(stderr,"         Si quieres toda la imagen, opcion 0,0,0,0\n");
    fprintf(stderr,"         Dimension de la imagen: %d x %d\n",
	    hin.naxis1,hin.naxis2);
    fprintf(stderr,"         Subcaja: [%d:%d,%d:%d]\n",xmin,xmax,ymin,ymax);
    exit(1);
  }
  dimsx=xmax-xmin+1;
  dimsy=ymax-ymin+1;
  
  /*  ---------------------------------------------------- */
  /*  Calculo los descriptores para las imagenes de salida */
  /*  ---------------------------------------------------- */
  hout = hin;   /*    Copio todos los descriptores y ahora calculare los cambios */
      
  hout.bitpix=16;
  
  /*  No. de cajas == No. de pixels de la imagen de salida */
  /*  ---------------------------------------------------- */
  hout.naxis1=(dimsx-1)/xscell+1;
  hout.naxis2=(dimsy-1)/yscell+1;
  
  /*     xbin=(dimsx-xscell)/(hout.naxis1-1); */
  /*   ybin=dimsy/hout.naxis2; */
  /*   printf(" xbin %f ybin %f\n",xbin,ybin); */
  /*  Coordenadas del primer pixel */
  /*  ---------------------------- */
  hout.crval1=hin.crval1+(xmin-1+(xscell-1.0)/2.0)*hin.cdelt1;
  hout.crval2=hin.crval2+(ymin-1+(yscell-1.0)/2.0)*hin.cdelt2;
  
  /*  Espaciado de las cajas == Tamano del pixel != Tamano de la caja (Solapam.) */
  /*  -------------------------------------------------------------------------- */
  /*    Esto era antes: */
  hout.cdelt1=(float)(dimsx-xscell)/(hout.naxis1-1)*hin.cdelt1;
  hout.cdelt2=(float)(dimsy-yscell)/(hout.naxis2-1)*hin.cdelt2;
  /*     printf(" CDELT Antes %f %f\n",hout.cdelt1,hout.cdelt2); */
  /*    Haciendo un binning: */
  /*   hout.cdelt1=hin.cdelt1*xscell; */
  /*   hout.cdelt2=hin.cdelt2*yscell; */
  /*   printf(" CDELT Ahora %f %f\n",hout.cdelt1,hout.cdelt2); */
  
  
  
  /*  Descriptores COMMENT y HISTORY */
  /*  ------------------------------ */
  f_clrl(hout.comment);
  f_clrl(hout.history);


  /*  Impresion de reslutados */
  /*  ----------------------- */
  printf("\nOutput image carateristics               \n");
  printf("Dimensions: %d x %d pixels\n",hout.naxis1,hout.naxis2);
  printf("Coord. 1st. pixel  (%f,%f)\n",hout.crval1,hout.crval2);
  printf("Pixel size         %f x %f\n",hout.cdelt1,hout.cdelt2);
  printf("--------------------------------------\n");
  /*   exit(1); */
    
  /*  -------------------------------------- */
  /*  Dimensionado de las imagenes de salida */
  /*  -------------------------------------- */
  if((imasky=malloc(npixels*sizeof(float))) == NULL) {
    printf("CreateSky: ERROR. No puedo dimensionar imasky de %ld bytes\n",
	   npixels*sizeof(float));
    exit(1);
  }
  if((oldimasky=malloc(hout.naxis1*hout.naxis2*sizeof(float))) == NULL) {
    printf("CreateSky: ERROR. No puedo dimensionar imasky de %d bytes\n",
	   hout.naxis1*hout.naxis2*sizeof(float));
    exit(1);
  }
  if((imasig=malloc(hout.naxis1*hout.naxis2*sizeof(float))) == NULL) {
    printf("CreateSky: ERROR. No puedo dimensionar imasig de %d bytes\n",
	   hout.naxis1*hout.naxis2*sizeof(float));
    exit(1);
    }
  if((imaerrsky=malloc(hout.naxis1*hout.naxis2*sizeof(float))) == NULL) {
    printf("CreateSky: ERROR. No puedo dimensionar imaerrsky de %d bytes\n",
	   hout.naxis1*hout.naxis2*sizeof(float));
    exit(1);
  }
  if((oldcaja=malloc(xscell*yscell*sizeof(float))) == NULL) {
    printf("CreateSky: ERROR. No puedo dimensionar caja de %d bytes\n",
	   xscell*yscell*sizeof(float));
    exit(1);
  }
  if((caja=malloc(xscell*yscell*sizeof(float))) == NULL) {
    printf("CreateSky: ERROR. No puedo dimensionar caja de %d bytes\n",
	   xscell*yscell*sizeof(float));
    exit(1);
  }
  
  /*  Apertura de la imagen de entrada (No reviso pues he leido la cabecera) */
  /*  ---------------------------------------------------------------------- */
  
  
  /*  --------------- */
  /*  Bucle Principal */
  /*  --------------- */
  
  /*  Espaciado en pixels (real) */
  /*  -------------------------- */
  deltpixx=(float)(xmax-xmin+1-xscell)/(hout.naxis1-1);
  deltpixy=(float)(ymax-ymin+1-yscell)/(hout.naxis2-1);

  if(DEBUG) printf("Aqui nm a\n");
  ntot=hout.naxis1*hout.naxis2;
  for(j=0; j<hout.naxis2; j++)
    for(i=0; i<hout.naxis1; i++) {
      nact=i+j*hout.naxis1+1;
      
      /*  Calculo y lectura de la caja  */
      /*  ---------------------------- */
      x1=xmin+(int)(i*deltpixx+.5);
      x2=x1+xscell-1;
      y1=ymin+(int)(j*deltpixy+.5);
      y2=y1+yscell-1;
      /*       f_rfb(fp,&(hin),oldcaja,x1,x2,y1,y2); Esta es la forma antigua*/
      /* Leo con FITSIO */
      blc[0]=x1;blc[1]=y1;trc[0]=x2;trc[1]=y2;incr[0]=1;incr[1]=1;
      if (fits_read_subset_flt(infits,1,2,naxes,blc,trc,incr,0.,caja,&anynull,&status)) fits_report_error(stderr,status);
      
      
      /*  Calculo del fondo de cielo y sigma */
      /*  ---------------------------------- */
      /*       Pongo nbin=128 en el calculo del fondo de cielo */
      imin=imamin; imax=imamax;
      if(DEBUG) printf(" Computing Sky for %s[%d:%d,%d:%d]\n",fitsin,x1,x2,y1,y2);
      fondo=Sky(xscell*yscell,caja,nbin,&imin,&imax,&sigma,ftol,p);
      if(DEBUG) printf(" Mimino %f maximo %f \n",imin,imax); 
      
      
      imasky[i+j*hout.naxis1]=fondo;
      imasig[i+j*hout.naxis1]=sigma;
      imaerrsky[i+j*hout.naxis1]=sigma/sqrt(0.95*xscell*yscell); /* El 0.95 es porque estoy cogiendo el 95% de la distribucion gaussiana para calcular el histograma (linea 385) */
      printf("Sky for %s[%d:%d,%d:%d] = %f +/- %f (%d/%d:%4.1f)\n",
	     fitsin,x1,x2,y1,y2,fondo,sigma,nact,ntot,nact*100.0/ntot);
      /*  Dibujo */
      /*  ------ */
      if(graphmod == 'F') {
        printf(" Coordenadas mundo de la imagen de cielo %f %f\n",(i)*hout.cdelt1+hout.crval1,(j)*hout.cdelt2+hout.crval2);
        cpgsvpXY(0.10,0.6,0.10,0.90,0.0,xscell+1.0,0.0,yscell+1.0);
        cpgswin(0.0,xscell+1.0,0.0,yscell+1.0);
        cpgbox("bctns",0.0,0,"bctns",0.0,0);
        sprintf(chl,"%s[%d:%d,%d:%d]",fitsin,x1,x2,y1,y2);
        cpglab("X (pix)","Y (pix)",chl);
	/* 	        cpggray(caja,xscell,yscell,1,xscell,1,yscell,fondo+3*sigma,fondo-sigma,tr); */
	cpggray(caja,xscell,yscell,1,xscell,1,yscell,imax,imin,tr);
        cpgsvp(0.65,0.90,0.12,0.90);
        GraphHisto(xscell*yscell,caja,nbin,imin,imax,fondo,sigma,p);
	if(DEBUG) {
	  Quartil(xscell*yscell,caja,&first,&median,&third);
	  cpgsci(3);
	  cpgmove(first,0);
	  cpgdraw(first,1000.);
	  cpgmove(median,0);
	  cpgdraw(median,1000.);
	  cpgmove(third,0);
	  cpgdraw(third,1000.);
	  cpgsci(1);
		  }
        cpglab("Pixel Value","N",chl);
        cpgpage();

        }
      if(graphmod == 'C') {
        cpgenv(0.,1.,0.,1.,0,-2);
        GraphHisto(xscell*yscell,caja,128,imin,imax,fondo,sigma,p);
        sprintf(chl,"%s[%d:%d,%d:%d]",fitsin,x1,x2,y1,y2);
        cpglab("Pixel Value","N",chl);
      }
    }

  if(DEBUG) printf(" SALIOOOOOOOO\n");
  free(caja);
  /*  Salvamos las imagenes */
  /*  --------------------- */
  for(j=0; j<hout.naxis2; j++)
    for(i=0; i<hout.naxis1; i++) {
      /*      printf(" Valor %d %d: %f\n",i,j,imasig[i+j*hout.naxis1]); */
    }
  sprintf(chc,"Imagen de fondo de cielo de %s",fitsin);
  f_actlab(hout.history,chc);
  sprintf(chc,"Subcaja [%d:%d,%d:%d]",xmin,xmax,ymin,ymax);
  f_actlab(hout.history,chc);
  sprintf(chc,"Caja para el calculo del fondo de %dx%d pixels",xscell,yscell);
  f_actlab(hout.history,chc);
  sprintf(chc,"Imagen de desviacion del fondo en %s",fitssig);
  f_actlab(hout.history,chc);
  sprintf(chc,"Imagen de error de fondo de cielo en %s",fitserrsky);
  f_actlab(hout.history,chc);
  /*   fitswf(fitssky,&hout,imasky,'T'); */
  
  /*Grabo con FITSIO */
  if(DEBUG) printf(" ante sde lasd \n"); 
  status=0;
  ffinit(&skyfits,fitssky,&status);
  if(status) fits_report_error(stderr,status);
  naxes[0]=hout.naxis1;
  naxes[1]=hout.naxis2;
  fits_create_img(skyfits,-32,2,naxes,&status);
  fits_write_img(skyfits,TFLOAT,1,naxes[0]*naxes[1],imasky,&status);
  if(fits_write_key(skyfits,TSTRING,"BUNIT",hout.bunit,"pixel units",&status))  fits_report_error(stderr,status);
  if(fits_write_key(skyfits,TSTRING,"OBJECT",hout.object,"Object name",&status))  fits_report_error(stderr,status);
  if(fits_write_key(skyfits,TSTRING,"ORIGIN",hout.origin,"Origin of the image",&status))  fits_report_error(stderr,status);
  if(fits_write_key(skyfits,TSTRING,"INSTRUME",hout.instrume,"Instrument",&status))  fits_report_error(stderr,status);
  if(fits_write_key(skyfits,TFLOAT,"CRVAL1",&(hout.crval1),"First pixel of X coordinate",&status))  fits_report_error(stderr,status);
  if(fits_write_key(skyfits,TFLOAT,"CRVAL2",&(hout.crval2),"First pixel of Y coordinate",&status))  fits_report_error(stderr,status);
  if(fits_write_key(skyfits,TFLOAT,"CDELT1",&(hout.cdelt1),"Pixel size in X",&status))  fits_report_error(stderr,status);
  if(fits_write_key(skyfits,TFLOAT,"CDELT2",&(hout.cdelt2),"Pixel size in Y",&status))  fits_report_error(stderr,status);
  if(fits_write_key(skyfits,TSTRING,"CTYPE1",hout.ctype1,"Units in X axis",&status))  fits_report_error(stderr,status);
  if(fits_write_key(skyfits,TSTRING,"CTYPE2",hout.ctype2,"Units in Y axis",&status))  fits_report_error(stderr,status);
  if(fits_write_key(skyfits,TSTRING,"DATE",hout.date,"Date",&status))  fits_report_error(stderr,status);
  if(fits_write_key(skyfits,TSTRING,"FILENAME",hout.filename,"Name of this file",&status))  fits_report_error(stderr,status);
  if(fits_write_key(skyfits,TFLOAT,"DATAMIN",&(hout.datamin),"Minimum data value",&status))  fits_report_error(stderr,status);
  if(fits_write_key(skyfits,TFLOAT,"DATAMAX",&(hout.datamax),"Maximum data value",&status))  fits_report_error(stderr,status);
  fits_close_file(skyfits,&status);
  if(status) fits_report_error(stderr,status);




  f_clrl(hout.history);
  sprintf(chc,"Imagen de desviacion del cielo de %s",fitsin);
  f_actlab(hout.history,chc);
  sprintf(chc,"Subcaja [%d:%d,%d:%d]",xmin,xmax,ymin,ymax);
  f_actlab(hout.history,chc);
  sprintf(chc,"Caja para el calculo del fondo de %dx%d pixels",xscell,yscell);
  f_actlab(hout.history,chc);
  sprintf(chc,"Imagen del fondo de cielo en %s",fitssky);
  f_actlab(hout.history,chc);
  sprintf(chc,"Imagen de error de fondo de cielo en %s",fitserrsky);
  f_actlab(hout.history,chc);
  /*   fitswf(fitssig,&(hout),imasig,'T'); */

  /*Grabo con FITSIO */
  status=0;
  ffinit(&sigfits,fitssig,&status);
  if(status) fits_report_error(stderr,status);
  naxes[0]=hout.naxis1;
  naxes[1]=hout.naxis2;
  fits_create_img(sigfits,-32,2,naxes,&status);
  fits_write_img(sigfits,TFLOAT,1,naxes[0]*naxes[1],imasig,&status);
  if(fits_write_key(sigfits,TSTRING,"BUNIT",hout.bunit,"pixel units",&status))  fits_report_error(stderr,status);
  if(fits_write_key(sigfits,TSTRING,"OBJECT",hout.object,"Object name",&status))  fits_report_error(stderr,status);
  if(fits_write_key(sigfits,TSTRING,"ORIGIN",hout.origin,"Origin of the image",&status))  fits_report_error(stderr,status);
  if(fits_write_key(sigfits,TSTRING,"INSTRUME",hout.instrume,"Instrument",&status))  fits_report_error(stderr,status);
  if(fits_write_key(sigfits,TFLOAT,"CRVAL1",&(hout.crval1),"First pixel of X coordinate",&status))  fits_report_error(stderr,status);
  if(fits_write_key(sigfits,TFLOAT,"CRVAL2",&(hout.crval2),"First pixel of Y coordinate",&status))  fits_report_error(stderr,status);
  if(fits_write_key(sigfits,TFLOAT,"CDELT1",&(hout.cdelt1),"Pixel size in X",&status))  fits_report_error(stderr,status);
  if(fits_write_key(sigfits,TFLOAT,"CDELT2",&(hout.cdelt2),"Pixel size in Y",&status))  fits_report_error(stderr,status);
  if(fits_write_key(sigfits,TSTRING,"CTYPE1",hout.ctype1,"Units in X axis",&status))  fits_report_error(stderr,status);
  if(fits_write_key(sigfits,TSTRING,"CTYPE2",hout.ctype2,"Units in Y axis",&status))  fits_report_error(stderr,status);
  if(fits_write_key(sigfits,TSTRING,"DATE",hout.date,"Date",&status))  fits_report_error(stderr,status);
  if(fits_write_key(sigfits,TSTRING,"FILENAME",hout.filename,"Name of this file",&status))  fits_report_error(stderr,status);
  if(fits_write_key(sigfits,TFLOAT,"DATAMIN",&(hout.datamin),"Minimum data value",&status))  fits_report_error(stderr,status);
  if(fits_write_key(sigfits,TFLOAT,"DATAMAX",&(hout.datamax),"Maximum data value",&status))  fits_report_error(stderr,status);
  fits_close_file(sigfits,&status);
  if(status) fits_report_error(stderr,status);







  f_clrl(hout.history);
  sprintf(chc,"Imagen de desviacion del cielo de %s",fitsin);
  f_actlab(hout.history,chc);
  sprintf(chc,"Subcaja [%d:%d,%d:%d]",xmin,xmax,ymin,ymax);
  f_actlab(hout.history,chc);
  sprintf(chc,"Caja para el calculo del fondo de %dx%d pixels",xscell,yscell);
  f_actlab(hout.history,chc);
  sprintf(chc,"Imagen del fondo de cielo en %s",fitssky);
  f_actlab(hout.history,chc);
  sprintf(chc,"Imagen de error de fondo de cielo en %s",fitserrsky);
  f_actlab(hout.history,chc);
  /*   fitswf(fitserrsky,&(hout),imaerrsky,'T'); */
  
  /*Grabo con FITSIO */
  status=0;
  ffinit(&errskyfits,fitserrsky,&status);
  if(status) fits_report_error(stderr,status);
  naxes[0]=hout.naxis1;
  naxes[1]=hout.naxis2;
  fits_create_img(errskyfits,-32,2,naxes,&status);
  fits_write_img(errskyfits,TFLOAT,1,naxes[0]*naxes[1],imaerrsky,&status);
  if(fits_write_key(errskyfits,TSTRING,"BUNIT",hout.bunit,"pixel units",&status))  fits_report_error(stderr,status);
  if(fits_write_key(errskyfits,TSTRING,"OBJECT",hout.object,"Object name",&status))  fits_report_error(stderr,status);
  if(fits_write_key(errskyfits,TSTRING,"ORIGIN",hout.origin,"Origin of the image",&status))  fits_report_error(stderr,status);
  if(fits_write_key(errskyfits,TSTRING,"INSTRUME",hout.instrume,"Instrument",&status))  fits_report_error(stderr,status);
  if(fits_write_key(errskyfits,TFLOAT,"CRVAL1",&(hout.crval1),"First pixel of X coordinate",&status))  fits_report_error(stderr,status);
  if(fits_write_key(errskyfits,TFLOAT,"CRVAL2",&(hout.crval2),"First pixel of Y coordinate",&status))  fits_report_error(stderr,status);
  if(fits_write_key(errskyfits,TFLOAT,"CDELT1",&(hout.cdelt1),"Pixel size in X",&status))  fits_report_error(stderr,status);
  if(fits_write_key(errskyfits,TFLOAT,"CDELT2",&(hout.cdelt2),"Pixel size in Y",&status))  fits_report_error(stderr,status);
  if(fits_write_key(errskyfits,TSTRING,"CTYPE1",hout.ctype1,"Units in X axis",&status))  fits_report_error(stderr,status);
  if(fits_write_key(errskyfits,TSTRING,"CTYPE2",hout.ctype2,"Units in Y axis",&status))  fits_report_error(stderr,status);
  if(fits_write_key(errskyfits,TSTRING,"DATE",hout.date,"Date",&status))  fits_report_error(stderr,status);
  if(fits_write_key(errskyfits,TSTRING,"FILENAME",hout.filename,"Name of this file",&status))  fits_report_error(stderr,status);
  if(fits_write_key(errskyfits,TFLOAT,"DATAMIN",&(hout.datamin),"Minimum data value",&status))  fits_report_error(stderr,status);
  if(fits_write_key(errskyfits,TFLOAT,"DATAMAX",&(hout.datamax),"Maximum data value",&status))  fits_report_error(stderr,status);
  fits_close_file(errskyfits,&status);
  if(status) fits_report_error(stderr,status);  


  free(imasky);
  free(imasig);
  free(imaerrsky);

  printf("\n CreateSky finished\n");

  return(0);
}


/* // -------------------------------------------------------------------------- */
/* // -------------------------------------------------------------------------- */


float Sky(int n, float *ima, int nbin, float *min, float *max, float *sig,
          float ftol, float *p)

{

  float moda,hmax,delt;
  int *h,i;
  float *xp,*yp;
/*   float wxmin,wxmax,wymin,wymax; */
  int e;
/*   int emin,emax,np; */
  int iter;
  float p0[5],sigp0[5];
  float first,third;
  float median;
  float factorquartil;
  float decimin,decimax;
/*   int acum; */

/*   printf(" Calculating Sky for this subset\n"); */
  /*     printf("Antes min %f  max %f    \n",*min,*max); */
  
  /*  Calculo la moda con nbins y quartiles */
  /*  ------------------------------------- */
  /*   if((moda = StModa(n,ima,nbin,min,max)) == NULL) exit(1); */
  if(DEBUG) printf(" imamain %f imamaix %f\n",*min,*max); 
  factorquartil=15.;
  moda = StModa(n,ima,nbin,min,max);
  Quartil(n,ima,&first,&median,&third);
  if(median-(median-first)*factorquartil > *min ) *min=median-(median-first)*factorquartil;  /* Con esto se coge el 95% de la distribucion si es gaussiana */
  if(median+(third-median)*factorquartil < *max)  *max=median+(third-median)*factorquartil;  /* Si pongo 3.46 cojo el 99% */
  if(DEBUG)    printf("PASO 1 moda %f median %f first %f third %f min %f max %f\n",moda,median,first,third,*min,*max);   
  factorquartil=2.44;
  moda = StModa(n,ima,nbin,min,max);
  Decimil(n,ima,moda,&decimin,&decimax);
  if(DEBUG) printf(" Moda %f Decimin %f Decimax %f\n",moda,decimin,decimax);
  
  /*   Quartil(n,ima,&first,&median,&third); */
  /*   if(moda-(median-first)*factorquartil > *min ) *min=moda-(median-first)*factorquartil;  */ /* Con esto se coge el 95% de la distribucion si es gaussiana  */
  /*   if(moda+(third-median)*factorquartil < *max)  *max=moda+(third-median)*factorquartil;  */ /* Si pongo 3.46 cojo el 99% */ 

  *min=moda-(moda-decimin)*7;
  *max=moda+(decimax-moda)*7;

  if(DEBUG)    printf("PASO 2 moda %f median %f first %f third %f min %f max %f\n",moda,median,first,third,*min,*max);
  if(DEBUG) {
/*     *min=382728; */
/*     *max=382732; */
  }
  if(DEBUG) printf("INICIAL  Moda %f Mediana  %f Fiart %f Thisrf %f\n       Min %f   max %f\n",moda,median,first,third,*min,*max); 
  if(DEBUG) {
/*     for(i=0;i<n;i++) printf(" Valor %d : %f\n",i,ima[i]); */
  }
  
  
/*   // Calculo del histograma */
/*   // ---------------------- */
/*   //if((h=StHisto(n,ima,nbin,min,max)) == NULL) exit(1); */
  h=StHisto2(n,ima,nbin,min,max);
/*   //h=StHisto(n,ima,n,min,max); */
  delt=(*max-*min)/(nbin-1.0);
  e=(int)((moda-*min)/delt+.5);
/*   //hmax=StSuma1(n,ima,1)/nbin/2.47/(third-first)*2.0/2.3/1.5; */
  hmax=(h[(int)(nbin/2)]+h[(int)(nbin/2)+1]+h[(int)(nbin/2)-1])/3;
/*   //cpgpage(); */
  
/*   //cpgsvp(0.1,0.90,0.10,0.90); */
/*   //if(graphmod == 'F')  GraphHisto(n,ima,nbin,*min,*max,moda,*max-*min,p); */
/*   //cpglab("Pachanga","N"," " ); */
  
/*   //  printf("HOOOOOOOO min %f  max %f Moda %f Mediana   %f  Fiart %f Thisrf %f \n",*min,*max,moda,median,first,third); */
  
/*   //exit(1); */
/*   // Busco los puntos a una altura 1/2 del maximo */
/*   // -------------------------------------------- */
  /*
    i=e;
    while(i > -1 && h[i]>0.5*hmax) i--;
    emin=i;
    i=e;
    while(i < nbin && h[i]>0.5*hmax) i++;
    emax=i;
    //  printf(" emin %d emax %d\n",emin,emax);
    // Lo hago diferente. Meto tambien lo de los cuartiels aqui 
    acum=0;
    
    i=0;
    while((acum+0.)<(n/4.)) {
    acum+=h[i];
    i++;
    }
    emin=i;
    while((acum+0.)<(3*n/4.)) {
    acum+=h[i];
    i++;
    }
   */
    
  
/*   //  printf("MIO  emin %d emax %d acum %d\n",emin,emax,acum); */
  
  
  
  
/*   // Meto los datos del intervalo emin - emax (x3) en la matriz de ajuste */
/*   // -------------------------------------------------------------------- */
  
  /*
    emax=e+(emax-e)*3;
    emin=e-(e-emin)*3;
    if(emax>nbin-1) emax=nbin-1;
    if(emin<0) emin=0;
    np=emax-emin+1;
    if(np < 10) {
    fprintf(stderr,"Sky: ERROR. Muy pocos puntos para el ajuste\n");
    *sig=0;
    p[0]=0;
    p[1]=0;
    p[2]=0;
    p[3]=0;
    p[4]=0;
    return(0);
    }
   */
  
/*   // Dimensionando las matrices para el ajuste */
/*   // -----------------------------------------. */
  if((xp=malloc(nbin*sizeof(float))) == NULL) {
    fprintf(stderr,"Sky: ERROR. No puedo dimensionar la matriz xp de ");
    fprintf(stderr,"%d bytes\n",nbin*sizeof(float));
    exit(1);
    };
  if((yp=malloc(nbin*sizeof(float))) == NULL) {
    fprintf(stderr,"Sky: ERROR. No puedo dimensionar la matriz yp de ");
    fprintf(stderr,"%d bytes\n",nbin*sizeof(float));
    exit(1);
  };
  
  


/* // Solucion inicial */
/* // ---------------- */
  iter=0;
/*  // while(iter<5) { */
/* //    p0[1]=moda;				// Posicion */
/* //    p0[2]=(emax-emin)*delt/3.0;		// FWHM */
    p0[0]=median;			/* // Posicion */
    p0[1]=(third-first)/2.0/1.5;	/* // FWHM */
    p0[2]=hmax;				/* // Altura */
    p0[3]=0;				/* // a */
    p0[4]=0;				/* // b */
    sigp0[0]=p0[1]/5.;
    sigp0[1]=p0[1]/5.;
    sigp0[2]=hmax/40.;
    sigp0[3]=sqrt(n);
    sigp0[4]=2.;

/*     //cpgsci(4); */
/* //	printf("NBIN %d\n",nbin); */
  for(i=0;i<nbin; i++) { 
      xp[i] = *min + i*delt;
      yp[i]=p0[2]*exp(-M_LN2*4*(xp[i]-p0[0])*(xp[i]-p0[0])/(p0[1]*p0[1]))+p0[4]*xp[i];
  }
/*   //cpgline(nbin,xp,yp); */

/*   //cpgsci(2); */
/*   // MEto todo, se acabo. */
  /*  for(i=0; i<np; i++) {
      xp[i] = *min + (emin+i)*delt;
      yp[i] = (float)h[i+emin];
      };
   */
  for(i=0; i<nbin; i++) {
    xp[i] = *min + i*delt;
    yp[i] = (float)h[i];
/*     //cpgpt1(xp[i],yp[i],3); */
  };
/*   //cpgsci(1); */

/* //    printf("Solucion inicial.\n"); */
/* //    printf("               Maximo ... %f %f \n",p0[0],sigp0[0]); */
/* //    printf("               Posicion.. %f %f \n",p0[1],sigp0[1]); */
/* //    printf("               FWHM ..... %f %f \n",p0[2],sigp0[2]); */
/* //    printf("               a ........ %f %f \n",p0[3],sigp0[3]); */
/* //    printf("               b ........ %f %f \n",p0[4],sigp0[4]); */
/* //    printf("               Error..... %f\n",Amoe_Funk(nbin,xp,yp,p0)); */



/* // Llamada a Amoeba */
/* // ---------------- */
/*   //  iter=Amoeba(nbin,xp,yp,5,p0,sigp0,ftol,MAXITER); */
    iter=Amoeba(nbin,xp,yp,3,p0,sigp0,ftol,MAXITER,Amoe_Funk);
    if(DEBUG) {
      

    }
    printf(" Iterations: %d\n",iter); 
/*     //  } */
    *sig=0.4247*p0[1];

    
    if(iter == 0 || (p0[0] > *max || p0[0] < *min) || p0[1]>(*max-*min) ) {
/*       printf(" SIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII min %f max %f\n",*min,*max); */
      p0[0]=median;
      p0[1]=(third-first)/2.0*2.3*1.5;
    }

  
/* //  printf("Ajuste final en %d iteraciones\n",iter); */
/* //  printf("               Maximo ... %f\n",p0[0]); */
/* //  printf("               Posicion.. %f\n",p0[1]); */
/* //  printf("               FWHM ..... %f\n",p0[2]); */
/* //  printf("               a ........ %f\n",p0[3]); */
/* //  printf("               b ........ %f\n",p0[4]); */
/* //  printf("               Error..... %f\n",Amoe_Funk(nbin,xp,yp,p0)); */
  

/* //  printf("\n"); */
/* //  printf("Resultados:    Moda ..... %f\n",p0[1]); */
/* //  printf("               sigma .... %f\n",*sig); */
/* //  printf("\n"); */

/*   //cpgsci(3); */
  for(i=0;i<nbin; i++) yp[i]=p0[2]*exp(-M_LN2*4*(xp[i]-p0[0])*(xp[i]-p0[0])/(p0[1]*p0[1]))+p0[4]*xp[i];
/*   //  for(i=0;i<nbin; i++) yp[i]=p0[0]*exp(-M_LN2*4*(xp[i]-p0[1])*(xp[i]-p0[1])/(p0[2]*p0[2]))+p0[3]+p0[4]*xp[i]; */
/*   //cpgline(nbin,xp,yp); */
/*   //cpgsci(1); */




  memcpy(p,p0,5*sizeof(float));
  free(xp);
  free(yp);
  free(h);
/*   printf(" Al salir min %f max %f\n",*min,*max); */
  return(p0[0]);
} 


/* //-------------------------------------------------------------------------- */

float Amoe_Funk(int n, float *x, float *y, float *p)

{

  int i;
  float f,s;

  s=0.0;
  for(i=0; i<n; i++) {
    f=p[2]*exp(-M_LN2*4*(x[i]-p[0])*(x[i]-p[0])/(p[1]*p[1]));
/*     //    f += p[3]+p[4]*x[i]; */
/*     //f += p[4]*x[i]; */
    s += (f-y[i])*(f-y[i]);
    };
  if(DEBUG) printf(" s %f p0 %f p1 %f p2 %f\n",s,p[0],p[1],p[2]);

  return(s);
}

/* //-------------------------------------------------------------------------- */

void GraphHisto(int n, float *ima, int nbin, float min, float max, float moda, float sigma, float *p)

{

  int *h,i;
  float *xp,*yp;
  int e,emin,emax,np;
  float delt;

/*   //if((h=StHisto(n,ima,nbin,&min,&max)) == NULL) exit(1); */
/*   printf(" min %f max %f moda %f\n",min,max,moda); */
  h=StHisto(n,ima,nbin,&min,&max);
  delt=(max-min)/(nbin-1.0);
  e=(int)((moda-min)/delt+.5);
  emin=(int)((moda-3*sigma-min)/delt+.5);
  emax=(int)((moda+3*sigma-min)/delt+.5);
  if(emin < 0) emin=0;
  if(emax > nbin-1) emax=nbin-1;
  np=emax-emin+1;

  if((xp=malloc(nbin*sizeof(float))) == NULL) {
    fprintf(stderr,"Sky: ERROR. No puedo dimensionar la matriz xp de ");
    fprintf(stderr,"%d bytes\n",nbin*sizeof(float));
    exit(1);
  }
  if((yp=malloc(nbin*sizeof(float))) == NULL) {
    fprintf(stderr,"Sky: ERROR. No puedo dimensionar la matriz yp de ");
    fprintf(stderr,"%d bytes\n",nbin*sizeof(float));
    exit(1);
  }

  for(i=0; i<np; i++) {
    xp[i] = min + (emin+i)*delt;
    yp[i] = (float)h[i+emin];
  }
  for(i=0; i<nbin; i++) {
    xp[i] = min + (i)*delt;
    yp[i] = (float)h[i];
  }
  if(DEBUG) printf("AAAA Mimino %f maximo %f \n",min,max); 
  
  if(DEBUG) printf(" 222 pasa for moda %f sigma %f e %d nbin %d delt %f\n",moda,sigma,e,nbin,delt);
  
  cpgswin(moda-3.2*sigma, moda+3.2*sigma,-h[e]*.1,h[e]*1.2);
  cpgswin(min, max,-h[e]*.1,h[e]*1.2);
  if(DEBUG) printf(" Y limits: %f %f\n",-h[e]*.1,h[e]*1.2);
  if(DEBUG) cpgswin(min, max,-h[e]*.1,h[e]*1.2);

/*   cpgbox("BCTNSV",0.0,0,"bctns",0.0,0); */
  cpgbox("VBCTNS",0.0,0,"Vbctns",0.0,0);
/*   cpgbin(np  ,xp,yp,1); */
  cpgbin(nbin,xp,yp,1); 
  
  cpgsci(2);

  cpgmove(p[0],-h[e]*.1);
  cpgdraw(p[0],h[e]*1.1);

  

  for(i=0;i<np; i++) yp[i]=
      p[2]*exp(-M_LN2*4*(xp[i]-p[0])*(xp[i]-p[0])/(p[1]*p[1]));
  cpgline(np,xp,yp);
  for(i=0;i<nbin; i++) yp[i]=
      p[2]*exp(-M_LN2*4*(xp[i]-p[0])*(xp[i]-p[0])/(p[1]*p[1]));
  cpgline(nbin,xp,yp);

  if(DEBUG) printf(" AKAKA p0 %f p1 %f p2 %f \n",p[0],p[1],p[2]);

  cpgsci(1);

  free(h);
  free(xp);
  free(yp);
}
 

void Decimil(int n,float *x,float moda, float *decmin, float *decmax)
{

  float *y;
  int i, imoda;

  if((y=malloc(n*sizeof(float))) == NULL)  printf("Quartil: ERROR. Cannot dimension y of %d bytes\n ",n*sizeof(float)) ;
 
 
  memcpy(y,x,n*sizeof(float));

  qsort(y,n,sizeof(float),icompare);

  for(i=0;i<n;i++) {
    if(y[i]>moda) break;
  }
  imoda=i;

  if(DEBUG) printf(" KK imoda %d \n",imoda);

  *decmin=y[(int)(imoda-n/10)];
  *decmax=y[(int)(imoda+n/10)];
  free(y);
}


