#include "modulos.h"




/* //Variables del dibujo */
float xbmin,xbmax,ybmin,ybmax;
float xfmin,xfmax,yfmin,yfmax;

/* Varaibales para leer el FITS */
char inputfile[51];
int status=0;
int nfound, anynull;
fitsfile *image;
float datamin, datamax, nullval;
float *buffer;
long naxes[2], fpixel, nbuffer, npixels, ii,jj;

/* //Otras variables */
int *logi;
float *xp,*yp,*xep,*yep,*pa,*isoflux,*area,*aperflux;
float *iobj;
int num=0;
int elipse=1;
int  nobj=0;

/* Variables para leer el BUHO */

int xcol=14,ycol=15,xepcol=6,yepcol=7,pacol=8,fluxcol=4;
int apfluxcol=20,areacol=5,idcol=1;

/* //Declaracion de funciones */
void MainPlot();

int main()
{
  int pg1,pg2;
  long s_naxes[1], s_pixels;
  float especmin,especmax;
  float close;
  char rootspec[51];
  fitsfile *spec;
  char specfile[51];

  float *spectrum,*x;

  float xcu,ycu;



  int i;
/*   FILE *filebuho; */
  char snul[50];
  char buhofile[51];
  char opt[10];
/*   char option='C'; */
/*   float flux; */
  char cnul;
/*   float imin,imax,jmin,jmax; */
  /* Variables del dibujo */
/*   float factorarr=1; */
  float tr[6];
  float mean, sigma;
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;


  printf(" Input FITS file: ");
  scanf("%s",inputfile);
  
    
/* Leo la imagen FITS de entrada */
  if( ffopen(&image, inputfile, READONLY, &status)) fits_report_error(stderr,status);
  if(fits_read_keys_lng(image, "NAXIS", 1, 2, naxes, &nfound, &status)) fits_report_error(stderr,status);
  if(nfound!=2) {
    printf(" FITS file does not contain 2 dimensions\n");
    exit(1);
  }
  npixels=naxes[0]*naxes[1];
  if((buffer=malloc(npixels*sizeof(float)))==NULL) printf("I cannot dimension buffer of %ld elements \n",npixels);
  fpixel=1;
  nullval=0;
  datamin=1.0e30;
  datamax=-1.0e30;
  printf("...Reading image %s \n",inputfile);
  if(fits_read_img(image, TFLOAT, fpixel, npixels, &nullval, buffer, &anynull, &status )) fits_report_error(stderr,status);
  printf("...Computing datamin and datamax \n");

  for (ii=0 ;ii<npixels;ii++) {
/*     //printf(" buffer %d %e \n",ii,buffer[ii]); */
    if( buffer[ii]< datamin) datamin=buffer[ii];
    if(buffer[ii]> datamax) datamax=buffer[ii];
  }
  printf(" Datamin %f Datamax %f \n",datamin,datamax);

  mean=StMedia(npixels,buffer,&sigma);
  datamin=mean-sigma*3;
  datamax=mean+sigma*3;
  
  /* Dibujo la imagen  */

  setvbuf(stdin,"",_IOLBF,0);

/* //  pg2=cpgopen("?"); */
  cpgask(0);

  pg1=cpgopen("?");
  cpgscf(2);
  cpgask(0);
/*   //cpglab("pixel","pixel",""); */
  cpgsvp(0.05,0.95,0.05,0.95);
  cpgwnad(0.,naxes[0],0.,naxes[1]);

  
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  cpggray(buffer,naxes[0],naxes[1],1,naxes[0],1,naxes[1],datamax,datamin,tr);
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
  
  /* Leo el fichero BUHO */

  printf(" Input BUHO file: ");
  scanf("%s",buhofile);
  printf("Input column with x (pixel ) coordinate: "); 
  xcol=readi(xcol);
  printf("Input column with y (pixel ) coordinate: ");
  ycol=readi(ycol);
  printf("Input column with flux: ");
  fluxcol=readi(fluxcol);
  printf("Input column with x elip fit : ");
  xepcol=readi(xepcol);
  printf("Input column with y elip fit : ");
  yepcol=readi(yepcol);
  printf("Input column with PA elip fit : ");
  pacol=readi(pacol);
  printf("Input column with aperture flux : ");
  apfluxcol=readi(apfluxcol);
  printf("Input column with area : ");
  areacol=readi(areacol);
  printf("Input column with identification number : ");
  idcol=readi(idcol);


  nobj=FileNLin(buhofile);
  if((xp=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension xp of %d elements \n",nobj);
  if((yp=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension xp of %d elements \n",nobj);
  if((xep=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension xp of %d elements \n",nobj);
  if((yep=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension xp of %d elements \n",nobj);
  if((pa=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension xp of %d elements \n",nobj);
  if((logi=malloc(nobj*sizeof(int)))==NULL) printf("I cannot dimension xp of %d elements \n",nobj);
  if((isoflux=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension isoflux of %d elements \n",nobj);
  if((aperflux=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension aperflux of %d elements \n",nobj);
  if((area=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension area of %d elements \n",nobj);
  if((iobj=malloc(nobj*sizeof(float)))==NULL) printf("I cannot dimension iobj of %d elements \n",nobj);


  printf("Number of objects %d\n",nobj);

  printf("Root name for espectra: ");
  scanf("%s",rootspec);




  ReadNumcol(buhofile,xcol     ,xp,logi,&nobj);
  ReadNumcol(buhofile,ycol     ,yp,logi,&nobj);
  ReadNumcol(buhofile,xepcol   ,xep,logi,&nobj);
  ReadNumcol(buhofile,yepcol   ,yep,logi,&nobj);
  ReadNumcol(buhofile,pacol    ,pa,logi,&nobj);
  ReadNumcol(buhofile,fluxcol  ,isoflux,logi,&nobj);
  ReadNumcol(buhofile,apfluxcol,aperflux,logi,&nobj);
  ReadNumcol(buhofile,areacol  ,area,logi,&nobj);
  ReadNumcol(buhofile,idcol    ,iobj,logi,&nobj);
  

  /* Termino de leerlo */
  cpgsci(2);
  
  for(i=0;i<nobj;i++) {
    if(logi[i]) {
/*       //     printf(" Drawing object %f in position %f,%f\n",iobj[i]); */
      cpgelip(xp[i],yp[i],xep[i],yep[i],pa[i]);
/*       //printf("%d %f %f %f %f\n",i,xp[i],yp[i],xep[i],yep[i]); */
    }
  }
  cpgsci(1);

      xbmin=1;xbmax=naxes[0];ybmin=1;ybmax=naxes[1];

  
  
/*   //exit(1); */
 /* Menu principal */

  while(opt[0]!='E') {
    
    printf("\n\n  Z Zoom image with cursor\n");
    printf("  O Zoom out\n");
    printf("  N Write numbers [%d]\n",num);
    printf("  P Draw ellipses [%d]\n",elipse);
    printf("  S Get spectrum to nearest position\n");
    printf("  I Information about clicked object\n");
    printf("  C Change cuts max and min\n");
    printf("  E Exit\n");
    fflush(NULL);
    fflush(stdin);
    fflush(stdout);
    scanf("%s",opt);
    switch (opt[0]) {
    case 'Z' :
    case 'z' :
      cpgsci(2);
      printf(" Press bottom left square with mouse...\n");
      cpgcurs(&xbmin,&ybmin,&cnul);
      printf(" Press top right square with mouse...\n");
      cpgband(2,1,xbmin,ybmin,&xbmax,&ybmax,&cnul);
      cpgsci(1);
      break;
    case 'O' : 
    case 'o' : 
      xbmin=1;xbmax=naxes[0];ybmin=1;ybmax=naxes[1];
      break;
    case 'C' :
    case 'c' :
      printf(" Current Maximum: %f Minimum %f \n",datamax,datamin);
      printf(" New minimum: ");
      datamin=readf(datamin);
      printf(" New maximum: ");
      datamax=readf(datamax);
      break;
    case 'E' :
    case 'e' :
      cpgclos();
      exit(0);
      break;
    case 'N' :
    case 'n' :
      if(num) num=0;	
      else num=1;
      break;
    case 'P' :
    case 'p' :
      if(elipse) elipse=0;	
      else elipse=1;
      break;
    case 'I' :
    case 'i' :
      while(cnul!='X') {
	cpgcurs(&xcu,&ycu,&cnul);
	printf(" Cursor position %f,%f \n",xcu,ycu); 
	close=naxes[0]*naxes[1];
	ii=-1;
	for(i=0;i<nobj;i++) {
	  if(logi[i]) {
	    if(((xcu-xp[i]+1)*(xcu-xp[i]+1)*(ycu-yp[i]+1)*(ycu-yp[i]+1))<close) {
	      close=(xcu-xp[i]+1)*(xcu-xp[i]+1)*(ycu-yp[i]+1)*(ycu-yp[i]+1);
	      ii=i;
	      
	    }
	  }
	}
	printf(" Nearest  object: %d in position %f, %f\n",(int)(iobj[ii]),xp[ii],yp[ii]); 
	printf(" Ellipse parameters: Semimajor axis: %f Semiminor axis: %f PA: %f\n",xep[ii],yep[ii],pa[ii]);
	printf(" Isophotal flux: %f Aperture flux: %f Area: %f \n",isoflux[ii],aperflux[ii],area[ii]);
      }
      break;
    case 'S' :
    case 's' :
      while(cnul!='X') {
	cpgcurs(&xcu,&ycu,&cnul);
	printf(" Cursor position %f,%f \n",xcu,ycu); 
	close=naxes[0]*naxes[1];
	ii=-1;
	for(i=0;i<nobj;i++) {
	  if(logi[i]) {
	    if(((xcu-xp[i]+1)*(xcu-xp[i]+1)*(ycu-yp[i]+1)*(ycu-yp[i]+1))<close) {
	      close=(xcu-xp[i]+1)*(xcu-xp[i]+1)*(ycu-yp[i]+1)*(ycu-yp[i]+1);
	      ii=i;
	      
	    }
	  }
	}
	printf(" Nearest  object: %d in position %f, %f\n",(int)(iobj[ii]),xp[ii],yp[ii]); 
	strcpy(specfile,rootspec);
	sprintf(snul,"%06d",(int)iobj[ii]);
	strcat(specfile,snul);
	strcat(specfile,".fits");
	
	
	if( ffopen(&spec, specfile, READONLY, &status)) fits_report_error(stderr,status);
	if(fits_read_keys_lng(spec, "NAXIS", 1, 1, s_naxes, &nfound, &status)) fits_report_error(stderr,status);
	if(nfound!=1) {
	  printf(" FITS file does not contain 1 dimension\n");
	  exit(1);
	}
	
	s_pixels=s_naxes[0];
	if((spectrum=malloc(s_pixels*sizeof(float)))==NULL) printf("I cannot dimension spectrum of %ld elements \n",s_pixels);
	if((x=malloc(s_pixels*sizeof(float)))==NULL) printf("I cannot dimension spectrum of %ld elements \n",s_pixels);
	fpixel=1;
	nullval=0;
	printf("...Reading spectrum %s \n",specfile);
	if(fits_read_img(spec, TFLOAT, fpixel, s_pixels, &nullval, spectrum, &anynull, &status )) fits_report_error(stderr,status);
	cpgslct(pg2);
	cpgpage();
	cpgscf(2);
	
	especmin=1.0e30;
	especmax=-1.0e30;
	for (ii=0 ;ii<s_pixels;ii++) {
	  if(spectrum[ii]< especmin) especmin=spectrum[ii];
	  if(spectrum[ii]> especmax) especmax=spectrum[ii];
	}
	printf(" Especmin %f Especmax %f \n",especmin,especmax);
	for(i=0;i<s_pixels;i++) x[i]=(float)i;
	cpgswin(0.,s_pixels+1,especmin-(especmax-especmin)*.2,especmax+(especmax-especmin)*.2);
	cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
	cpgbin(s_pixels,x,spectrum,1);
	cpglab("Pixel","Flux",specfile);

	cpgslct(pg1);
	cpgcurs(&xcu,&ycu,&cnul);
	free(spectrum);free(x);
	ffclos(spec,&status);
	MainPlot();
      }
      cnul='F';
      break;
    }
    MainPlot();
  }
  




  
  return(0);
}



void MainPlot()
{
  int i;
  float tr[6];
/*   float mean, sigma; */
  char snul[50];

  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;

  
  cpgpage();
  cpgsvp(0.05,0.95,0.05,0.95);
  cpgwnad(xbmin,xbmax,ybmin,ybmax);
/*   //cpgswin(xbmin,xbmax,ybmin,ybmax); */
  
  printf(" Drawing subimage [%f:%f,%f:%f]\n",xbmin,xbmax,ybmin,ybmax);
  
  
  cpggray(buffer,naxes[0],naxes[1],xbmin,xbmax,ybmin,ybmax,datamax,datamin,tr);
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);    
  cpgsci(2);
  
  for(i=0;i<nobj;i++) {
    if(logi[i]) {
      if(elipse) {
        cpgelip(xp[i],yp[i],xep[i],yep[i],pa[i]);
        cpgpt1(xp[i],yp[i],2);
      }
      sprintf(snul,"%d",(int)iobj[i]);
      if(num) cpgtext(xp[i],yp[i],snul);
    }

  }
  cpgsci(1);
  
  
}
