#include "modulos.h"


void LoadParam_kbd();
void LoadParam_file();
void SaveParam();
void CheckDist ();
void PlotCheck(float as);

/*Variables fichero de parametros */
char parfilename[100];
int interact=0;
/* Parametros */
char file1[51],file2[51],filematch[51];
int colxp1,colxp2,colyp1,colyp2;
int ispix=0;/*  // Esta variable comprueba si la comparacion se va a hacer utilizando astrometria (ispix=0) o las coordenadas en pixeles( ispix=1) */

int isalfadel1=0;/*  //Esta variable comprueba si la comparacion se hace directamente con las coordenadas alfa delta de los ficheros de comparacion. */
int isalfadel2=0; /* //Hacen falta dos, una para cada fichero. */
char pgdevice[51];


/* //Variables para el plot de comprobacion */
int nplot;
float ra[10000],dec[10000];
float dra[10000],ddec[10000];

double ra1,ra2,dec1,dec2;
double xp1,xp2,yp1,yp2;
double minra2,mindec2;
double maxradius;
double  dist,mindist;
/* //double distwcs; */
/* // Variables para la astrometria */
char astrfits1[51],astrfits2[51];
struct WorldCoor *wcscat1=0,*wcscat2=0;
int lhead,nbfits;  
char *header1,*header2;


int main(int argc, char **argv)
{
 

  float as=0;

  char snul[32];

  FILE *cat1,*cat2;
/*   FILE *outnomatch1,*outnomatch2; */
  FILE *outmatch;

  char c;

  
/*   char filenomatch1[51],filenomatch2[51]; */
  char nul1[1000],nul2[1000];


  int nobj1,nobj2;
  char wordra[200],worddec[200];
  int i,j;
  int jsel;
  unsigned int len1=0,len2=0;
  char record1[1000],record2[1000];
  int maxnlin;
  int foundflag;

  char option;

  /* FITS variables for FITSIO*/
  fitsfile *fitsimage,*ftemp;
  long naxes[2];
  int nfound;
  int status=0;


  if(argc<2) {
    LoadParam_kbd();
    interact=1;
  }
  else {
    strcpy(parfilename,argv[1]);
    LoadParam_file(); 
  }

/*    printf("ISPIX %d\n",ispix); */

  printf("MAXRADIUS %f\n",maxradius);

  if(!isalfadel1) {
    if(strcmp(astrfits1,"NONE")) {
      ispix=0;
      /* Open the FITS file and create a dummy file to read the header with fitsrhead */
      if( ffopen(&fitsimage, astrfits1, READONLY, &status)) fits_report_error(stderr,status);
      if(fits_read_keys_lng(fitsimage, "NAXIS", 1, 2, naxes, &nfound, &status)) fits_report_error(stderr,status);
      unlink("temp_header.fits");
      ffinit(&ftemp,"temp_header.fits",&status);
      fits_copy_header(fitsimage, ftemp, &status);
      ffclos(fitsimage,&status);
      ffclos(ftemp,&status);
      
      if ((header1 = fitsrhead ("temp_header.fits", &lhead, &nbfits)) == NULL) {
	fprintf (stderr, "Cannot read FITS header of file %s\n", astrfits1);
	exit(1);
      }
      if((wcscat1=wcsinit(header1))==NULL) {
	printf(" No WCS information found in header\n Exiting");  
	exit(1);
      }
      else {
	printf(" WCS information from header:\n");
	PrintWCS(header1,1);
      }
      unlink("temp_header.fits");
    }
    else ispix=1;
  }
  
  if(!isalfadel2) {
    if(!ispix) {
      /* Open the FITS file and create a dummy file to read the header with fitsrhead */
      if( ffopen(&fitsimage, astrfits2, READONLY, &status)) fits_report_error(stderr,status);
      if(fits_read_keys_lng(fitsimage, "NAXIS", 1, 2, naxes, &nfound, &status)) fits_report_error(stderr,status);
      unlink("temp_header.fits");
      ffinit(&ftemp,"temp_header.fits",&status);
      fits_copy_header(fitsimage, ftemp, &status);
      ffclos(fitsimage,&status);
      ffclos(ftemp,&status);
      if ((header2 = fitsrhead ("temp_header.fits", &lhead, &nbfits)) == NULL) {
	fprintf (stderr, "Cannot read FITS header of file %s\n", astrfits2);
	exit(1);
      }
      if((wcscat2=wcsinit(header2))==NULL) {
	printf(" No WCS information found in header\n Exiting");  
	exit(1);
      }
      else {
	printf(" WCS information from header:\n");
	PrintWCS(header2,1);
      }  
    }
  }


/*    printf(" ISPIX %d\n",ispix); */

  /* Reading catalog */

  nobj1=FileNLin(file1);
  nobj2=FileNLin(file2);

  if((cat1=fopen(file1,"r"))==NULL) {
    printf("\nERROR: Can't open file %s\n",file1);
    exit(1);
  };
  if((cat2=fopen(file2,"r"))==NULL) {
    printf("\nERROR: Can't open file %s\n",file2);
    exit(1);
  };
  
  outmatch=fopen(filematch,"w");
/*    outnomatch1=fopen(filenomatch1,"w"); */
/*    outnomatch2=fopen(filenomatch2,"w"); */
  fgets(nul1,1000,cat1); 
  fgets(nul2,1000,cat2); 
  len1=strlen(nul1);
  len2=strlen(nul2);
  strncpy(record1,nul1,len1-1);
  strncpy(record2,nul2,len2-1);
  
  fprintf(outmatch,"#%s          ra             dec             %s            ra             dec\n",record1,record2);
/*    fprintf(outnomatch1,"%s\n",record1); */
/*    fprintf(outnomatch2,"%s\n",record2); */
  
  if(nobj1> nobj2) maxnlin=nobj1;
  else maxnlin=nobj2;

  for (i=0;i<maxnlin;i++) { 
    fgets(nul1,1000,cat1); 
    fgets(nul2,1000,cat2); 
    if(strlen(nul1)>len1) len1=strlen(nul1);
    if(strlen(nul2)>len2) len2=strlen(nul2);
  }
  printf("\n maximas longitudes de registro: %d %d\n",len1,len2);
  rewind(cat2);
  rewind(cat1);

  /* Beginning crosscorrelation */

  nplot=0;
  printf(" Linea 000000/%6d      000000 Matchings\b\b\b\b\b\b\b\b\b\b",nobj1);
  
  for (i=0;i<nobj1;i++) { 
    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%6d/%6d      %6d",i,nobj1,nplot);
    fflush(stdout);
    fgets(nul1,1000,cat1); 
    c=nul1[0]; 
    if(c=='#') continue;
    LeeWord(nul1,colxp1,wordra);
    LeeWord(nul1,colyp1,worddec);
    if(!isalfadel1) {
      xp1=atof(wordra);
      yp1=atof(worddec);
      if(!ispix) pix2wcs(wcscat1,xp1, yp1 ,&ra1, &dec1);
      else {
	ra1=xp1;
	dec1=yp1;
      }
    }
    else {
      ra1=str2ra(wordra);
      dec1=str2dec(worddec);
    }
/*     //El -1 es porque la he cagado con el asunto de la columna 2 y 3 */

/*      printf(" xp1 %f yp1 %f ra1 %f dec1 %f\n",xp1,yp1,ra1,dec1);    */

    rewind(cat2);

    mindist=1.e30;
    foundflag=0;
    //    if(i==331) printf(" Para el 332: xp1 %f yp1 %f wordra %s wordec %s ra %f dec %f\n",xp1,yp1,wordra,worddec,ra1,dec1);
    for (j=0;j<nobj2;j++) { 
/*       printf("n1\n"); */
      fgets(nul2,1000,cat2); 
/*       printf("n2\n"); */
      c=nul2[0]; 
      if(c=='#') continue;
      LeeWord(nul2,colxp2,wordra);
      LeeWord(nul2,colyp2,worddec);
      if(!isalfadel2) {
	xp2=atof(wordra);
	yp2=atof(worddec);
	if(!ispix) pix2wcs(wcscat2,xp2, yp2 ,&ra2, &dec2);
	else {
	  ra2=xp2;
	  dec2=yp2;
	}
      }
      else {
/*       printf("n4\n"); */
	ra2=str2ra(wordra);
	dec2=str2dec(worddec);
      }

      //     if(j==78) printf(" Para el 332: xp1 %f yp1 %f wordra %s wordec %s ra %f dec %f\n",xp1,yp1,wordra,worddec,ra1,dec1);     

/*        printf(" xp2 %f yp2 %f ra2 %f dec2 %f\n",xp2,yp2,ra2,dec2);  */

/*       //Pongo la buena */
/*       //Aqui tanto ra como dec van en grados */
/*       //dist=((ra2-ra1)*(ra2-ra1)*cos((dec2+dec1)*3.1415/180./2.)+(dec2-dec1)*(dec2-dec1)); */
/*       //dist=sqrt(dist); */
      if(!ispix) dist=wcsdist(ra1,dec1,ra2,dec2);
      else dist=sqrt((ra1-ra2)*(ra1-ra2)+(dec1-dec2)*(dec1-dec2))/3600;
      
/*       printf("n5\n"); */
/*       // Esta es la buena */
/*       //      dist=((ar2-ar1)*(ar2-ar1)*cos((dec2+dec1)/2)+(dec2-dec1)*(dec2-dec1)); */
/*       //dist=sqrt(((ra2-ra1+5)*(ra2-ra1+5))+((dec2-dec1-3)*(dec2-dec1-3))); */
/*        printf(" obj 1 %f %f obj 2  %f %f  dist %f \n",ra1,dec1,ra2,dec2,3600.*dist); */
      if(3600*dist<=maxradius && 3600*dist<=mindist) {
/*       printf("n6\n"); */
	jsel=j;
/* 	//printf(" %d %d CONCINDEN %f %f %f %f\n",i+1,j+1,dist,dist*3600,maxradius,mindist); */
	/* ra2str(snul,32,ra1,3); */
/*       printf("KKKKK ra1 %s",snul); */
/*       dec2str(snul,32,dec1,3); */
/*       printf(" dec1 %s",snul); */
/*       ra2str(snul,32,ra2,3); */
/*       printf(" ra2 %s",snul); */
/*       dec2str(snul,32,dec2,3); */
/*       printf(" dec2 %s\n",snul); */

	minra2=ra2;mindec2=dec2;
	mindist=3600*dist;
        foundflag=1;
/* 	//   fprintf(outcross,"%d  %f %f %d %f %f \n",i+1,ra1,dec1,j+1,ra2,dec2); */
	len1=strlen(nul1);
	len2=strlen(nul2);
/* 	printf("n7\n"); */
	strncpy(record1,"\0",1000);
	strncpy(record2,"\0",1000);
	strncpy(record1,nul1,len1-1);
	strncpy(record2,nul2,len2-1);
/* 	//printf(" NUL1 <<%s>>\n LEN1 %d RECORD1 <<%s>>\n",nul1,len1,record1); */
/* 	//printf(" NUL2 <<%s>>\n LEN2 %d RECORD2 <<%s>>\n",nul2,len2,record2); */

/* 	//fprintf(outcross,"%s %s\n",record1,record2); */
/* 	  printf("n8\n"); */
      }
/*       printf("n9\n"); */

    }
    
    if(foundflag)  {
/*       printf("n10\n"); */
/*       //printf(" Objeto %d y %d coinciden con separacion %f\n",i+1,jsel+1,mindist); */
      if(!ispix) ra2str(snul,32,ra1,3);
/*       //printf(" ra1 %s",snul); */
      if(!ispix) dec2str(snul,32,dec1,3);
/*       //printf(" dec1 %s",snul); */
      if(!ispix) ra2str(snul,32,minra2,3);
/*       //printf(" ra2 %s",snul); */
      if(!ispix) dec2str(snul,32,mindec2,3);
/*       //printf(" dec2 %s\n",snul); */
      
/*       //fprintf(outcross," %d %d",i,j); */
      if(!ispix) {
/*       printf("n11\n"); */
        fprintf(outmatch,"%s %12.8f %12.8f %s %12.8f %12.8f  %g %g %g\n",record1,ra1/15,dec1,record2,minra2/15,mindec2,mindist,(ra1-minra2)/15.,dec1-mindec2);
        ra[nplot]=(float)(ra1/15.);dec[nplot]=(float)dec1;
        dra[nplot]=(float)(ra1-minra2)*3600*cos(dec1/180*3.1415);
        ddec[nplot]=(float)(dec1-mindec2)*3600.;
        nplot++;
      }
      else {
/*       printf("n12\n"); */
	fprintf(outmatch,"%s %12.8f %12.8f %s %12.8f %12.8f  %g %g %g\n",record1,ra1,dec1,record2,minra2,mindec2,mindist,(ra1-minra2),dec1-mindec2);
        ra[nplot]=(float)(ra1);dec[nplot]=(float)dec1;
        dra[nplot]=(float)(ra1-minra2);
        ddec[nplot]=(float)(dec1-mindec2);
        nplot++;
      }
    }

  }
  if(interact) {
    printf(" Input radius of plots (arcsec, 0=auto,<0 exit): ");
    while(as>=0) {
      as=readf(as);
      PlotCheck(as);
      printf(" Input radius of plots (arcsec, 0=auto,<0 exit): ");
    }
  }

  if(interact) {
    printf(" Do you want to save parameters in a file? [y/n]: ");
    option=readc('y'); 
    if(option=='y')   SaveParam();
  }



  return 0;
}


void CheckDist () {
  /* Aqui ar y dec van en radianes!!!!!!!!*/
  dist=((ra2-ra1)*(ra2-ra1)*cos((dec2+dec1)/2)+(dec2-dec1)*(dec2-dec1));
  if(dist<=maxradius ) {
  }
}


void PlotCheck(float as) {


  
  float drasta[10000],ddecsta[10000];
  int nsta;
  char snul[100];
  float xmin,xmax,ymin,ymax;
  int pgid;
  int i;
  float meandra,sigmadra;
  float meanddec,sigmaddec;
  pgid=cpgopen("?");

  if(as==0){
    meandra=StMedia(nplot,dra,&sigmadra);
    meanddec=StMedia(nplot,ddec,&sigmaddec);
  }
  else{
    nsta=0;
    for(i=1;i<nplot;i++) {
      if(fabs(dra[i])<=as && fabs(ddec[i])<=as) {
	drasta[nsta]=dra[i];
	ddecsta[nsta]=ddec[i];
	nsta++;
      }
    }
    meandra=StMedia(nsta,drasta,&sigmadra);
    meanddec=StMedia(nsta,ddecsta,&sigmaddec);
  }
 

  cpgask(1);
  cpgsch(.8);
  if(as==0) {
    pgLimits(nplot,dra,&xmin,&xmax);
    pgLimits(nplot,ddec,&ymin,&ymax);
  }
  else {
    xmin=-as;xmax=as;ymin=-as;ymax=as;
  }
  
  cpgsvp(0.1,0.4,0.3,0.9);
  cpgwnad(xmin,xmax,ymin,ymax);

  strcpy(snul,astrfits1);
  strcat(snul," & ");
  strcat(snul,astrfits2);
  
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpglab("RA difference (arcsec)","DEC difference (arcsec)",snul);
  cpgpt(nplot,dra,ddec,1);
  cpgsvp(0.1,0.5,0.0,0.3);
  cpgswin(0.,1.,0.,1.);

  sprintf(snul,"RA difference *cos(dec) is distributed with mean %g and stddev %g  (arcsec)",meandra,sigmadra);
  cpgptxt(0.0,0.15,0.0,0.0,snul);
  sprintf(snul,"DEC difference is distributed with mean %g and stddev %g  (arcsec)",meanddec,sigmaddec);
  cpgptxt(0.0,0.35,0.0,0.0,snul);



  if(as==0) {
    xmin=meandra-5*sigmadra;xmax=meandra+5*sigmadra;
  }

  cpgsvp(0.5,0.9,0.3,0.55);
  cpgswin(xmin,xmax,0.,nplot/2.);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpglab("RA difference * cos(dec) (arcsec)","","");
  cpghist(nplot,dra,xmin,xmax,20,1);
  printf("\n RA difference *cos(dec) is distributed with mean %g and stddev %g  (arcsec)\n",meandra,sigmadra);

  if(as==0) {
    xmin=meanddec-5*sigmaddec;xmax=meanddec+5*sigmaddec;
  }
  cpgsvp(0.5,0.9,0.65,0.9);
  cpgswin(xmin,xmax,0.,nplot/2.);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpglab("DEC difference (arcsec)","","");
  cpghist(nplot,ddec,xmin,xmax,20,1);
  printf("\n DEC difference is distributed with mean %g and stddev %g  (arcsec)\n",meanddec,sigmaddec);
  printf(" Enter <CR>\n");
  reads("",snul);
  cpgpage();  
  
  cpgsvp(0.05,0.95,0.1,0.5);
  pgLimits(nplot,ra,&xmin,&xmax);
  if(as==0)   pgLimits(nplot,dra,&ymin,&ymax);
  else {
    ymin=-as;ymax=as;
  }
  cpgswin(xmin,xmax,ymin,ymax);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpglab("RA (hours)","RA difference *cos(dec) (arcsec)","");
  cpgpt(nplot,ra,dra,3);

  cpgsvp(0.05,0.95,0.6,1.);
  pgLimits(nplot,dec,&xmin,&xmax);
  if(as==0) pgLimits(nplot,ddec,&ymin,&ymax);
  else {
    ymin=-as;ymax=as;
  }
  cpgswin(xmin,xmax,ymin,ymax);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpglab("DEC (degrees)","DEC difference (arcsec)","");
  cpgpt(nplot,dec,ddec,3);


  
  cpgclos();

}




void LoadParam_kbd() {
  printf("\n Input file with catalog 1: ");
  reads(file1,file1);
  printf(" Input column with X pixel in catalog 1 or -column for RA format: ");
  colxp1=readi(14);
  if(colxp1<0) {
    isalfadel1=1;
    colxp1=-colxp1;
  }
  if(isalfadel1) 
    printf(" Input column with DEC in catalog 1: ");
  else
    printf(" Input column with Y pixel in catalog 1: ");
  colyp1=readi(15);
  
  strcpy(astrfits1,"NONE");
  if(!isalfadel1) {
    printf(" Input FITS file with astrometric solution for catalog 1 (NONE for coordinate comparision, without WCS) : ");
    reads(astrfits1,astrfits1);
    if(strcmp(astrfits1,"NONE"))  ispix=0;
    else ispix=1;
  }
  printf("\n Input file with catalog 2: ");
  reads(file2,file2);
  printf(" Input column with X pixel in catalog 2 or -column for RA format: ");
  colxp2=readi(14);
  if(colxp2<0) {
    isalfadel2=1;
    colxp2=-colxp2;
  }
  if(isalfadel2) 
    printf(" Input column with DEC in catalog 2 (format gg mm ss): ");
  else
    printf(" Input column with Y pixel in catalog 2: ");
  colyp2=readi(15);
  strcpy(astrfits2,"NONE");
  if(!isalfadel2) {
    if(!ispix) {
      printf(" Input FITS file with astrometric solution for catalog 2: ");
      reads(astrfits2,astrfits2);
    }
  }
  if(!ispix)  printf(" Input maximum radius within the two objects can match (arcseconds): ");
  else printf(" Input maximum radius within the two objects can match (pixels): ");
  maxradius=readi(0.);
  printf(" Input output file with matching entries: ");
  reads(filematch,filematch);
/*    printf(" Input output file with non-matching entries in catalog 1: "); */
/*    reads(filenomatch1,filenomatch1); */
/*    printf(" Input output file with non-matching entries in catalog 2: "); */
/*    reads(filenomatch2,filenomatch2); */


}



void SaveParam() {
  
  FILE *fp;
  int nc=0,nt;
  char ch51[51];
  printf("Name of parameter file: ");
  reads(parfilename,parfilename);
  if((fp=fopen(parfilename,"w")) ==NULL) {
    printf("ERROR: Can't open file\n");
    return;
  }
  fprintf(fp,"COMMENT  Parameter file for CrossCat                                           \n");
  nc++;
  sprintf(ch51,"'%s'",file1    );
  fprintf(fp,"CAT1    = %-51.51s / Cat 1 to cross \n",ch51);
  sprintf(ch51,"'%s'",file2    );
  fprintf(fp,"CAT2    = %-51.51s / Cat 2 to cross \n",ch51);
  sprintf(ch51,"'%s'",filematch);
  fprintf(fp,"CROSSCAT= %-51.51s / Cat 2 to cross \n",ch51);
  fprintf(fp,"RADIUS  =%21f / Maximum radius for matching. Pixel or arcsec. \n",maxradius   );

  fprintf(fp,"XCOL1   =%21d / Column with X pix or RA (WCS1='NONE')  in Cat1\n",colxp1      );
  fprintf(fp,"YCOL1   =%21d / Column with Y pix or Dec(WCS1='NONE')  in Cat1\n",colyp1      );
  fprintf(fp,"ISRADEC1=%21d / Are cat 1 RADEC? Set =1 else pixel =0         \n",isalfadel1  );
  fprintf(fp,"XCOL2   =%21d / Column with X pix or RA (WCS2='NONE')  in Cat2\n",colxp2      );
  fprintf(fp,"YCOL2   =%21d / Column with Y pix or Dec(WCS2='NONE')  in Cat2\n",colyp2      );
  fprintf(fp,"ISRADEC2=%21d / Are cat 1 RADEC? Set =1 else pixel =0         \n",isalfadel2  );
  sprintf(ch51,"'%s'",astrfits1);
  fprintf(fp,"WCS1    = %-51.51s / WCS Image 1    \n",ch51);
  sprintf(ch51,"'%s'",astrfits2);
  fprintf(fp,"WCS2    = %-51.51s / WCS Image 2    \n",ch51);
  sprintf(ch51,"'%s'","?" );
  fprintf(fp,"DEVICE  = %-51.51s / PGPLOT device  \n",ch51);
  nc+=13;
  fprintf(fp,"COMMENT  Blank lines to assure the last block of 2880 caracters is complete    \n");
  fprintf(fp,"COMMENT                                                                        \n");
  nc+=2;
  for(nt=nc;nt<(1+(int)(nc/36))*36-1;nt++) {
    fprintf(fp,"COMMENT                                                                        \n");
  }
  fprintf(fp,"END                                                                            \n");
  fclose(fp);
  if(!isalfadel1) {
    if(strcmp(astrfits1,"NONE")) ispix=0;
    else ispix=1;
  }
}

void LoadParam_file() {

  int status=0;
  char comment[51];
  fitsfile *parfile;
  printf("Using parameter file <<%s>>\n",parfilename);
  if( ffopen2(&parfile,parfilename, READONLY, &status)) fits_report_error(stderr,status);

  ffgky(parfile,TSTRING,"CAT1"    ,file1    ,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"CAT2"    ,file2    ,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"CROSSCAT",filematch,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TDOUBLE,"RADIUS",&maxradius,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"XCOL1",&colxp1,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"YCOL1",&colyp1,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"ISRADEC1",&isalfadel1,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"XCOL2",&colxp2,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"YCOL2",&colyp2,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"ISRADEC2",&isalfadel2,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"WCS1",astrfits1,comment,&status); 
  fits_report_error(stderr,status); 
  ffgky(parfile,TSTRING,"WCS2",astrfits2,comment,&status); 
  fits_report_error(stderr,status); 
  ffgky(parfile,TSTRING,"DEVICE",pgdevice,comment,&status);
  fits_report_error(stderr,status);
  
  if(!isalfadel1 && !isalfadel2) {
    if(!strcmp(astrfits1,"NONE") && !strcmp(astrfits2,"NONE"))     ispix=1;
  }
  

  
  if(status!=0) {
    fits_report_error(stderr,status);
    printf("ERROR: Reading parameter file. Not enough keywords.\n");
    exit(1);
  }
  fits_close_file(parfile,&status);
  
  cpgopen(pgdevice);
  cpgask(0);
}




