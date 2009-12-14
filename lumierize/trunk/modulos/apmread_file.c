/*
 * Get data from the APM catalogues server. 
 *
 * This is a standalone program which makes a request of the
 * APM catalogues system based on command line arguments and returns
 * either a list or a postscript finding chart on standard output.
 *
 * VMS Notes:
 *    This program assumes Multinet TCP/IP support.
 *    This program expects to run as a foreign command.
 *    If standard output (SYS$OUTPUT) is redirected to a file (using
 *      ASSIGN or DEFINE) the resulting file will have control
 *      characters before each record.  
 *   
 * This version should compile under Ultrix and VMS/Multinet.
 *
 * Created: 1-February-1995 by T. McGlynn
 *          Goddard Space Flight Center
 *          Universities Space Research Association
 *          Code 668.1
 * Modified by Geraint Lewis and Mike Irwin to support APM online catalogues
 *          1-April-1996
 *
 */
/*#define APMCAT "131.111.68.247" - alternative form if name resolver is crap*/
/* //#define APMCAT "www.aao.gov.au" */
#define ERROR  -1

#include "modulos.h"



int apmread_file(file,cra,cdec,dra,ddec,drad,sysout,eqout,epout,mag1,mag2,classd,nstarmax,
	    gnum,gra,gdec,gmag,gtype,nlog)

     

     char file[51];
     double  cra;            /* Search center J2000 right ascension in degrees */
     double  cdec;           /* Search center J2000 declination in degrees */
     double  dra;            /* Search half width in right ascension in degrees */
     double  ddec;           /* Search half-width in declination in degrees */
     double  drad;           /* Limiting separation in degrees (ignore if 0) */
     int     sysout;         /* Search coordinate system */
     double  eqout;          /* Search coordinate equinox */
     double  epout;          /* Proper motion epoch (0.0 for no proper motion) */
     double  mag1,mag2;      /* Limiting magnitudes (none if equal) */
     int     classd;         /* Desired object class (-1=all, 0=stars, 3=nonstars) */
     int     nstarmax;       /* Maximum number of stars to be returned */
     double  *gnum;          /* Array of APM numbers (returned) */
     double  *gra;           /* Array of right ascensions (returned) */
     double  *gdec;          /* Array of declinations (returned) */
     double  *gmag;          /* Array of magnitudes (returned) */
     int     *gtype;         /* Array of object types (returned) */
     int     nlog;           /* 1 for diagnostics */
{
  
  char  w1[1000],w2[1000],w3[1000];
  int   ns=0;
  FILE *fapm;
  char line[1000];
  int nlist;
  int i;
/*   float fnul; */

/*   float rah,ram,ras; */
/*   float decg,decm; */
  float decs;
  float mag,num;
  int class;
  double ra,dec;

  nlist=FileNLin(file);
  printf(" Number of lines %d\n",nlist);

  if((fapm=fopen(file,"r")) == NULL) {
    fprintf(stderr,"APMREAD_FILE: Cannot open file %s\n", file);
    exit(1);
  }
  
  fgets(line,1000,fapm);
  fgets(line,1000,fapm);
  fgets(line,1000,fapm);
  fgets(line,1000,fapm);
  fgets(line,1000,fapm);
  fgets(line,1000,fapm);
  fgets(line,1000,fapm);
  
  ns=0;
  for(i=8;i<=nlist;i++) {
    /* scanf(" %f %f %f %f %f %f %f %d %f %f %f %f %f %f %f %f %f %f %d", */
/*         &rah,&ram,&ras,&decg,&decm,&decs,&mag,&class,&fnul,&fnul,&fnul,&fnul, */
/*         &fnul,&fnul,&fnul,&fnul,&fnul,&fnul,&num); */
    
/*     //ra=rah+ram/60+ras/3600; */
  
    fgets(line,1000,fapm);

    LeeWord(line,1,w1);LeeWord(line,2,w2);LeeWord(line,3,w3);
    ra=(atof(w1)+atof(w2)/60.+atof(w3)/3600.)*15;  
    LeeWord(line,4,w1);LeeWord(line,5,w2);LeeWord(line,6,w3);
    if(w1[0]=='-') decs=-1;  
    else decs=+1;  
/*     //	  printf(" w1 %s w2 %s w3 %s dec %d\n",w1,w2,w3,decs); */
    dec=decs*(fabs(atof(w1))+atof(w2)/60.+atof(w3)/3600.);  
    LeeWord(line,7,w1);
    mag=atof(w1);
    LeeWord(line,8,w1);
    class=(int)atof(w1);
    LeeWord(line,19,w1);
    num=atof(w1);
    if(nlog) printf(" star %d ar %f dec %f mag %f num %f type %d \n",
		      ns,gra[ns],gdec[ns],gmag[ns],gnum[ns],gtype[ns]);  
     printf(" star %d ar %f dec %f mag %f num %f type %d \n",
		     ns,ra,dec,mag,num,class);

     printf(" Criteros : ");    
     if(class == classd) {
       printf(" 111");
       if((mag>mag1 && mag<mag2) || mag1==mag2) {
	 printf(" 222");
	 printf(" cra-dra %f + %f \n",cra-dra,cra+dra);
	 if(ra>cra-dra && ra< cra+dra && dec > cdec-ddec && dec< cdec+ddec ) {
	   printf(" 333");
	   gra[ns]=ra;
	   gdec[ns]=dec;
	   gnum[ns]=num;
	   gtype[ns]=class;
	   gmag[ns]=mag;
	   ns++;  
	   if(ns>=nstarmax) {
	     fclose(fapm);
	     return (ns);
	   }
	 }
       }
     }
  } 
  fclose(fapm);
  return(ns);
}

	
