#include "modulos.h"
    

/* Parametros del programa */

char dbfile[100];
char pgdevice[100];

struct SurveyDB sdb;


void LoadParam_kbd(void);
void PlotSkyMap(void);
void SelectImages(void);
void Statistics(void);
void Info(void);

int main(void) {

  char opt='E';

  LoadParam_kbd();
  ReadSurDB(dbfile, &sdb);

  
  do{
     
    printf("\n A Plot survey sky map\n"); 
    printf(" I Select images with a given criteria\n");
    printf(" S Statistics of the survey\n");
    printf(" P Print info for each image\n");
    
    printf(" E Exit\n"); 
    opt=readc(opt); 
    switch (opt) { 
    case 'A':
    case 'a':
      PlotSkyMap();
      break;
    case 'I':
    case 'i':
      SelectImages();
      break;
    case 'S':
    case 's':
      Statistics();
      break;
    case 'P':
    case 'p':
      Info();
      break;
    case 'V':
    case 'v':
      break;
    }

  }while(opt!='E' && opt!='e');

   cpgclos();

  return(0);
}

void LoadParam_kbd(void) {

  printf(" Input name of database file: ");
  reads("",dbfile);
  
  cpgopen("?");
  cpgask(0);
  cpgsch(1.5);
  cpgscf(2);

  
}



void PlotSkyMap(void) {

  float *ra,*dec;
  int i;

  static float ramin=0,ramax=0,decmin=0,decmax=0;

  double px[4],py[4];
  double pxrot[4],pyrot[4];

  double raproj[4],decproj[4];

  char cnul;
  char opt='E';

  static int labflag=0;
  
  ra =vector_f(sdb.nitems);
  dec=vector_f(sdb.nitems);

  for(i=0;i<sdb.nitems;i++) {
    ra[i] =(sdb.si[i]).alfac;
    dec[i]=(sdb.si[i]).deltac;
/*     printf(" Ta %f de %f %f %f\n",ra[i],dec[i],(sdb.si[i]).xdim,(sdb.si[i]).ydim); */
  }
  
   do{

     if(ramin==ramax && decmin==decmax) {
       MinMax(sdb.nitems,ra,&ramin,&ramax);
       MinMax(sdb.nitems,dec,&decmin,&decmax);
       ramin=ramin-(sdb.si[0]).xdim/1./3600./15.;
       ramax=ramax+(sdb.si[0]).xdim/1./3600./15.;
       decmin=decmin-(sdb.si[0]).ydim/1./3600.;
       decmax=decmax+(sdb.si[0]).ydim/1./3600.;    
     }
      
     cpgpage();
     
     cpgvstd();
     cpgwnad(ramax*3600*15.,ramin*3600*15.,decmin*3600,decmax*3600);
     cpgswin(ramax*3600,ramin*3600,decmin*3600,decmax*3600);
     cpgswin(ramax*3600,ramin*3600,decmin*3600,decmax*3600);
     /*   cpgswin(ramax,ramin,decmin,decmax); */
     cpgtbox("ZXBCTNSHG",0.0,0,"ZBCTNSDGV",0.0,0);
     cpgswin(ramax,ramin,decmin,decmax);
     
     for(i=0;i<sdb.nitems;i++) {
       
       px[0]=-(sdb.si[i]).xdim/3600./2./180*M_PI;py[0]=-(sdb.si[i]).ydim/3600./2./180*M_PI;
       px[1]=+(sdb.si[i]).xdim/3600./2./180*M_PI;py[1]=-(sdb.si[i]).ydim/3600./2./180*M_PI;
       px[2]=+(sdb.si[i]).xdim/3600./2./180*M_PI;py[2]=+(sdb.si[i]).ydim/3600./2./180*M_PI;
       px[3]=-(sdb.si[i]).xdim/3600./2./180*M_PI;py[3]=+(sdb.si[i]).ydim/3600./2./180*M_PI;
       pxrot[0]=(sdb.si[i]).alfac + (px[0]*cos((sdb.si[i]).rot/180*M_PI)-py[0]*sin((sdb.si[i]).rot/180*M_PI));
       pyrot[0]=(sdb.si[i]).deltac+ px[0]*sin((sdb.si[i]).rot/180*M_PI)+py[0]*cos((sdb.si[i]).rot/180*M_PI);
       pxrot[1]=(sdb.si[i]).alfac + (px[1]*cos((sdb.si[i]).rot/180*M_PI)-py[1]*sin((sdb.si[i]).rot/180*M_PI));
       pyrot[1]=(sdb.si[i]).deltac+ px[1]*sin((sdb.si[i]).rot/180*M_PI)+py[1]*cos((sdb.si[i]).rot/180*M_PI);
       pxrot[2]=(sdb.si[i]).alfac + (px[2]*cos((sdb.si[i]).rot/180*M_PI)-py[2]*sin((sdb.si[i]).rot/180*M_PI));
       pyrot[2]=(sdb.si[i]).deltac+ px[2]*sin((sdb.si[i]).rot/180*M_PI)+py[2]*cos((sdb.si[i]).rot/180*M_PI);
       pxrot[3]=(sdb.si[i]).alfac + (px[3]*cos((sdb.si[i]).rot/180*M_PI)-py[3]*sin((sdb.si[i]).rot/180*M_PI));
       pyrot[3]=(sdb.si[i]).deltac+ px[3]*sin((sdb.si[i]).rot/180*M_PI)+py[3]*cos((sdb.si[i]).rot/180*M_PI);
       Plac2Ecu(pxrot[0],pyrot[0],(sdb.si[i]).alfac/180*M_PI*15.,(sdb.si[i]).deltac/180*M_PI,raproj+0,decproj+0);
       Plac2Ecu(pxrot[1],pyrot[1],(sdb.si[i]).alfac/180*M_PI*15.,(sdb.si[i]).deltac/180*M_PI,raproj+1,decproj+1);
       Plac2Ecu(pxrot[2],pyrot[2],(sdb.si[i]).alfac/180*M_PI*15.,(sdb.si[i]).deltac/180*M_PI,raproj+2,decproj+2);
       Plac2Ecu(pxrot[3],pyrot[3],(sdb.si[i]).alfac/180*M_PI*15.,(sdb.si[i]).deltac/180*M_PI,raproj+3,decproj+3);
       raproj[0]*=180/M_PI/15; decproj[0]*=180/M_PI;
       raproj[1]*=180/M_PI/15; decproj[1]*=180/M_PI;
       raproj[2]*=180/M_PI/15; decproj[2]*=180/M_PI;
       raproj[3]*=180/M_PI/15; decproj[3]*=180/M_PI;
       
       /*        printf(" ra %f dec %f p0 %f %f p1 %f %f p2 %f %f p3 %f %f\n",(sdb.si[i]).alfac,(sdb.si[i]).deltac,pxrot[0],pyrot[0],pxrot[1],pyrot[1],pxrot[2],pyrot[2],pxrot[3],pyrot[3]); */


       
       cpgsfs(1);
       cpgsci(0);
       cpgpoly_d(4,raproj,decproj);
       cpgsci(2);
       cpgsfs(2);
       cpgpoly_d(4,raproj,decproj);
       cpgsci(1);
       if(labflag) cpgtext((pxrot[0]+pxrot[1]+pxrot[2]+pxrot[3])/4,(pyrot[0]+pyrot[1]+pyrot[2]+pyrot[3])/4,sdb.si[i].image);
/*        labflag=readi(labflag); */
       
     }
     
     cpglab("RA","","");
     
     cpgptxt(ramin-(ramax-ramin)*0.08,decmin+(decmax-decmin)/2.,90,0.5,"Dec");
     
     printf("\n Z Zoom with cursor\n"); 
     printf(" O Zoom out\n");
     printf(" L Print Labels for each image\n");
     printf(" E Exit\n");
     
     opt=readc(opt); 
     switch (opt) { 
     case 'Z':
     case 'z':
       cpgsci(2);
       printf(" Press bottom left square with mouse...\n");
       cpgcurs(&ramin,&decmin,&cnul);
       printf(" Press top right square with mouse...\n");
       cpgband(2,1,ramin,decmin,&ramax,&decmax,&cnul);
       cpgsci(1);
       break;
     case 'O':
     case 'o':
       ramin=0;
       ramax=0;
       decmin=0;
       decmax=0;
       break;
     case 'L':
     case 'l':
       labflag=1-labflag;
       break;
       
     }
     
   }while(opt!='E' && opt!='e');
   
   
   free(ra);
   free(dec);
}



void SelectImages(void) {

/*   int i; */

  char opt='E';
  char outputfile[101];
/*   FILE *of; */
  
  float minseeing,maxseeing;
  
  static int flagseeing=0;
/*   char wcsar[16],wcsdec[16]; */
 
  
  printf("\n S Select acording with seeing\n"); 
  printf(" E Exit\n");
  
  opt=readc(opt); 
  
  switch (opt) { 
  case 'S':
  case 's':
    printf(" Input minimum seeing: ");
    minseeing=readf(0.);
    printf(" Input maximum seeing: ");
    maxseeing=readf(0.);
    flagseeing=1;
    break;
  case 'E':
  case 'e':
    return;
    break;
  }
  
  printf(" Output register file: ");
  reads("",outputfile);

  SaveSurDB(outputfile,sdb);
}



void Statistics(void) 
{
  float area;

  printf(" Number of images: %d \n",sdb.nitems); 
  area=Surveyrad(sdb);
  printf(" Square degrees surveyed: %f    square rads: %f\n",area*180*180/M_PI/M_PI,area);
   
} 

void Info(void)
{
  int i;
  printf("%-51s %-12s  %-12s\n","Image","Alpha","Delta");
  for(i=0;i<sdb.nitems;i++) {
    printf("%-51s %12f  %12f\n",sdb.si[i].image,sdb.si[i].alfac,sdb.si[i].deltac);
  }

  

}
