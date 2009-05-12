#include "modulos.h"



void CheckDist ();
void PlotCheck(float as);

//Variables para el plot de comprobacion
int nplot;
float ra[10000],dec[10000];
float dra[10000],ddec[10000];

double ra1,ra2,dec1,dec2;
double xp1,xp2,yp1,yp2;
double minra2,mindec2;
double maxradius;
double  dist,mindist;
//double distwcs;
// Variables para la astrometria
char astrfits1[51],astrfits2[51];
struct WorldCoor *wcscat1=0,*wcscat2=0;
int lhead,nbfits;  
char *header1,*header2;
int isalfadel1=0,isalfadel2=0;

int  main()
{
 
  int ispix=1; // Esta variable comprueba si la comparacion se va a hacer utilizando astrometria (ispix=0) o las coordenadas en pixeles( ispix=1)



  float as=0;

  char snul[32];

  FILE *cat1,*cat2;
  FILE *outcross;

  char c;

  
  char file1[51],file2[51],fileout[51];
  char nul1[1000],nul2[1000];

  int colxp1,colxp2,colyp1,colyp2;
  int nobj1,nobj2;
  char wordra[200],worddec[200];
  int i,j;
  int jsel;
  unsigned int len1=0,len2=0;
  char record1[1000],record2[1000];
  int maxnlin;
  int foundflag;



/*   printf(" input coufla\b"); */
/*   maxradius=readd(10); */
/*   printf(" maxradius %g\n"); */
/*   maxradius=readd(10.23); */
/*   printf(" maxradius %g\n"); */
/*   maxradius=readd(0); */
/*   printf(" maxradius %g\n"); */
/*   exit(1); */

  
  //  strcpy(nul2,"ESTOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n");
  //  printf("nul2 <<%s>>\n");
  //  strncpy(record2,nul2,strlen(nul2)-1);

  //  printf("record2 <<%s>>\n");
  //  exit(1);


  // CUIDADO!! En realidad, hay que leer las columnas 14 y 15 de BuhoFits.
  // La astrometria estara hecha (WFC cagada) con 2 y 3. Esto hay que cambiarlo




  printf("\n Input file with catalog 1: ");
  colxp1=14;
  colyp1=15;
  printf("\n Input file with catalog 2: ");
  colxp2=14;
  colyp2=15;
  printf(" Input tolerance radius within the two objects can match (arcseconds): ");
  scanf("%lf",&maxradius);
  printf(" Input output file with offsets to apply: ");
  scanf("%s",fileout);

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
  
  outcross=fopen(fileout,"w");
  fgets(nul1,1000,cat1); 
  fgets(nul2,1000,cat2); 
  len1=strlen(nul1);
  len2=strlen(nul2);
  strncpy(record1,nul1,len1-1);
  strncpy(record2,nul2,len2-1);
  
  fprintf(outcross,"%s          ra             dec             %s            ra             dec\n",record1,record2);
  
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

  nplot=0;
  printf(" Linea 000000/%6d      000000 Matchings\b\b\b\b\b\b\b\b\b\b",nobj1);
  
  for (i=0;i<nobj1;i++) { 
    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%6d/%6d      %6d",i,nobj1,nplot);
    fflush(stdout);
    fgets(nul1,1000,cat1); 
    c=nul1[0]; 
    if(c=='#') continue;
    //LeeWord(nul1,colra1,word);
    LeeWord(nul1,colxp1,wordra);
    //ra1=atof(word);
    xp1=atof(wordra);
    //LeeWord(nul1,coldec1,word);
    LeeWord(nul1,colyp1,worddec);
    //dec1=atof(word);
    if(!isalfadel1) {
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
      //printf(" Coor 1: ra %g  dec %g\n",ra1,dec1);
    }
    //El -1 es porque la he cagado con el asunto de la columna 2 y 3

    rewind(cat2);

    mindist=1.e30;
    foundflag=0;
    for (j=0;j<nobj2;j++) { 
      fgets(nul2,1000,cat2); 
      c=nul2[0]; 
      if(c=='#') continue;
/*       LeeWord(nul2,colra2,word); */
/*       ra2=atof(word); */
/*       LeeWord(nul2,coldec2,word); */
/*       dec2=atof(word); */
      LeeWord(nul2,colxp2,wordra);
      xp2=atof(wordra);
      LeeWord(nul2,colyp2,worddec);
      yp2=atof(worddec);
      if(!isalfadel2) {
	if(!ispix) pix2wcs(wcscat2,xp2, yp2 ,&ra2, &dec2);
	else {
	  ra2=xp2;
	  dec2=yp2;
	}
      }
      else {
	ra2=str2ra(wordra);
	dec2=str2dec(worddec);
      }
      //Pongo la buena
      //Aqui tanto ra como dec van en grados
      //dist=((ra2-ra1)*(ra2-ra1)*cos((dec2+dec1)*3.1415/180./2.)+(dec2-dec1)*(dec2-dec1));
      //dist=sqrt(dist);
      if(!ispix) dist=wcsdist(ra1,dec1,ra2,dec2);
      else dist=sqrt((ra1-ra2)*(ra1-ra2)+(dec1-dec2)*(dec1-dec2))/3600;
      // Esta es la buena
      //      dist=((ar2-ar1)*(ar2-ar1)*cos((dec2+dec1)/2)+(dec2-dec1)*(dec2-dec1));
      //dist=sqrt(((ra2-ra1+5)*(ra2-ra1+5))+((dec2-dec1-3)*(dec2-dec1-3)));
      //printf(" obj 1 %f %f obj 2  %f %f  dist %f \n",ra1,dec1,ra2,dec2,dist);
      if(3600*dist<=maxradius && 3600*dist<=mindist) {
	jsel=j;
	//printf(" %d %d CONCINDEN %f %f %f %f\n",i+1,j+1,dist,dist*3600,maxradius,mindist);
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
	//   fprintf(outcross,"%d  %f %f %d %f %f \n",i+1,ra1,dec1,j+1,ra2,dec2);
	len1=strlen(nul1);
	len2=strlen(nul2);
	strncpy(record1,"\0",1000);
	strncpy(record2,"\0",1000);
	strncpy(record1,nul1,len1-1);
	strncpy(record2,nul2,len2-1);
	//printf(" NUL1 <<%s>>\n LEN1 %d RECORD1 <<%s>>\n",nul1,len1,record1);
	//printf(" NUL2 <<%s>>\n LEN2 %d RECORD2 <<%s>>\n",nul2,len2,record2);

	//fprintf(outcross,"%s %s\n",record1,record2);
      }


    }
    
    if(foundflag)  {
      //printf(" Objeto %d y %d coinciden con separacion %f\n",i+1,jsel+1,mindist);
      if(!ispix) ra2str(snul,32,ra1,3);
      //printf(" ra1 %s",snul);
      if(!ispix) dec2str(snul,32,dec1,3);
      //printf(" dec1 %s",snul);
      if(!ispix) ra2str(snul,32,minra2,3);
      //printf(" ra2 %s",snul);
      if(!ispix) dec2str(snul,32,mindec2,3);
      //printf(" dec2 %s\n",snul);
      
      //fprintf(outcross," %d %d",i,j);
      if(!ispix) {
        fprintf(outcross,"%s %12.8f %12.8f %s %12.8f %12.8f  %f %f %f\n",record1,ra1/15,dec1,record2,minra2/15,mindec2,mindist,(ra1-minra2)/15,dec1-mindec2);
        ra[nplot]=(float)(ra1/15);dec[nplot]=(float)dec1;
        dra[nplot]=(float)(ra1-minra2)*3600*cos(dec1/180*3.1415);
        ddec[nplot]=(float)(dec1-mindec2)*3600;
        nplot++;
      }
      else {
	fprintf(outcross,"%s %12.8f %12.8f %s %12.8f %12.8f  %f %f %f\n",record1,ra1,dec1,record2,minra2,mindec2,mindist,(ra1-minra2),dec1-mindec2);
        ra[nplot]=(float)(ra1);dec[nplot]=(float)dec1;
        dra[nplot]=(float)(ra1-minra2);
        ddec[nplot]=(float)(dec1-mindec2);
        nplot++;
      }
    }

  }
  while(as>=0) {
    printf(" Input radius of plots (arcsec, 0=auto,<0 exit): ");
    as=readf(as);
    PlotCheck(as);
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

  

  cpgask(1);
  if(as==0) {
    pgLimits(nplot,dra,&xmin,&xmax);
    pgLimits(nplot,ddec,&ymin,&ymax);
  }
  else {
    xmin=-as;xmax=as;ymin=-as;ymax=as;
  }
  
  cpgsvp(0.1,0.4,0.1,0.9);
  cpgwnad(xmin,xmax,ymin,ymax);

  strcpy(snul,astrfits1);
  strcat(snul," & ");
  strcat(snul,astrfits2);
  
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpglab("RA difference (arcsec)","DEC difference (arcsec)",snul);
  
  cpgpt(nplot,dra,ddec,1);
  //for(i=0;i<nplot;i++)    cpgarro(0.,0.,dra[i],ddec[i]);
  
  
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
    meandra=StMedia(nplot,drasta,&sigmadra);
    meanddec=StMedia(nplot,ddecsta,&sigmaddec);
  }
 
  if(as==0) {
    xmin=meandra-5*sigmadra;xmax=meandra+5*sigmadra;
  }

  cpgsvp(0.5,0.9,0.1,0.5);
  cpgswin(xmin,xmax,0.,nplot/2.);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpglab("RA difference * cos(dec) (arcsec)","","");
  cpghist(nplot,dra,xmin,xmax,20,1);
  printf("\n RA difference *cos(dec) is distributed with mean %f and stddev %f  (arcsec)\n",meandra,sigmadra);
  if(as==0) {
    xmin=meanddec-5*sigmaddec;xmax=meanddec+5*sigmaddec;
  }
  cpgsvp(0.5,0.9,0.6,1.);
  cpgswin(xmin,xmax,0.,nplot/2.);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpglab("DEC difference (arcsec)","","");
  cpghist(nplot,ddec,xmin,xmax,20,1);
  printf("\n DEC difference is distributed with mean %f and stddev %f  (arcsec)\n",meanddec,sigmaddec);

  printf(" Enter <CR>\n");
  reads(snul,snul);
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
