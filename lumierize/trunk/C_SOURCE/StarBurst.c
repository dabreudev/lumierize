#include "modulos.h"
#define NMAXDATAF 10
#define DEBUG 0
#define DEBUG2 0

struct CMHKgrid {
  float fOIII;
  float fOII;
  float fHa;
  float fHb;
  float fNII;
  float fOI;
  float ewHb;
};

struct CMHKmodel {
  
  int nmet;
  int nage;
  
  float *age;
  float *z;
  struct CMHKgrid  Ginst[20][1000];
  struct CMHKgrid Gconst[20][1000];

};



struct gal {
  int   dataset;

  float fOIII,efOIII;
  float fOII,efOII;
  float fHa,efHa;
  float fHb,efHb;
  float fNII,efNII;
  float fOI,efOI;
  float fHeI,efHeI;
  float fSII_1,efSII_1;
  float fSII_2,efSII_2;
  float ewHa,eewHa;
  float ewHb,eewHb;
  char  name[51];
  int   ifOIII,iefOIII;
  int   ifOII,iefOII;
  int   ifHa,iefHa;
  int   ifHb,iefHb;
  int   ifNII,iefNII;
  int   ifOI,iefOI;
  int   ifHeI,iefHeI;
  int   ifSII_1,iefSII_1;
  int   ifSII_2,iefSII_2;
  int   iewHa,ieewHa;
  int   iewHb,ieewHb;

  float x,y;
};

struct data {
  char  datafile[51];
  char  label[51];

  char  *name;
  float *fOIII,*efOIII;
  float *fOII,*efOII;
  float *fHa,*efHa;
  float *fHb,*efHb;
  float *fNII,*efNII;
  float *fOI,*efOI;
  float *fHeI,*efHeI;
  float *fSII_1,*efSII_1;
  float *fSII_2,*efSII_2;
  float *ewHa,*eewHa;
  float *ewHb,*eewHb;
  int   *iname;
  int   *ifOIII,*iefOIII;
  int   *ifOII,*iefOII;
  int   *ifHa,*iefHa;
  int   *ifHb,*iefHb;
  int   *ifNII,*iefNII;
  int   *ifOI,*iefOI;
  int   *ifHeI,*iefHeI;
  int   *ifSII_1,*iefSII_1;
  int   *ifSII_2,*iefSII_2;
  int   *iewHa,*ieewHa;
  int   *iewHb,*ieewHb;
  
  int   nobj;
  int   ncat;
};
struct extlaw
{
  float *ldo;
  float *y;
  int n;
};



/* Parameters */

int interact=1; 
int errflag=0;
int CMHKflag=0;
int modelflag=0;

char parfilename[100];
char extinfile[100];
char pgdevice[51];
char CMHKdir[51];


int fOIIIcol=26,fOIIIerrcol=27;
int fOIIcol=4,fOIIerrcol=5;
int fHacol=32,fHaerrcol=33;
int fHbcol=22,fHberrcol=23;
int fNIIcol=34,fNIIerrcol=35;
int fOIcol=30,fOIerrcol=31;
int fHeIcol=28,fHeIerrcol=29;
int fSII_1col=38,fSII_1errcol=39;
int fSII_2col=40,fSII_2errcol=41;


int ewHacol=75,ewHaerrcol=76;
int ewHbcol=65,ewHberrcol=66;


int objcol;

/*Data variables */

struct gal *G;
struct data D[NMAXDATAF];
struct extlaw EX;
struct CMHKmodel CMHK;

int ngal;
int nfiles;


void MainMenu(void);
void PlotData(int nplot);
void PlotCMHK(int nplot);
void ReadCat(void);
void SaveParam(void);
void LoadParam_file(void);
void LoadParam_kbd(void);
void ReadExtlaw(void);
void ReadCMHK(void);
float EXT(float ldo);


int main(int argc, char **argv)
{
  if(argc<2) {
    LoadParam_kbd();
    interact=1;
  }
  else {
    strcpy(parfilename,argv[1]);
    LoadParam_file(); 
  }
  printf("...Reading data files\n");
  ReadCat();
  printf("...DONE\n");
  if(CMHKflag) {
    printf("...Reading models\n");
    ReadCMHK();
  }
  printf("...DONE\n");

  ReadExtlaw();

  MainMenu();
  
  cpgclos();  
  if(interact)  SaveParam();
  
  return 0;
}

void MainMenu(void) {
  
  int nplot=1;
  
  while(nplot!=-1) {
    
    printf("-1  Exit\n");
    printf(" 0  Plotting options\n");
    
    printf(" 1  [OIII]5007/Hb             vs      [NII]6583/Ha\n");
    printf(" 2  [OIII]5007/Hb             vs      [OI]6300/Ha\n");
    printf(" 3  [OIII]5007/Hb             vs      [OII]3727/[OIII]\n");
    printf(" 4  [OIII]5007/Hb             vs      [OII]3727/Hb\n");
    printf(" 5  [OIII]5007/Hb             vs      [SII]6717+6731/Hb\n");
    printf(" 6  [NII]6583/[OII]3727       vs      ([OII]3727+[OIII]5007)/Hb\n");
    
    printf(" 7  [OI]6300/Hb               vs      EW(Hb)\n");
    printf(" 8  [HeI]5876/Hb              vs      EW(Hb)\n");
    printf(" 9  [OIII]5007/Hb             vs      EW(Hb)\n");
    printf(" 10 [OII]3727/Hb              vs      EW(Hb)\n");
    printf(" 11 ([OII]3727+[OIII]5007)/Hb vs      EW(Hb)\n");
    printf(" 12 [NII]6583/[OII]3727       vs      EW(Hb)\n");
    
    
    
    printf(" Type number of diagnostic diagram to plot: ");
    nplot=readi(nplot);
    
    if(nplot==-1 && nplot<13);
    else    {
      PlotData(nplot);
      if(CMHKflag) PlotCMHK(nplot);
    }
    
  }
  
  
}




void PlotData(int nplot) {

  int i;
  int j;

  float dum1,dum2,dum3;
  
  float *x,*errx;
  float *y,*erry;
  int   *id;
  int   *ids;
  char  *labels;
  float xmin=1,xmax=1;
  float ymin=0,ymax=0;
  int nlines=0;

  float logxextin=0,logyextin=0;
  
  char xlabel[100];
  char ylabel[100];
  
  static int labflag[NMAXDATAF]={0,0,0,0,0,0,0,0,0,0};
  static int   color[NMAXDATAF]={1,1,1,1,1,1,1,1,1,1};
  static int  symbol[NMAXDATAF]={6,7,8,9,10,11,12,13,14,15};
  static int ds=0;
  static int mo=1;
  static int ic;
  
  if(DEBUG) printf(" Dentro de main\n");
  
  if((x   =malloc(ngal*2*sizeof(float    )))==NULL) printf("I cannot dimension x       of %d elements \n",ngal  );
  if((y   =malloc(ngal*2*sizeof(float    )))==NULL) printf("I cannot dimension y       of %d elements \n",ngal  );
  if((errx=malloc(ngal*2*sizeof(float    )))==NULL) printf("I cannot dimension errx    of %d elements \n",ngal  );
  if((erry=malloc(ngal*2*sizeof(float    )))==NULL) printf("I cannot dimension erry    of %d elements \n",ngal  );
  if((id  =malloc(ngal*2*sizeof(int      )))==NULL) printf("I cannot dimension id      of %d elements \n",ngal );
  if((ids =malloc(ngal*2*sizeof(int      )))==NULL) printf("I cannot dimension ids     of %d elements \n",ngal );
  if((labels=malloc(ngal*2*51*sizeof(char  )))==NULL) printf("I cannot dimension labels     of %d elements \n",ngal );
  
  if(DEBUG) printf(" Antes swtcih\n");
  switch(nplot) {
  case 0:
    printf("\n");
    printf(" 1 Change data   plot settings\n");
    printf(" 2 Change models plot settings\n\n");
    

    ic=readi(ic);
    switch(ic) {
    case 1:
      for(i=0;i<nfiles;i++)  printf(" %-2d %s\n",i+1,D[i].datafile);      
      printf("\n Enter number of dataset to change (1-%d)\n",nfiles);
      ds=readi(ds+1);
      ds--;
      printf(" 1 Change label on/off\n");
      printf(" 2 Change symbols\n");
      printf(" 3 Change color\n");
      printf("\n What do you want to change?\n");
      ic=readi(ic);
      switch(ic) {
      case 1:
	printf(" Labels?\n");
	labflag[ds]=readi(labflag[ds]);
	break;
      case 2:
	printf(" Input symbol\n");
	symbol[ds]=readi(symbol[ds]);
	break;
      case 3:
	printf(" Input color\n");
	color[ds]=readi(color[ds]);
	break;
      }
      break;
    case 2:
      printf(" 1 CMHK models\n");
      
      printf("\n Input model to change\n");
      mo=readi(mo);
      switch(mo) {
      case 1:
	printf(" 1 Plot on/off (%d)\n",CMHKflag);
	
	printf("\n Input option\n");
	ic=readi(ic);
	switch(ic) {
	case 1:
	  printf(" Plot model CMHK? (0=no, 1=yes)\n");
	  CMHKflag=readi(CMHKflag);
	  break;
	}
	break;
      }
      break;
    }
    break;
    
  case 1:
    xmin=-2.3;xmax=1;
    ymin=-1.0;ymax=1.4;
    strcpy(xlabel,"log([NII]6583/H\\ga)");
    strcpy(ylabel,"log([OIII]5007/H\\gb)");
    strcpy(xlabel,"[NII]6583/H\\ga");
    strcpy(ylabel,"[OIII]5007/H\\gb");
    if(DEBUG) printf(" Estoy en uno\n");
    nlines=0;
    for(j=0;j<ngal;j++) {
      if((G[j]).ifOIII && (G[j]).ifHb && (G[j]).ifHa && (G[j]).ifNII && (G[j]).fHb!=0 && (G[j]).fHa!=0) {
	G[j].x=(G[j]).fNII/(G[j]).fHa;
	G[j].y=(G[j]).fOIII/(G[j]).fHb;
	x[nlines]=G[j].x;
	y[nlines]=G[j].y;
	if(DEBUG) printf(" antes err\n");
	if(errflag) {
	  errx[nlines]=sqrt((G[j]).efNII/(G[j]).fHa*(G[j]).efNII/(G[j]).fHa+(G[j]).fNII/((G[j]).fHa*(G[j]).fHa)*(G[j]).fNII/((G[j]).fHa*(G[j]).fHa)*(G[j]).efHa*(G[j]).efHa);
	  erry[nlines]=sqrt((G[j]).efOIII/(G[j]).fHb*(G[j]).efOIII/(G[j]).fHb+(G[j]).fOIII/((G[j]).fHb*(G[j]).fHb)*(G[j]).fOIII/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb);
 	  if(DEBUG) printf(" error %s   %g %g   %g %g %g   %g  %g    %g %g %g \n",G[j].name,(G[j]).fOII,(G[j]).fOIII,(G[j]).efOII,(G[j]).efOIII,(G[j]).efHb,(G[j]).efNII,(G[j]).efHa,(G[j]).efOIII/(G[j]).fHb*(G[j]).efOIII/(G[j]).fHb,(G[j]).fOIII/((G[j]).fHb*(G[j]).fHb)*(G[j]).fOIII/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb,erry[nlines]); 
	}
	if(DEBUG) printf(" desd err\n");
	id[nlines]=j;	  ids[nlines]=G[j].dataset;
	if(DEBUG) printf(" lo id pas\n");
	strcpy(labels+nlines*51,G[j].name);
	if(DEBUG) printf(" y altr\n");
/* 	if(DEBUG) printf("G %d   x %g y %g    %s %s\n",j,G[j].x,G[j].y,G[j].name,labels[nlines*51]); */
/* 	if(DEBUG) printf(" OIII %g Hb %g Ha %g NII %g\n",G[j].fOIII,G[j].fHb,G[j].fHa,G[j].fNII);	   */
	nlines++;
	logxextin=1.47*(EXT(6583.)-EXT(6563.));
	logyextin=1.47*(EXT(5007.)-EXT(4861.));
      }
      
    }
    if(DEBUG) printf(" Ya vio y %d\n",nlines);
    break;
  case 2:
    xmin=-4.0;xmax=0;
    ymin=-1.0;ymax=1.4;
    strcpy(xlabel,"log([OI]6300/H\\ga)");
    strcpy(ylabel,"log([OIII]5007/H\\gb)");
    strcpy(xlabel,"[OI]6300/H\\ga");
    strcpy(ylabel,"[OIII]5007/H\\gb");
    if(DEBUG) printf(" Estoy en uno\n");
    nlines=0;
    for(j=0;j<ngal;j++) {
      if((G[j]).ifOIII && (G[j]).ifHb && (G[j]).ifHa && (G[j]).ifOI && (G[j]).fHb!=0 && (G[j]).fHa!=0) {
	G[j].x=(G[j]).fOI/(G[j]).fHa;
	G[j].y=(G[j]).fOIII/(G[j]).fHb;
	x[nlines]=G[j].x;
	y[nlines]=G[j].y;
	if(errflag) {
	  errx[nlines]=sqrt((G[j]).efOI/(G[j]).fHa*(G[j]).efOI/(G[j]).fHa+(G[j]).fOI/((G[j]).fHa*(G[j]).fHa)*(G[j]).fOI/((G[j]).fHa*(G[j]).fHa)*(G[j]).efHa*(G[j]).efHa);
	  erry[nlines]=sqrt((G[j]).efOIII/(G[j]).fHb*(G[j]).efOIII/(G[j]).fHb+(G[j]).fOIII/((G[j]).fHb*(G[j]).fHb)*(G[j]).fOIII/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb);
	}
	id[nlines]=j;	  ids[nlines]=G[j].dataset;
	strcpy(labels+nlines*51,G[j].name);
	nlines++;
	logxextin=1.47*(EXT(6300.)-EXT(6563.));
	logyextin=1.47*(EXT(5007.)-EXT(4861.));
      }
      
    }
    break;
  case 3:
    xmin=-1.5;xmax=1.8;
    ymin=-1.0;ymax=1.4;
    strcpy(xlabel,"log([OII]3727/[OII]5007)");
    strcpy(ylabel,"log([OIII]5007/H\\gb)");
    strcpy(xlabel,"[OII]3727/[OIII]5007");
    strcpy(ylabel,"[OIII]5007/H\\gb");
    if(DEBUG) printf(" Estoy en uno\n");
    nlines=0;
    for(j=0;j<ngal;j++) {
      if((G[j]).ifOIII && (G[j]).ifHb && (G[j]).ifOIII && (G[j]).ifOII && (G[j]).fHb!=0 && (G[j]).fOIII!=0) {
	G[j].x=(G[j]).fOII/(G[j]).fOIII;
	G[j].y=(G[j]).fOIII/(G[j]).fHb;
	x[nlines]=G[j].x;
	y[nlines]=G[j].y;
	if(errflag) {
	  errx[nlines]=sqrt((G[j]).efOII/(G[j]).fOIII*(G[j]).efOII/(G[j]).fOIII+(G[j]).fOII/((G[j]).fOIII*(G[j]).fOIII)*(G[j]).fOII/((G[j]).fOIII*(G[j]).fOIII)*(G[j]).efOIII*(G[j]).efOIII);
	  erry[nlines]=sqrt((G[j]).efOIII/(G[j]).fHb*(G[j]).efOIII/(G[j]).fHb+(G[j]).fOIII/((G[j]).fHb*(G[j]).fHb)*(G[j]).fOIII/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb);
	}
	id[nlines]=j;	  ids[nlines]=G[j].dataset;
	strcpy(labels+nlines*51,G[j].name);
	nlines++;
	logxextin=1.47*(EXT(3727.)-EXT(5007.));
	logyextin=1.47*(EXT(5007.)-EXT(4861.));
      }
      
    }
    break;
  case 4:
    xmin=-1.0;xmax=1.8;
    ymin=-1.0;ymax=1.4;
    strcpy(xlabel,"log([OII]3727/H\\gb)");
    strcpy(ylabel,"log([OIII]5007/H\\gb)");
    strcpy(xlabel,"[OII]3727/H\\gb");
    strcpy(ylabel,"[OIII]5007/H\\gb");
    if(DEBUG) printf(" Estoy en uno\n");
    nlines=0;
    for(j=0;j<ngal;j++) {
      if((G[j]).ifOII && (G[j]).ifHb && (G[j]).ifOIII && (G[j]).fHb!=0) {
	G[j].x=(G[j]).fOII/(G[j]).fHb;
	G[j].y=(G[j]).fOIII/(G[j]).fHb;
	x[nlines]=G[j].x;
	y[nlines]=G[j].y;
	if(errflag) {
	  errx[nlines]=sqrt((G[j]).efOII/(G[j]).fHb*(G[j]).efOII/(G[j]).fHb+(G[j]).fOII/((G[j]).fHb*(G[j]).fHb)*(G[j]).fOII/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb);
	  erry[nlines]=sqrt((G[j]).efOIII/(G[j]).fHb*(G[j]).efOIII/(G[j]).fHb+(G[j]).fOIII/((G[j]).fHb*(G[j]).fHb)*(G[j]).fOIII/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb);
	  if(DEBUG) printf(" %s fOII %g +/- %g %g   efHb %g +/- %g %g     x %g +/- %g   %g %g\n",G[j].name,(G[j]).fOII,(G[j]).efOII,(G[j]).fOII/(G[j]).efOII,(G[j]).fHb,(G[j]).efHb,(G[j]).fHb/(G[j]).efHb,x[nlines],errx[nlines],(G[j]).efOII/(G[j]).fHb*(G[j]).efOII/(G[j]).fHb,(G[j]).fOII/((G[j]).fHb*(G[j]).fHb)*(G[j]).fOII/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb);
	}
	id[nlines]=j;	  ids[nlines]=G[j].dataset;
	strcpy(labels+nlines*51,G[j].name);
	nlines++;
	logxextin=1.47*(EXT(3727.)-EXT(4861.));
	logyextin=1.47*(EXT(5007.)-EXT(4861.));
      }
      
    }
    break;
  case 5:
    xmin=-1.8;xmax=1.5;
    ymin=-1.0;ymax=1.4;
    strcpy(xlabel,"log(([SII]6717+[SII]6731)/H\\gb)");
    strcpy(ylabel,"log([OIII]5007/H\\gb)");
    strcpy(xlabel,"([SII]6717+[SII]6731)/H\\gb");
    strcpy(ylabel,"[OIII]5007/H\\gb");
    if(DEBUG) printf(" Estoy en uno\n");
    nlines=0;
    for(j=0;j<ngal;j++) {
      if((G[j]).ifOIII && (G[j]).ifHb && (G[j]).ifSII_1 && (G[j]).ifSII_2 && (G[j]).fHb!=0) {
	G[j].x=((G[j]).fSII_1+(G[j]).fSII_2)/(G[j]).fHb;
	G[j].y=(G[j]).fOIII/(G[j]).fHb;
	x[nlines]=G[j].x;
	y[nlines]=G[j].y;
	if(errflag) {
	  errx[nlines]=sqrt((G[j]).efSII_2/(G[j]).fHb*(G[j]).efSII_2/(G[j]).fHb+(G[j]).efSII_1/(G[j]).fHb*(G[j]).efSII_1/(G[j]).fHb+((G[j]).fSII_1+(G[j]).fSII_2)/((G[j]).fHb*(G[j]).fHb)*((G[j]).fSII_1+(G[j]).fSII_2)/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb);
	  erry[nlines]=sqrt((G[j]).efOIII/(G[j]).fHb*(G[j]).efOIII/(G[j]).fHb+(G[j]).fOIII/((G[j]).fHb*(G[j]).fHb)*(G[j]).fOIII/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb);
	  if(DEBUG) printf(" %s fSII %g +/- %g %g   efHb %g +/- %g %g     x %g +/- %g   %g %g\n",G[j].name,(G[j]).fSII_1,(G[j]).efSII_1,(G[j]).fSII_1/(G[j]).efSII_1,(G[j]).fHb,(G[j]).efHb,(G[j]).fHb/(G[j]).efHb,x[nlines],errx[nlines],(G[j]).efSII_1/(G[j]).fHb*(G[j]).efSII_1/(G[j]).fHb,(G[j]).fSII_1/((G[j]).fHb*(G[j]).fHb)*(G[j]).fSII_1/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb);
	}
	id[nlines]=j;	  ids[nlines]=G[j].dataset;
	strcpy(labels+nlines*51,G[j].name);
	nlines++;
	logxextin=1.47*((EXT(6717.)+EXT(6731.))/2.-EXT(4861.));  /* Esto es aproximado...*/
	logyextin=1.47*(EXT(5007.)-EXT(4861.));
      }
      
    }
    break;

  case 6:
    xmin=-0.7;xmax=1.5;
    ymin=-2.3;ymax=0.4;
    strcpy(xlabel,"log(([OII]3727+[OIII]5007)/H\\gb)");
    strcpy(ylabel,"log([NII]6583/[OII]3727)");
    strcpy(xlabel,"([OII]3727+[OIII]5007)/H\\gb");
    strcpy(ylabel,"[NII]6583/[OII]3727");
    if(DEBUG) printf(" Estoy en uno\n");
    nlines=0;
    for(j=0;j<ngal;j++) {
      if((G[j]).ifOIII && (G[j]).ifHb && (G[j]).ifOII && (G[j]).ifNII && (G[j]).fOII!=0 && (G[j]).fHb!=0) {
	G[j].x=((G[j]).fOIII+(G[j]).fOII)/(G[j]).fHb;
	G[j].y=(G[j]).fNII/(G[j]).fOII;
	x[nlines]=G[j].x;
	y[nlines]=G[j].y;
	if(errflag) {
	  errx[nlines]=sqrt((G[j]).efOIII/(G[j]).fHb*(G[j]).efOIII/(G[j]).fHb+(G[j]).efOII/(G[j]).fHb*(G[j]).efOII/(G[j]).fHb+((G[j]).fOII+(G[j]).fOIII)/((G[j]).fHb*(G[j]).fHb)*((G[j]).fOII+(G[j]).fOIII)/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb);
	  erry[nlines]=sqrt((G[j]).efNII/(G[j]).fOII*(G[j]).efNII/(G[j]).fOII+(G[j]).fNII/((G[j]).fOII*(G[j]).fOII)*(G[j]).fNII/((G[j]).fOII*(G[j]).fOII)*(G[j]).efOII*(G[j]).efOII);
	}
	id[nlines]=j;	  ids[nlines]=G[j].dataset;
	strcpy(labels+nlines*51,G[j].name);
	nlines++;
	logxextin=1.47*((EXT(3727.)+EXT(5007.))/2.-EXT(4861.)); /*Aproximado */
	logyextin=1.47*(EXT(6583.)-EXT(3727.));
      }
      
    }
    break;




  case 7:
    xmin=0;xmax=3;
    ymin=-4;ymax=0;
    strcpy(xlabel,"log(EW(H\\gb))");
    strcpy(ylabel,"log([OI]6300/H\\gb)");
    strcpy(xlabel,"EW(H\\gb)");
    strcpy(ylabel,"[OI]6300/H\\gb");
    if(DEBUG) printf(" Estoy en uno\n");
    nlines=0;
    for(j=0;j<ngal;j++) {
      if((G[j]).ifOI && (G[j]).ifHb && (G[j]).iewHb && (G[j]).fHb!=0 ) {
	G[j].x=(G[j]).ewHb;
	G[j].y=(G[j]).fOI/(G[j]).fHb;
	x[nlines]=G[j].x;
	y[nlines]=G[j].y;
	if(errflag) {
	  errx[nlines]=(G[j]).eewHb;
	  erry[nlines]=sqrt((G[j]).efOI/(G[j]).fHb*(G[j]).efOI/(G[j]).fHb+(G[j]).fOI/((G[j]).fHb*(G[j]).fHb)*(G[j]).fOI/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb);
	}
	id[nlines]=j;	  ids[nlines]=G[j].dataset;
	strcpy(labels+nlines*51,G[j].name);
	nlines++;
	logxextin=0.;
	logyextin=1.47*(EXT(6300.)-EXT(4861.));
      }
      
    }
    break;
  case 8:
    xmin=0;xmax=3;
    ymin=-3;ymax=0;
    strcpy(xlabel,"log(EW(H\\gb))");
    strcpy(ylabel,"log([HeI]5876/H\\gb)");
    strcpy(xlabel,"EW(H\\gb)");
    strcpy(ylabel,"[HeI]5876/H\\gb");
    if(DEBUG) printf(" Estoy en uno\n");
    nlines=0;
    for(j=0;j<ngal;j++) {
      if((G[j]).ifHeI && (G[j]).ifHb && (G[j]).iewHb && (G[j]).fHb!=0 ) {
	G[j].x=(G[j]).ewHb;
	G[j].y=(G[j]).fHeI/(G[j]).fHb;
	x[nlines]=G[j].x;
	y[nlines]=G[j].y;
	if(errflag) {
	  errx[nlines]=(G[j]).eewHb;
	  erry[nlines]=sqrt((G[j]).efHeI/(G[j]).fHb*(G[j]).efHeI/(G[j]).fHb+(G[j]).fHeI/((G[j]).fHb*(G[j]).fHb)*(G[j]).fHeI/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb);
	}
	id[nlines]=j;	  ids[nlines]=G[j].dataset;
	strcpy(labels+nlines*51,G[j].name);
	nlines++;
	logxextin=0.;
	logyextin=1.47*(EXT(5876.)-EXT(4861.));
      }
      
    }
    break;
  case 9:
    xmin=0;xmax=3;
    ymin=-1.0;ymax=1.4;
    strcpy(xlabel,"log(EW(H\\gb))");
    strcpy(ylabel,"log([OIII]5007/H\\gb)");
    strcpy(xlabel,"EW(H\\gb)");
    strcpy(ylabel,"[OIII]5007/H\\gb");
    if(DEBUG) printf(" Estoy en uno\n");
    nlines=0;
    for(j=0;j<ngal;j++) {
      if((G[j]).ifOIII && (G[j]).ifHb && (G[j]).iewHb && (G[j]).fHb!=0 ) {
	G[j].x=(G[j]).ewHb;
	G[j].y=(G[j]).fOIII/(G[j]).fHb;
	x[nlines]=G[j].x;
	y[nlines]=G[j].y;
	if(errflag) {
	  errx[nlines]=(G[j]).eewHb;
	  erry[nlines]=sqrt((G[j]).efOIII/(G[j]).fHb*(G[j]).efOIII/(G[j]).fHb+(G[j]).fOIII/((G[j]).fHb*(G[j]).fHb)*(G[j]).fOIII/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb);
	}
	id[nlines]=j;	  ids[nlines]=G[j].dataset;
	strcpy(labels+nlines*51,G[j].name);
	nlines++;
	logxextin=0.;
	logyextin=1.47*(EXT(5007.)-EXT(4861.));
      }
      
    }
    break;
  case 10:
    xmin=0;xmax=3;
    ymin=-1.3;ymax=1.4;
    strcpy(xlabel,"log(EW(H\\gb))");
    strcpy(ylabel,"log([OII]3727/H\\gb)");
    strcpy(xlabel,"EW(H\\gb)");
    strcpy(ylabel,"[OII]3727/H\\gb");
    if(DEBUG) printf(" Estoy en uno\n");
    nlines=0;
    for(j=0;j<ngal;j++) {
      if((G[j]).ifOII && (G[j]).ifHb && (G[j]).iewHb && (G[j]).fHb!=0 ) {
	G[j].x=(G[j]).ewHb;
	G[j].y=(G[j]).fOII/(G[j]).fHb;
	x[nlines]=G[j].x;
	y[nlines]=G[j].y;
	if(errflag) {
	  errx[nlines]=(G[j]).eewHb;
	  erry[nlines]=sqrt((G[j]).efOII/(G[j]).fHb*(G[j]).efOII/(G[j]).fHb+(G[j]).fOII/((G[j]).fHb*(G[j]).fHb)*(G[j]).fOII/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb);
	}
	id[nlines]=j;	  ids[nlines]=G[j].dataset;
	strcpy(labels+nlines*51,G[j].name);
	nlines++;
	logxextin=0.;
	logyextin=1.47*(EXT(3727.)-EXT(4861.));
      }
      
    }
    break;
  case 11:
    xmin=0;xmax=3;
    ymin=-0.7;ymax=1.5;
    strcpy(xlabel,"log(EW(H\\gb))");
    strcpy(ylabel,"log(([OII]3727+[OIII]5007)/H\\gb)");
    strcpy(xlabel,"EW(H\\gb)");
    strcpy(ylabel,"([OII]3727+[OIII]5007)/H\\gb");
    if(DEBUG) printf(" Estoy en uno\n");
    nlines=0;
    for(j=0;j<ngal;j++) {
      if((G[j]).ifOII && (G[j]).ifOIII && (G[j]).ifHb && (G[j]).iewHb && (G[j]).fHb!=0 ) {
	G[j].x=(G[j]).ewHb;
	G[j].y=((G[j]).fOII+(G[j]).fOIII)/(G[j]).fHb;
	x[nlines]=G[j].x;
	y[nlines]=G[j].y;
	if(errflag) {
	  errx[nlines]=(G[j]).eewHb;
	  erry[nlines]=sqrt(((G[j]).efOII+(G[j]).efOIII)/(G[j]).fHb*((G[j]).efOII+(G[j]).efOIII)/(G[j]).fHb+((G[j]).fOII+(G[j]).efOIII)/((G[j]).fHb*(G[j]).fHb)*((G[j]).fOII+(G[j]).efOIII)/((G[j]).fHb*(G[j]).fHb)*(G[j]).efHb*(G[j]).efHb);
	}
	id[nlines]=j;	  ids[nlines]=G[j].dataset;
	strcpy(labels+nlines*51,G[j].name);
	nlines++;
	logxextin=0.;
	logyextin=1.47*((EXT(3727.)+EXT(5007.))/2.-EXT(4861.)); /*Aproximado */
/*  	logyextin=1.47*(EXT(3727.)+EXT(5007.)-2*EXT(4861.));  */
      }
      
    }
    break;
  case 12:
    xmin=0;xmax=3;
    ymin=-2.3;ymax=0.4;
    strcpy(xlabel,"log(EW(H\\gb))");
    strcpy(ylabel,"log([NII]6583/[OII]3727)");
    strcpy(xlabel,"EW(H\\gb)");
    strcpy(ylabel,"[NII]6583/[OII]3727");
    if(DEBUG) printf(" Estoy en uno\n");
    nlines=0;
    for(j=0;j<ngal;j++) {
      if((G[j]).ifOII && (G[j]).ifNII && (G[j]).iewHb && (G[j]).fOII!=0 ) {
	G[j].x=(G[j]).ewHb;
	G[j].y=(G[j]).fNII/(G[j]).fOII;
	x[nlines]=G[j].x;
	y[nlines]=G[j].y;
	if(errflag) {
	  errx[nlines]=(G[j]).eewHb;
	  erry[nlines]=sqrt((G[j]).efNII/(G[j]).fOII*(G[j]).efNII/(G[j]).fOII+(G[j]).fNII/((G[j]).fOII*(G[j]).fOII)*(G[j]).fNII/((G[j]).fOII*(G[j]).fOII)*(G[j]).efOII*(G[j]).efOII);
	}
	id[nlines]=j;	  ids[nlines]=G[j].dataset;
	strcpy(labels+nlines*51,G[j].name);
	nlines++;
	logxextin=0.;
	logyextin=1.47*(EXT(6583.)-EXT(3727.));
      }
      
    }
    break;
    
    if(DEBUG) printf(" Por aquipa\n");
    nlines=0;
    
  }
  
  /*Dibujo el diagrama */
  
  if(xmin!=xmax) {
    
    
    cpgpage();
    cpgsci(1);
    cpgsls(1);
    
    cpgsch(1.5); 
    cpgvstd();  
    cpglab(xlabel,ylabel,""); 
    cpgsch(2.);
    cpgvstd();
    cpgsch(1.3);
    cpgswin(xmin,xmax,ymin,ymax);
/*     cpglab(xlabel,ylabel,""); */
    if(DEBUG) printf(" ante box\n");
    cpgbox("BCTNSL",0,0,"VBCTNSL",0,0); 
    if(DEBUG) printf(" Antes pintar %d\n",nlines);
    printf(" %d objects plotted\n",nlines);
    for(j=0;j<nlines;j++) {
      if(DEBUG2) printf("digujo %d\n",j);
      if(DEBUG) printf("%d  x %f  %f y %f  %f name %s\n",id[j],x[j],log10(x[j]),y[j],log10(y[j]),labels+j*51);
      cpgsci(color[ids[j]]);
      cpgsch(1.3);
      cpgpt1(log10(x[j]),log10(y[j]),symbol[ids[j]]); 
      cpgsch(0.8);
      if(labflag[ids[j]]) cpgtext(log10(x[j]),log10(y[j]),labels+j*51);
      if(errflag) {
	dum1=(log10(x[j]-errx[j]));dum2=(log10(x[j]+errx[j]));
	dum3=log10(y[j]);
	if(x[j]-errx[j]<0) {
	  dum1=log10(x[j])-errx[j]/x[j]/log(10);
	  dum2=log10(x[j])+errx[j]/x[j]/log(10);
	}      
	if(DEBUG) printf(" errx %f interva x %f %f\n",errx[j],dum1,dum2);
	cpgerrx(1,&dum1,&dum2,&dum3,0.25);
	dum1=(log10(y[j]-erry[j]));dum2=(log10(y[j]+erry[j]));
	dum3=log10(x[j]);
	if(y[j]-erry[j]<0) {
	  dum1=log10(y[j])-erry[j]/y[j]/log(10.);
	  dum2=log10(y[j])+erry[j]/y[j]/log(10.);
	}       
	if(DEBUG) printf(" erry %f interva y %f %f\n",erry[j],dum1,dum2);
	cpgerry(1,&dum3,&dum1,&dum2,0.25);
      }
    }
    cpgsch(1.);
    cpgsci(1);
    if(DEBUG) printf(" Extin x %f y %f\n",logxextin,logyextin);
    if(xmin+(xmax-xmin)*0.20+logxextin>xmin && ymin+(ymax-ymin)*0.17+logyextin > ymin) {
      cpgptxt(xmin+(xmax-xmin)*0.20,ymin+(ymax-ymin)*0.15,0.0,0.0,"E(B-V)=1");
      cpgsch(0.8);
      cpgarro(xmin+(xmax-xmin)*0.18,ymin+(ymax-ymin)*0.15,xmin+(xmax-xmin)*0.18+logxextin,ymin+(ymax-ymin)*0.15+logyextin);
    }
    else if(xmin+(xmax-xmin)*0.20+logxextin/2>xmin && ymin+(ymax-ymin)*0.17+logyextin/2 > ymin){
      cpgptxt(xmin+(xmax-xmin)*0.20,ymin+(ymax-ymin)*0.15,0.0,0.0,"E(B-V)=0.5");
      cpgsch(0.8);
      cpgarro(xmin+(xmax-xmin)*0.18,ymin+(ymax-ymin)*0.15,xmin+(xmax-xmin)*0.18+logxextin/2.,ymin+(ymax-ymin)*0.15+logyextin/2.);
    }
    else {
      cpgptxt(xmin+(xmax-xmin)*0.20,ymin+(ymax-ymin)*0.15,0.0,0.0,"E(B-V)=0.25");
      cpgsch(0.8);
      cpgarro(xmin+(xmax-xmin)*0.18,ymin+(ymax-ymin)*0.15,xmin+(xmax-xmin)*0.18+logxextin/4.,ymin+(ymax-ymin)*0.15+logyextin/4.);
    }
    
    cpgsci(1);
  }
  
  
  free(x);free(y);free(errx);free(erry);free(id);free(ids);free(labels);  
  
}


void PlotCMHK(int nplot) {

  int iz,ia;
/*   float x,y; */

  if(DEBUG) printf(" entra mod\n");

  switch(nplot) {
  case 0:
    break;
  case 1:
    for(iz=0;iz<CMHK.nmet;iz++) {
      cpgsls(iz+1); 
      cpgmove(CMHK.Ginst[iz][0].fNII-CMHK.Ginst[iz][0].fHa,CMHK.Ginst[iz][0].fOIII-CMHK.Ginst[iz][0].fHb);
      for(ia=0;ia<CMHK.nage;ia++) {
	if(DEBUG) printf("%d z %f age %f  %f %f MOD %f %f \n",iz,CMHK.z[iz],CMHK.age[ia],CMHK.Ginst[iz][ia].fNII,CMHK.Ginst[iz][ia].fHa,CMHK.Ginst[iz][ia].fNII-CMHK.Ginst[iz][ia].fHa,CMHK.Ginst[iz][ia].fOIII-CMHK.Ginst[iz][ia].fHb);
 	cpgdraw(CMHK.Ginst[iz][ia].fNII-CMHK.Ginst[iz][ia].fHa,CMHK.Ginst[iz][ia].fOIII-CMHK.Ginst[iz][ia].fHb); 
      }
    }
    break;
  case 2:
    for(iz=0;iz<CMHK.nmet;iz++) {
      cpgsls(iz+1); 
      cpgmove(CMHK.Ginst[iz][0].fOI-CMHK.Ginst[iz][0].fHa,CMHK.Ginst[iz][0].fOIII-CMHK.Ginst[iz][0].fHb);
      for(ia=0;ia<CMHK.nage;ia++) {
 	if(DEBUG) printf(" %f %f \n",CMHK.Ginst[iz][ia].fOI-CMHK.Ginst[iz][ia].fHa,CMHK.Ginst[iz][ia].fOIII-CMHK.Ginst[iz][ia].fHb); 
 	cpgdraw(CMHK.Ginst[iz][ia].fOI-CMHK.Ginst[iz][ia].fHa,CMHK.Ginst[iz][ia].fOIII-CMHK.Ginst[iz][ia].fHb); 
      }
    }
    break;
  case 3:
    for(iz=0;iz<CMHK.nmet;iz++) {
      cpgsls(iz+1); 
      cpgmove(CMHK.Ginst[iz][0].fOII-CMHK.Ginst[iz][0].fOIII,CMHK.Ginst[iz][0].fOIII-CMHK.Ginst[iz][0].fHb);
      for(ia=0;ia<CMHK.nage;ia++) {
 	cpgdraw(CMHK.Ginst[iz][ia].fOII-CMHK.Ginst[iz][ia].fOIII,CMHK.Ginst[iz][ia].fOIII-CMHK.Ginst[iz][ia].fHb); 
      }
    }
    break;
  case 4:
    for(iz=0;iz<CMHK.nmet;iz++) {
      cpgsls(iz+1); 
      cpgmove(CMHK.Ginst[iz][0].fOII-CMHK.Ginst[iz][0].fHb,CMHK.Ginst[iz][0].fOIII-CMHK.Ginst[iz][0].fHb);
      for(ia=0;ia<CMHK.nage;ia++) {
 	cpgdraw(CMHK.Ginst[iz][ia].fOII-CMHK.Ginst[iz][ia].fHb,CMHK.Ginst[iz][ia].fOIII-CMHK.Ginst[iz][ia].fHb); 
      }
    }
    break;
  case 5:
/*     for(iz=0;iz<CMHK.nmet;iz++) { */
/*       cpgsls(iz+1);  */
/*       cpgmove(CMHK.Ginst[iz][0].fNII-CMHK.Ginst[iz][0].fHa,CMHK.Ginst[iz][0].fOIII-CMHK.Ginst[iz][0].fHb); */
/*       for(ia=0;ia<CMHK.nage;ia++) { */
/* 	if(DEBUG) printf("%d z %f age %f  %f %f MOD %f %f \n",iz,CMHK.z[iz],CMHK.age[ia],CMHK.Ginst[iz][ia].fNII,CMHK.Ginst[iz][ia].fHa,CMHK.Ginst[iz][ia].fNII-CMHK.Ginst[iz][ia].fHa,CMHK.Ginst[iz][ia].fOIII-CMHK.Ginst[iz][ia].fHb); */
/*  	cpgdraw(CMHK.Ginst[iz][ia].fNII-CMHK.Ginst[iz][ia].fHa,CMHK.Ginst[iz][ia].fOIII-CMHK.Ginst[iz][ia].fHb);  */
/*       } */
/*     } */
    break;
  case 6:
    for(iz=0;iz<CMHK.nmet;iz++) {
      cpgsls(iz+1); 
      cpgmove(log10(pow(10.,CMHK.Ginst[iz][0].fOIII)+pow(10.,CMHK.Ginst[iz][0].fOII))-CMHK.Ginst[iz][0].fHb,CMHK.Ginst[iz][0].fNII-CMHK.Ginst[iz][0].fOII);
      for(ia=0;ia<CMHK.nage;ia++) {
	if(DEBUG) printf("%d  OII %f OIII %f Hb%f  coc %f\n",iz,CMHK.Ginst[iz][ia].fOII,CMHK.Ginst[iz][ia].fOIII,CMHK.Ginst[iz][ia].fHb,log10(pow(10.,CMHK.Ginst[iz][ia].fOIII)+pow(10.,CMHK.Ginst[iz][ia].fOII))-CMHK.Ginst[iz][ia].fHb);
 	cpgdraw(log10(pow(10.,CMHK.Ginst[iz][ia].fOIII)+pow(10.,CMHK.Ginst[iz][ia].fOII))-CMHK.Ginst[iz][ia].fHb,CMHK.Ginst[iz][ia].fNII-CMHK.Ginst[iz][ia].fOII); 
      }
    }
    break;
  case 7:
    for(iz=0;iz<CMHK.nmet;iz++) {
      cpgsls(iz+1);
      cpgmove(CMHK.Ginst[iz][0].ewHb,CMHK.Ginst[iz][0].fOI-CMHK.Ginst[iz][0].fHb);
      for(ia=0;ia<CMHK.nage;ia++) {
	cpgdraw(log10(CMHK.Ginst[iz][ia].ewHb),CMHK.Ginst[iz][ia].fOI-CMHK.Ginst[iz][ia].fHb);
      }
    }
    break;
  case 8:
/*     for(iz=0;iz<CMHK.nmet;iz++) { */
/*       cpgsls(iz+1); */
/*       cpgmove(CMHK.Ginst[iz][0].ewHb,CMHK.Ginst[iz][0].fOIII-CMHK.Ginst[iz][0].fHb); */
/*       for(ia=0;ia<CMHK.nage;ia++) { */
/* 	if(DEBUG) printf("%d z %f age %f  %f %f MOD %f %f \n",iz,CMHK.z[iz],CMHK.age[ia],CMHK.Ginst[iz][ia].fOIII,CMHK.Ginst[iz][ia].fHb,log10(CMHK.Ginst[iz][ia].ewHb),CMHK.Ginst[iz][ia].fOIII-CMHK.Ginst[iz][ia].fHb);	 */
/* 	cpgdraw(log10(CMHK.Ginst[iz][ia].ewHb),CMHK.Ginst[iz][ia].fOIII-CMHK.Ginst[iz][ia].fHb); */
/*       } */
/*     } */
    break;
  case 9:
    for(iz=0;iz<CMHK.nmet;iz++) {
      cpgsls(iz+1);
      cpgmove(CMHK.Ginst[iz][0].ewHb,CMHK.Ginst[iz][0].fOIII-CMHK.Ginst[iz][0].fHb);
      for(ia=0;ia<CMHK.nage;ia++) {
	if(DEBUG) printf("%d z %f age %f  %f %f MOD %f %f \n",iz,CMHK.z[iz],CMHK.age[ia],CMHK.Ginst[iz][ia].fOIII,CMHK.Ginst[iz][ia].fHb,log10(CMHK.Ginst[iz][ia].ewHb),CMHK.Ginst[iz][ia].fOIII-CMHK.Ginst[iz][ia].fHb);	
	cpgdraw(log10(CMHK.Ginst[iz][ia].ewHb),CMHK.Ginst[iz][ia].fOIII-CMHK.Ginst[iz][ia].fHb);
      }
    }
    break;
  case 10:
    for(iz=0;iz<CMHK.nmet;iz++) {
      cpgsls(iz+1);
      cpgmove(CMHK.Ginst[iz][0].ewHb,CMHK.Ginst[iz][0].fOII-CMHK.Ginst[iz][0].fHb);
      for(ia=0;ia<CMHK.nage;ia++) {
	cpgdraw(log10(CMHK.Ginst[iz][ia].ewHb),CMHK.Ginst[iz][ia].fOII-CMHK.Ginst[iz][ia].fHb);
      }
    }
    break;
  case 11:
    for(iz=0;iz<CMHK.nmet;iz++) {
      cpgsls(iz+1);
      cpgmove(CMHK.Ginst[iz][0].ewHb,log10(pow(10.,CMHK.Ginst[iz][0].fOIII)+pow(10.,CMHK.Ginst[iz][0].fOII))-CMHK.Ginst[iz][0].fHb);
      for(ia=0;ia<CMHK.nage;ia++) {
	cpgdraw(log10(CMHK.Ginst[iz][ia].ewHb),log10(pow(10.,CMHK.Ginst[iz][ia].fOIII)+pow(10.,CMHK.Ginst[iz][ia].fOII))-CMHK.Ginst[iz][ia].fHb);
      }
    }
    break;
  case 12:
    for(iz=0;iz<CMHK.nmet;iz++) {
      cpgsls(iz+1);
      cpgmove(CMHK.Ginst[iz][0].ewHb,CMHK.Ginst[iz][0].fNII-CMHK.Ginst[iz][0].fOII);
      for(ia=0;ia<CMHK.nage;ia++) {
	if(DEBUG) printf(" Z %f AGE %f NII %f OII %f x %f y %f \n",CMHK.z[iz],CMHK.age[ia],CMHK.Ginst[iz][ia].fNII,CMHK.Ginst[iz][ia].fOII,log10(CMHK.Ginst[iz][ia].ewHb),CMHK.Ginst[iz][ia].fNII-CMHK.Ginst[iz][ia].fOII);
	cpgdraw(log10(CMHK.Ginst[iz][ia].ewHb),CMHK.Ginst[iz][ia].fNII-CMHK.Ginst[iz][ia].fOII);
      }
    }
    break;
  
  }




}




void SaveParam(void) {

  FILE *fp;
  int nc=0,nt;
  char ch51[51];
  char opt='n';
  int i;
  char key[9]="";
  printf(" Do you want to save parameters in a file?: ");
  opt=readc(opt);
  if(opt=='y') {
    printf("Name of parameter file: ");
    reads(parfilename,parfilename);
    if((fp=fopen(parfilename,"w")) ==NULL) {
      printf("ERROR: Can't open file\n");
      return;
    }
    fprintf(fp,"COMMENT  Parameter file for StarBurst                                          \n");
    fprintf(fp,"COMMENT  If there are no errors, set those error columns to 0                  \n");
    nc+=2;
    for(i=0;i<nfiles;i++) {
      sprintf(key,"DATAF_%02d",i+1);
      sprintf(ch51,"'%s'",D[i].datafile);
      fprintf(fp,"%-9s= %-51.51s / F. /w flux data\n",key,ch51);
      sprintf(key,"LABEL_%02d",i+1);
      sprintf(ch51,"'%s'",D[i].label);
      fprintf(fp,"%-9s= %-51.51s / F. /w flux data\n",key,ch51);
    }
    sprintf(ch51,"'%s'",D[i].label);
    fprintf(fp,"EXTLAW  = %-51.51s / Extinction file\n",ch51);
    sprintf(ch51,"'%s'",CMHKdir);
    fprintf(fp,"CMHKMOD = %-51.51s / Dir CMHK models\n",ch51);
    fprintf(fp,"OBJCOL  =%21d / Column with object name in file DATAFILE      \n",objcol      );
    fprintf(fp,"FHALFA  =%21d / Column with Halfa flux          in DATAFILE   \n",fHacol      );
    fprintf(fp,"EFHALFA =%21d / Column with Halfa flux error    in DATAFILE   \n",fHaerrcol    );
    fprintf(fp,"EWHALFA =%21d / Column with Halfa EW            in DATAFILE   \n",ewHacol      );
    fprintf(fp,"EEWHALFA=%21d / Column with Halfa EW   error    in DATAFILE   \n",ewHaerrcol    );
    fprintf(fp,"FHBETA  =%21d / Column with Hbeta flux          in DATAFILE   \n",fHbcol      );
    fprintf(fp,"EFHBETA =%21d / Column with Hbeta flux error    in DATAFILE   \n",fHberrcol    );
    fprintf(fp,"EWHBETA =%21d / Column with Hbeta EW            in DATAFILE   \n",ewHbcol      );
    fprintf(fp,"EEWHBETA=%21d / Column with Hbeta EW   error    in DATAFILE   \n",ewHberrcol    );
    fprintf(fp,"FOIII   =%21d / Column with OIII5007 flux       in DATAFILE   \n",fOIIIcol      );
    fprintf(fp,"EFOIII  =%21d / Column with OIII5008 flux error in DATAFILE   \n",fOIIIerrcol    );
    fprintf(fp,"FOII    =%21d / Column with OII3727  flux       in DATAFILE   \n",fOIIcol      );
    fprintf(fp,"EFOII   =%21d / Column with OII3727  flux error in DATAFILE   \n",fOIIerrcol    );
    fprintf(fp,"FNII    =%21d / Column with NII6583  flux       in DATAFILE   \n",fNIIcol      );
    fprintf(fp,"EFNII   =%21d / Column with NII6583  flux error in DATAFILE   \n",fNIIerrcol    );
    fprintf(fp,"FOI     =%21d / Column with OI6300   flux       in DATAFILE   \n",fOIcol      );
    fprintf(fp,"EFOI    =%21d / Column with OI6300   flux error in DATAFILE   \n",fOIerrcol    );
    fprintf(fp,"FHEI    =%21d / Column with HeI5876  flux       in DATAFILE   \n",fHeIcol      );
    fprintf(fp,"EFHEI   =%21d / Column with HeI5876  flux error in DATAFILE   \n",fHeIerrcol    );
    fprintf(fp,"FSII_17 =%21d / Column with SII6717  flux       in DATAFILE   \n",fSII_1col      );
    fprintf(fp,"EFSII_17=%21d / Column with SII6717  flux error in DATAFILE   \n",fSII_1errcol    );
    fprintf(fp,"FSII_31 =%21d / Column with SII6731  flux       in DATAFILE   \n",fSII_2col      );
    fprintf(fp,"EFSII_31=%21d / Column with SII6731  flux error in DATAFILE   \n",fSII_2errcol    );
    sprintf(ch51,"'%s'","?" );
    fprintf(fp,"DEVICE  = %-51.51s / PGPLOT device  \n",ch51);
    nc+=26+nfiles*2;
    fprintf(fp,"COMMENT  Blank lines to assure the last block of 2880 caracters is complete    \n");
    fprintf(fp,"COMMENT                                                                        \n");
    nc+=2;
    for(nt=nc;nt<(1+(int)(nc/36))*36-1;nt++) {
      fprintf(fp,"COMMENT                                                                        \n");
    }
    fprintf(fp,"END                                                                            \n");
    fclose(fp);
  }
}



void ReadCat(void) {
  
  int i;
  int j;
/*   char name[2000]; */
/*   int iname[2000]; */
  int nwhole=0;
  int igal;
  for(i=0;i<nfiles;i++) {
    nwhole+=FileNLin(D[i].datafile);
  }
  if((G=malloc(nwhole*sizeof(struct gal)))==NULL) printf("I cannot dimension G         of %d elements \n",nwhole);

  
  igal=0;

  for(i=0;i<nfiles;i++) {
    printf(" Reading file %s ...\n",D[i].datafile);
    D[i].ncat=FileNLin(D[i].datafile);
    if(DEBUG) printf(" numer lin %d\n",D[i].ncat);
    if((D[i].name    =malloc(D[i].ncat*51*sizeof(char)))==NULL) printf("I cannot dimension name         of %d elements \n",D[i].ncat*51);
    if((D[i].fOIII   =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].fOIII     of %d elements \n",D[i].ncat);
    if((D[i].efOIII  =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].efOIII    of %d elements \n",D[i].ncat);
    if((D[i].fOII    =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].fOII      of %d elements \n",D[i].ncat);
    if((D[i].efOII   =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].efOII     of %d elements \n",D[i].ncat);
    if((D[i].fHa     =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].fHa       of %d elements \n",D[i].ncat);
    if((D[i].efHa    =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].efHa      of %d elements \n",D[i].ncat);
    if((D[i].fHb     =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].fHb       of %d elements \n",D[i].ncat);
    if((D[i].efHb    =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].efHb      of %d elements \n",D[i].ncat);
    if((D[i].fNII    =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].fNII      of %d elements \n",D[i].ncat);
    if((D[i].efNII   =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].efNII     of %d elements \n",D[i].ncat);
    if((D[i].fOI     =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].fOI       of %d elements \n",D[i].ncat);
    if((D[i].efOI    =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].efOI      of %d elements \n",D[i].ncat);
    if((D[i].fHeI    =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].fHeI      of %d elements \n",D[i].ncat);
    if((D[i].efHeI   =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].efHeI     of %d elements \n",D[i].ncat);
    if((D[i].fSII_1  =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].fSII_1    of %d elements \n",D[i].ncat);
    if((D[i].efSII_1 =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].efSII_1   of %d elements \n",D[i].ncat);
    if((D[i].fSII_2  =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].fSII_2    of %d elements \n",D[i].ncat);
    if((D[i].efSII_2 =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].efSII_2   of %d elements \n",D[i].ncat);
    if((D[i].ewHa    =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].ewHa      of %d elements \n",D[i].ncat);
    if((D[i].eewHa   =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].eewHa     of %d elements \n",D[i].ncat);
    if((D[i].ewHb    =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].ewHb      of %d elements \n",D[i].ncat);
    if((D[i].eewHb   =malloc(D[i].ncat*sizeof(float)))==NULL) printf("I cannot dimension D[i].eewHb     of %d elements \n",D[i].ncat);
    
    if((D[i].iname    =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension name    of %d elements \n",D[i].ncat);
    if((D[i].ifOIII   =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].iefOIII  =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].ifOII    =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].iefOII   =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].ifHa     =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].iefHa    =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].ifHb     =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].iefHb    =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].ifNII    =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].iefNII   =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].ifOI     =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].iefOI    =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].ifHeI    =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].iefHeI   =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].ifSII_1  =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].iefSII_1 =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].ifSII_2  =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].iefSII_2 =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].iewHa    =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].ieewHa   =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].iewHb    =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    if((D[i].ieewHb   =malloc(D[i].ncat*sizeof(int)))==NULL) printf("I cannot dimension ewo     of %d elements \n",D[i].ncat);
    
    for(j=0;j<D[i].ncat;j++) {
      (D[i]).iname[j]    =1;(D[i]).ifOIII[j]   =1;(D[i]).iefOIII[j]  =1;(D[i]).ifOII[j]    =1;(D[i]).iefOII [j]  =1;(D[i]).ifHa[j]     =1;(D[i]).iefHa[j]    =1;(D[i]).ifHb[j]     =1;(D[i]).iefHb[j]    =1;(D[i]).ifNII[j]    =1;(D[i]).iefNII[j]   =1;(D[i]).ifOI[j]     =1;(D[i]).iefOI[j]    =1;(D[i]).ifHeI[j]    =1;(D[i]).iefHeI[j]   =1;(D[i]).ifSII_1[j]  =1;(D[i]).iefSII_1[j] =1;(D[i]).ifSII_2[j]  =1;(D[i]).iefSII_2[j] =1;(D[i]).iewHa[j]    =1;(D[i]).ieewHa[j]   =1;(D[i]).iewHb[j]    =1;(D[i]).ieewHb[j]   =1; 
    }



    ReadCharcol((D[i]).datafile,objcol   ,(D[i]).name   ,(D[i]).iname,51,&((D[i]).ncat)); 
/*     ReadCharcol((D[i]).datafile,objcol   ,name   ,iname,51,&((D[i]).ncat));  */
    if(DEBUG) printf(" Names read\n");
    (D[i]).nobj=0; 
    for(j=0;j<(D[i]).ncat;j++) {
      if(DEBUG) printf(" Por el %d  ",j);
/*       if(DEBUG) printf(" name <%s> %d\n",name+j*51,iname[j]);     */
      if(DEBUG) printf(" name %s %d\n",(D[i]).name+j*51,((D[i])).iname[j]); 
     
      if((D[i]).iname[j]) ((D[i]).nobj)++;
    }

    if(DEBUG) printf(" Vale de aui\n");
    
    system("date");

/*     printf(" col III %d col II %d\n",fOIIIerrcol,fOIIerrcol); */

    ReadNumcol((D[i]).datafile,fOIIIcol       ,(D[i]).fOIII        ,(D[i]).ifOIII     ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fOIIIerrcol    ,(D[i]).efOIII       ,(D[i]).iefOIII    ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fOIIcol        ,(D[i]).fOII         ,(D[i]).ifOII      ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fOIIerrcol     ,(D[i]).efOII        ,(D[i]).iefOII     ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fHacol         ,(D[i]).fHa          ,(D[i]).ifHa       ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fHaerrcol      ,(D[i]).efHa         ,(D[i]).iefHa      ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fHbcol         ,(D[i]).fHb          ,(D[i]).ifHb       ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fHberrcol      ,(D[i]).efHb         ,(D[i]).iefHb      ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fNIIcol        ,(D[i]).fNII         ,(D[i]).ifNII      ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fNIIerrcol     ,(D[i]).efNII        ,(D[i]).iefNII     ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fOIcol         ,(D[i]).fOI          ,(D[i]).ifOI       ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fOIerrcol      ,(D[i]).efOI         ,(D[i]).iefOI      ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fHeIcol        ,(D[i]).fHeI         ,(D[i]).ifHeI      ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fHeIerrcol     ,(D[i]).efHeI        ,(D[i]).iefHeI     ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fSII_1col      ,(D[i]).fSII_1       ,(D[i]).ifSII_1    ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fSII_1errcol   ,(D[i]).efSII_1      ,(D[i]).iefSII_1   ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fSII_2col      ,(D[i]).fSII_2       ,(D[i]).ifSII_2    ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,fSII_2errcol   ,(D[i]).efSII_2      ,(D[i]).iefSII_2   ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,ewHacol        ,(D[i]).ewHa         ,(D[i]).iewHa      ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,ewHaerrcol     ,(D[i]).eewHa        ,(D[i]).ieewHa     ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,ewHbcol        ,(D[i]).ewHb         ,(D[i]).iewHb      ,&((D[i]).ncat));
    ReadNumcol((D[i]).datafile,ewHberrcol     ,(D[i]).eewHb        ,(D[i]).ieewHb     ,&((D[i]).ncat));
    system("date");

    for(j=0;j<(D[i]).ncat;j++) {
      if((D[i]).iname[j]) {
	G[igal].dataset=i;
	strcpy(G[igal].name,(D[i]).name+j*51);

	G[igal].fOIII   = (D[i]).fOIII[j]   ;G[igal].ifOIII   =  (D[i]).ifOIII[j]   ;
	G[igal].efOIII  = (D[i]).efOIII[j]  ;G[igal].iefOIII  =  (D[i]).iefOIII[j]  ;
	G[igal].fOII    = (D[i]).fOII[j]    ;G[igal].ifOII    =  (D[i]).ifOII[j]    ;
	G[igal].efOII   = (D[i]).efOII[j]   ;G[igal].iefOII   =  (D[i]).iefOII[j]   ;
	G[igal].fHa     = (D[i]).fHa[j]     ;G[igal].ifHa     =  (D[i]).ifHa[j]     ;
	G[igal].efHa    = (D[i]).efHa[j]    ;G[igal].iefHa    =  (D[i]).iefHa[j]    ;
	G[igal].fHb     = (D[i]).fHb[j]     ;G[igal].ifHb     =  (D[i]).ifHb[j]     ;
	G[igal].efHb    = (D[i]).efHb[j]    ;G[igal].iefHb    =  (D[i]).iefHb[j]    ;
	G[igal].fNII    = (D[i]).fNII[j]    ;G[igal].ifNII    =  (D[i]).ifNII[j]    ;
	G[igal].efNII   = (D[i]).efNII[j]   ;G[igal].iefNII   =  (D[i]).iefNII[j]   ;
	G[igal].fOI     = (D[i]).fOI[j]     ;G[igal].ifOI     =  (D[i]).ifOI[j]     ;
	G[igal].efOI    = (D[i]).efOI[j]    ;G[igal].iefOI    =  (D[i]).iefOI[j]    ;
	G[igal].fHeI    = (D[i]).fHeI[j]    ;G[igal].ifHeI    =  (D[i]).ifHeI[j]    ;
	G[igal].efHeI   = (D[i]).efHeI[j]   ;G[igal].iefHeI   =  (D[i]).iefHeI[j]   ;
	G[igal].fSII_1  = (D[i]).fSII_1[j]  ;G[igal].ifSII_1  =  (D[i]).ifSII_1[j]  ;
	G[igal].efSII_1 = (D[i]).efSII_1[j] ;G[igal].iefSII_1 =  (D[i]).iefSII_1[j] ;
	G[igal].fSII_2  = (D[i]).fSII_2[j]  ;G[igal].ifSII_2  =  (D[i]).ifSII_2[j]  ;
	G[igal].efSII_2 = (D[i]).efSII_2[j] ;G[igal].iefSII_2 =  (D[i]).iefSII_2[j] ;
	G[igal].ewHa    = (D[i]).ewHa[j]    ;G[igal].iewHa    =  (D[i]).iewHa[j]    ;
	G[igal].eewHa   = (D[i]).eewHa[j]   ;G[igal].ieewHa   =  (D[i]).ieewHa[j]   ;
	G[igal].ewHb    = (D[i]).ewHb[j]    ;G[igal].iewHb    =  (D[i]).iewHb[j]    ;
	G[igal].eewHb   = (D[i]).eewHb[j]   ;G[igal].ieewHb   =  (D[i]).ieewHb[j]   ;
	igal++;
      }
    }
    if(DEBUG) printf(" col %d %d \n",objcol,fHacol);
    
    printf(" Number of objects in file %s: %d \n",(D[i]).datafile,(D[i]).nobj);
        
    
  }

  if(DEBUG) printf(" Going out\n");

  ngal=igal;
}



void LoadParam_file(void) {

  int status=0;
  char comment[51];
  char key[9]="";
  int i;
  fitsfile *parfile;
  printf("Using parameter file <<%s>>\n",parfilename);
  if( ffopen2(&parfile,parfilename, READONLY, &status)) fits_report_error(stderr,status);

  nfiles=0;
  for(i=1;i<NMAXDATAF+1;i++) {
    sprintf(key,"DATAF_%02d",i);
/*     printf(" Looking key <%s>\n",key);  */
    if((ffgky(parfile,TSTRING,key,D[i-1].datafile,comment,&status))!=0 || !strcmp(D[-1].datafile,"NONE")) {
      break;
    }
    printf(" Data file %d : %s\n",i,D[i-1].datafile);
    status=0;
    sprintf(key,"LABEL_%02d",i);
    ffgky(parfile,TSTRING,key,D[i-1].label,comment,&status);
    fits_report_error(stderr,status);
    if(status) {
      printf(" Keyword  %s not found. Exiting\n",key);
      exit(1);
    }
    nfiles++;
    printf(" Total number of files: %d\n",nfiles);

  }
  if(nfiles==0) {
    printf(" No data files have been loaded. Exiting\n");
    exit(1);
  }
  status=0;
  
  fits_report_error(stderr,status);

  ffgky(parfile,TSTRING,"EXTLAW",extinfile,comment,&status); 
  fits_report_error(stderr,status); 
  ffgky(parfile,TSTRING,"CMHKMOD",CMHKdir,comment,&status); 
  fits_report_error(stderr,status); 
  ffgky(parfile,TINT,"OBJCOL",&objcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"FHALFA",&fHacol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"EFHALFA",&fHaerrcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"EWHALFA",&ewHacol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"EEWHALFA",&ewHaerrcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"FHBETA",&fHbcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"EFHBETA",&fHberrcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"EWHBETA",&ewHbcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"EEWHBETA",&ewHberrcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"FOIII",&fOIIIcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"EFOIII",&fOIIIerrcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"FOII",&fOIIcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"EFOII",&fOIIerrcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"FNII",&fNIIcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"EFNII",&fNIIerrcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"FOI",&fOIcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"EFOI",&fOIerrcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"FHEI",&fHeIcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"EFHEI",&fHeIerrcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"FSII_17",&fSII_1col,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"EFSII_17",&fSII_1errcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"FSII_31",&fSII_2col,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TINT,"EFSII_31",&fSII_2errcol,comment,&status);
  fits_report_error(stderr,status);
  ffgky(parfile,TSTRING,"DEVICE",pgdevice,comment,&status);
  fits_report_error(stderr,status);

  if(fHaerrcol!=0) errflag=1;

  if(strcmp(CMHKdir,"NONE")) CMHKflag=1;
  
  if(status!=0) {
    fits_report_error(stderr,status);
    printf("ERROR: Reading parameter file. Not enough keywords.\n");
    exit(1); 
  }

  fits_close_file(parfile,&status);
  
  cpgopen(pgdevice);
  cpgask(0);
  cpgsch(1.5);
  cpgsah(2,45.,1.0);
}


void LoadParam_kbd(void) {

  int i=0;
  printf(" Input file with data containing all information about lines: ");
  reads((D[i]).datafile,(D[i]).datafile);  
  cpgopen("?");
  cpgask(0);
  cpgsch(1.5);

}




void ReadExtlaw(void) {
  FILE *fp;
  
  int i;
  EX.n=FileNLin(extinfile);
  if((fp=fopen(extinfile,"r"))==NULL) {
    printf("Cannot open file %s\n",extinfile);
    exit(1);
  }
  
  if((EX.ldo=malloc(EX.n*sizeof(float)))== NULL) {
    printf("I can't dimension the vector EX.ldo of %d elements",EX.n);
    exit(1);
  }
  if((EX.y=malloc(EX.n*sizeof(float)))== NULL) {
    printf("I can't dimension the vector EX.y of %d elements",EX.n);
    exit(1);
  }
  for (i=0;i<EX.n;i++) {
    fscanf(fp," %f %f",EX.ldo+i,EX.y+i);
  }
}

float EXT(float ldo)
{
  return(Lagr2(EX.ldo,EX.y,EX.n,ldo));
}


void ReadCMHK(void) {

  char patternclolin[200];
  char patternnebem[200];
  glob_t pglob;
  int iz,ia;
  int idum;
  int err;
  char *snul;
  int nfile;
  
  float age[1000],fHa[1000],fHb[1000],fOI[1000],fOII[1000],fOIII[1000],fNII[1000],ewHb[1000];
  int   ilog[1000];

  snul=malloc(100*sizeof(char));
  for(ia=0;ia<1000;ia++) ilog[ia]=0;
  
  if(DEBUG) printf(" Entra\n");

  sprintf(patternclolin,"%s/clolin_*const*dat",CMHKdir);
  pglob.gl_offs = 0; 
  err=glob(patternclolin,GLOB_ERR,NULL,&pglob); 
  if(err) {
    printf(" There was an error 1 reading CMHK models. Exiting.\n");
    exit(1);
  }
  CMHK.nmet=pglob.gl_pathc;
  CMHK.nage=FileNLin(pglob.gl_pathv[0]);
  if((CMHK.z   =malloc(CMHK.nmet*sizeof(float    )))==NULL) printf("I cannot dimension CMHK.z    of %d elements \n",CMHK.nmet );
  if((CMHK.age =malloc(CMHK.nage*sizeof(float    )))==NULL) printf("I cannot dimension CMHK.age  of %d elements \n",CMHK.nage );
  for(iz=0;iz<CMHK.nmet;iz++) {
    if(DEBUG) printf("Readgin %s\n",pglob.gl_pathv[iz]);
    ReadNumcol(pglob.gl_pathv[iz],1 ,age,ilog,&nfile);
    ReadNumcol(pglob.gl_pathv[iz],2 ,fHb,ilog,&nfile);
    ReadNumcol(pglob.gl_pathv[iz],4 ,fHa,ilog,&nfile);
    ReadNumcol(pglob.gl_pathv[iz],10,fOI,ilog,&nfile);
    ReadNumcol(pglob.gl_pathv[iz],12,fOII,ilog,&nfile);
    ReadNumcol(pglob.gl_pathv[iz],14,fOIII,ilog,&nfile);
    ReadNumcol(pglob.gl_pathv[iz],16,fNII,ilog,&nfile);
    snul=strstr(pglob.gl_pathv[iz],"_z")+2;
    snul[3]='\0';
    CMHK.z[iz]=atof(snul)/100.;
    idum=0;
    for(ia=0;ia<nfile;ia++) {
      if(ilog[ia]) {
	CMHK.age[idum]=age[ia];
	(CMHK.Gconst[iz][idum]).fOIII=fOIII[ia];
	(CMHK.Gconst[iz][idum]).fOII=fOII[ia];
	(CMHK.Gconst[iz][idum]).fOI=fOI[ia];
	(CMHK.Gconst[iz][idum]).fHa=fHa[ia];
	(CMHK.Gconst[iz][idum]).fHb=fHb[ia];
	(CMHK.Gconst[iz][idum]).fNII=fNII[ia];
	idum++;
      }
    }
    CMHK.nage=idum;
  }
  globfree(&pglob);
  sprintf(patternnebem,"%s/nebem_*const*dat",CMHKdir);
  err=glob(patternnebem,GLOB_ERR,NULL,&pglob); 
  if(err) {
    printf(" There was an error 2 reading CMHK models. Exiting.\n");
    exit(1);
  }
  if(CMHK.nmet!=pglob.gl_pathc ) {
    printf(" There was an error 3 reading CMHK models. Exiting.\n");
    exit(1);
  }
  for(iz=0;iz<CMHK.nmet;iz++) {
    if(DEBUG) printf("Readgin %s\n",pglob.gl_pathv[iz]);
    ReadNumcol(pglob.gl_pathv[iz],4 ,ewHb,ilog,&nfile);
    idum=0;
    for(ia=0;ia<nfile;ia++) {
      if(ilog[ia]) {	
	(CMHK.Gconst[iz][idum]).ewHb=ewHb[ia];
	idum++;
      }
    }
    if(idum!=CMHK.nage) {
      printf(" There was an error 4 reading CMHK models. Exiting.\n");
      exit(1);
    }
  }  
  globfree(&pglob);
  sprintf(patternclolin,"%s/clolin_*inst*dat",CMHKdir);
  err=glob(patternclolin,GLOB_ERR,NULL,&pglob); 
  if(err) {
    printf(" There was an error 5 reading CMHK models. Exiting.\n");
    exit(1);
  }
  if(CMHK.nmet!=pglob.gl_pathc ) {
    printf(" There was an error 6 reading CMHK models. Exiting.\n");
    exit(1);
  }
  for(iz=0;iz<CMHK.nmet;iz++) {
    if(DEBUG) printf("Readgin %s\n",pglob.gl_pathv[iz]);
    ReadNumcol(pglob.gl_pathv[iz],1 ,age,ilog,&nfile);
    ReadNumcol(pglob.gl_pathv[iz],2 ,fHb,ilog,&nfile);
    ReadNumcol(pglob.gl_pathv[iz],4 ,fHa,ilog,&nfile);
    ReadNumcol(pglob.gl_pathv[iz],10,fOI,ilog,&nfile);
    ReadNumcol(pglob.gl_pathv[iz],12,fOII,ilog,&nfile);
    ReadNumcol(pglob.gl_pathv[iz],14,fOIII,ilog,&nfile);
    ReadNumcol(pglob.gl_pathv[iz],16,fNII,ilog,&nfile);
    idum=0;
    for(ia=0;ia<nfile;ia++) {
      if(ilog[ia]) {
	CMHK.age[idum]=age[ia];
	(CMHK.Ginst[iz][idum]).fOIII=fOIII[ia];
	(CMHK.Ginst[iz][idum]).fOII=fOII[ia];
	(CMHK.Ginst[iz][idum]).fOI=fOI[ia];
	(CMHK.Ginst[iz][idum]).fHa=fHa[ia];
	(CMHK.Ginst[iz][idum]).fHb=fHb[ia];
	(CMHK.Ginst[iz][idum]).fNII=fNII[ia];
	idum++;
      }
    }
    if(idum!=CMHK.nage) {
      printf(" There was an error 7 reading CMHK models. Exiting.\n");
      exit(1);
    }
  }  
  globfree(&pglob);
  sprintf(patternnebem,"%s/nebem_*inst*dat",CMHKdir);
  err=glob(patternnebem,GLOB_ERR,NULL,&pglob); 
  if(err) {
    printf(" There was an error 8 reading CMHK models. Exiting.\n");
    exit(1);
  }
  if(CMHK.nmet!=pglob.gl_pathc ) {
    printf(" There was an error 9 reading CMHK models. Exiting.\n");
    exit(1);
  }
  for(iz=0;iz<CMHK.nmet;iz++) {
    if(DEBUG) printf("Readgin %s\n",pglob.gl_pathv[iz]);
    ReadNumcol(pglob.gl_pathv[iz],4,ewHb,ilog,&nfile);
    idum=0;
    for(ia=0;ia<nfile;ia++) {
      if(ilog[ia]) {
	(CMHK.Ginst[iz][idum]).ewHb=ewHb[ia];
	idum++;
      }
    }
    if(idum!=CMHK.nage) {
      printf(" There was an error 10 reading CMHK models. Exiting.\n");
      exit(1);
    }
  }    


}
