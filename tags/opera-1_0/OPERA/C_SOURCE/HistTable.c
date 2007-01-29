#include "modulos.h"

#define  NGRAPH 10
#define  NBUF 6

struct lab {
  char xlabel[100];
  char ylabel[100];
  char title[100];
  int font;
  float height;
};

struct ley {
  int active;
  float x;
  float y;
  float height;
  int font;
};

struct lim {
  float xmin;
  float xmax;
  float ymin;
  float ymax;
};

struct dat {
  int status;
  char file[200];
  char label[100];
  int xcol;
  float *xdata;
  int xnlin;
  int *xtest;
};

struct opr {
  char xlog;
  char opt0;
  float xplus;
};

struct sty {
  int nbin;
  int filled;
  int color;
  int width;
  int listy;
};

struct dis {
  int distid; /* 0 = none, 1 = poisson, 2 = gaussian */
  float param1;
  float param2;
};


struct graph {
  int status;
  int pgid;
  struct lab labels;
  struct ley leyend;
  struct lim limits;
  struct dat data[NBUF];
  struct sty style[NBUF];
  struct opr operac[NBUF];
  struct dis disti;
};

char parfilename[100];

void LoadData(struct graph *grafica,int idat);
void PlotGraphs(struct graph grafica);
void InicData(struct graph *graficas);
void ChangeLimits(struct graph *grafica);
void ChangeLeyend(struct graph *grafica);
void ChangeLabels(struct graph *grafica);
/* //void Information(struct graph *grafica); */
void ChangeMarkers(struct graph *grafica);
void ReadOptions(struct graph *graficas,FILE *filepar, char file[51]);
void SaveOptions(struct graph *graficas);
void PerformOper(struct graph *graficas);
void TrueData(float *xdata,struct dat data,struct opr oper,int *ndata);
void CursorMode(struct graph grafica);
void OverPlotDist(struct graph *grafica);
int NumberData(int nj,struct dat data,struct opr oper);
int main(int argc, char **argv)
{
  struct graph graficas[NGRAPH];
  int i,current,current2;
/*   int active; */
  char option[1]="1";
  FILE *filepar;
  InicData(graficas);
  if(argc < 2) {
    LoadData(graficas,0);
    graficas[0].status=1;
  }
  else {
    if((filepar=fopen(argv[1],"r"))==NULL) {
      printf("ERROR: Can't open options file %s\n",argv[1]);
      exit(1);
    }
    strcpy(parfilename,argv[1]);
    ReadOptions(graficas,filepar, argv[1]);
  }
  cpgopen("?");
  cpgask(0);
  PlotGraphs(graficas[0]);

  current=0;

  /* Aqui empieza el menu */
  while(option[0]!='E') {
    printf("\n D   Load Data\n");
    printf(" M   Change style  \n");
    printf(" L   Change labels\n");
    printf(" I   Change limits\n");
    printf(" Y   Change leyend\n");
    printf(" O   Perform operations on data\n");
    printf(" H   Overplot distribution\n");
    printf(" C   Cursor mode\n");
    printf(" A   Change active graph\n");
    printf(" F   Print information\n");
    printf(" S   Save options\n");
    printf(" E   Exit\n");
    printf("\n\nEnter option: ");
    fflush(NULL);
    setvbuf(stdin,"",_IOLBF,0);
    setvbuf(stdout,"",_IOLBF,0);
    
    scanf("%s",&option[0]);
    current2=current;
  
    switch (option[0]) {
    case 'D' : 
    case 'd' :
      for(i=0;i<NBUF;i++) {
	if(!graficas[current].data[i].status)
	  {
	    printf("Input buffer:");
	    i=readi(i);
	    i--;
	    LoadData(graficas+current,i);
	    break;
	  }
      }
      break;
    case 'M' :
    case 'm' :
      ChangeMarkers(graficas+current);
      break;
    case 'L' :
    case 'l' :
      ChangeLabels(graficas+current);
      break;
    case 'I' : 
    case 'i' :
      ChangeLimits(graficas+current);
      break;
    case 'S' : 
    case 's' :
      SaveOptions(graficas);
      break;
    case 'E' : 
    case 'e' :
      cpgend();
      exit(0);
      break;
    case 'O' : 
    case 'o' :
      PerformOper(graficas+current);
      break;
    case 'A' : 
    case 'a' :
      printf("Input graph 1-%d :",NGRAPH);
      current=readi(current+1);
      current--;
      if(graficas[current].status==0) LoadData(graficas+current,0);
      graficas[current].status=1;
      current2=current;
      break;
    case 'C' :
    case 'c' :
      CursorMode(graficas[current]);
      break;
    case 'Y' : 
    case 'y' :
      ChangeLeyend(graficas+current);
      break;
    case 'H' :
    case 'h':
      OverPlotDist(graficas+current);
      break;
/* //    case 'F' : printf("Information about histogram %d",active); */
/* //      Information(graficas+active); */
/* //      break; */
    }
    PlotGraphs(graficas[current]);
  }

  return(0);  
}

void InicData(struct graph *graficas)
{
  int i,j;
  
  for(i=0;i<NGRAPH;i++) {
    graficas[i].status=0;
    graficas[i].labels.font=2;
    graficas[i].labels.height=1.8;
    graficas[i].leyend.active=0;
    graficas[i].limits.xmin=0;
    graficas[i].limits.xmax=0;
    graficas[i].limits.ymin=0;
    graficas[i].limits.ymax=0;
    strcpy(graficas[i].labels.xlabel,"\0");
    strcpy(graficas[i].labels.ylabel,"\0");
    strcpy(graficas[i].labels.title,"\0");
    for(j=0;j<NBUF;j++) {
      graficas[i].style[j].nbin=10;
      graficas[i].style[j].filled=0;
      graficas[i].style[j].width=1;
      graficas[i].style[j].color=1;
      graficas[i].style[j].listy=1;
      graficas[i].data[j].status=0;
      graficas[i].operac[j].xlog='N';
      graficas[i].operac[j].opt0='N';
      graficas[i].operac[j].xplus=0;

      strcpy(graficas[i].data[j].file,"\0");

    }
  }
  
}
   
      

void LoadData(struct graph *grafica,int idat)
{
/*   int i; */
  int nlin;
/*   float x[200]; */
/*   int lx[200]; */
  (*grafica).data[idat].status=1;
  
  printf("Enter file with data :");
  reads((*grafica).data[idat].file,(*grafica).data[idat].file);
  printf("Enter column with x data :");
  (*grafica).data[idat].xcol=readi((*grafica).data[idat].xcol);
  /* scanf("%d",&(*grafica).data[idat].xcol);   */

  
  nlin=FileNLin((*grafica).data[idat].file);
  
  (*grafica).data[idat].xdata=malloc(nlin*sizeof(float));
  (*grafica).data[idat].xtest=malloc(nlin*sizeof(int));
  
  ReadNumcol((*grafica).data[idat].file,(*grafica).data[idat].xcol,
	     (*grafica).data[idat].xdata,(*grafica).data[idat].xtest,
	     &(*grafica).data[idat].xnlin);
 
  

}


void PlotGraphs(struct graph grafica)
{
  struct data {
    float *xdat;
    int ndat;
  };
  float yley,yleyr,nul,xley;
  int i,j;
  int fill=1;
  float xmin,xmax,ymin,ymax,xminh,xmaxh;
  struct data data[NBUF];
  int nplotdist=200;
  float xplot;
  double distvalue;
  cpgpage();
  i=0;
  while(grafica.data[i].status){
    
    (data[i].xdat)=malloc((grafica.data[i].xnlin)*sizeof(float));
    TrueData(data[i].xdat,grafica.data[i],grafica.operac[i],&(data[i].ndat));    
    i++;
  }
  
  if((grafica.limits.xmin == 0. && grafica.limits.xmax == 0.)) {
    xmin=1e30;
    xmax=-1e30;
    i=0;
    while(grafica.data[i].status){
      for(j=0;j<data[i].ndat;j++) {       
	if(xmax<data[i].xdat[j]) xmax=data[i].xdat[j];
	if(xmin>data[i].xdat[j]) xmin=data[i].xdat[j];
      }
      i++;
    }
  }
  else {
    xmin=grafica.limits.xmin;
    xmax=grafica.limits.xmax;
  }
  if((grafica.limits.ymin == 0. && grafica.limits.ymax == 0.)) {
    i=0;
    ymax=0;
    ymin=0;
    while(grafica.data[i].status) {
      if (ymax<data[i].ndat) ymax=data[i].ndat/grafica.style[i].nbin*3;
      i++;
    }
  }
  else {
    ymin=grafica.limits.ymin;
    ymax=grafica.limits.ymax;
  }
  
  cpgsls(1);
  cpgslw(2);
  cpgscf(grafica.labels.font);
  cpgsch(grafica.labels.height);
  cpgvstd();
  cpgswin(xmin,xmax,ymin,ymax);
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
  cpglab(grafica.labels.xlabel,grafica.labels.ylabel,grafica.labels.title);
  printf("\n X limits: %f %f\n",xmin,xmax);
  printf(" Y limits: %f %f\n",ymin,ymax);
  xley=grafica.leyend.x+(xmax-xmin)/18; 
  i=0;
  while(grafica.data[i].status) {
    printf(" Buffer %d \n",i);
    
    cpgsci(grafica.style[i].color);
    cpgslw(grafica.style[i].width);
    cpgsls(grafica.style[i].listy);
    if(grafica.style[i].filled==2) {
      fill=3;
      cpgsfs(1); }
    if(grafica.style[i].filled==4) {
      fill=3;
      cpgsfs(3); }
    if(grafica.style[i].filled==5) {
      fill=3;
      cpgsfs(4); }
    if(grafica.style[i].filled==1) {
      fill=1; }
    if(grafica.style[i].filled==3) {
      fill=5; }
    
    
    if(xmin>xmax) {
      xminh=xmax;
      xmaxh=xmin; }
    else {
      xminh=xmin;
      xmaxh=xmax; }
    
    if(grafica.leyend.active) {
      cpgscf(grafica.leyend.font);
      cpgsch(grafica.leyend.height);
      cpglen(4,"H",&nul,&yleyr);
      cpgsci(1);
      yley=grafica.leyend.y-yleyr*1.8*i;
      cpgptxt(xley,yley-yleyr/2,0.0,0.0,grafica.data[i].label);
      cpgsci(grafica.style[i].color);
      cpgmove(grafica.leyend.x,yley);
      cpgdraw(grafica.leyend.x+(xmax-xmin)/20,yley);
    }

    cpghist(data[i].ndat,data[i].xdat,xminh,xmaxh,grafica.style[i].nbin,fill);
    
    if(grafica.disti.distid!=0) {
      printf(" PINTO\n");
      cpgmove(xmin,ymin);
      for(j=0;j<nplotdist;j++) {
	xplot=xmin+j*(xmax-xmin)/(nplotdist-1);
	if(grafica.disti.distid==1) distvalue=poidist(xplot,grafica.disti.param1);
	else if(grafica.disti.distid==2) distvalue=gaussian(xplot,grafica.disti.param1,grafica.disti.param2);
	else distvalue=0;
/* 	printf(" DISTRI %d\n",grafica.disti.distid); */
 	printf(" xplot %f dis %f yplot %f\n",xplot,distvalue,distvalue*data[i].ndat*(xmax-xmin)/grafica.style[i].nbin);
	cpgdraw(xplot,distvalue*data[i].ndat*(xmax-xmin)/grafica.style[i].nbin);
      }
    }

    cpgsci(1);

    
    printf("   File: %s   Data column: %d\n",grafica.data[i].file,grafica.data[i].xcol);
    printf("   Number of data: %d  Color: %d  Width: %d\n",data[i].ndat,grafica.style[i].color,grafica.style[i].width);
    free(data[i].xdat);
    i++;
  }
}
void ChangeLeyend(struct graph *grafica)
{
/*   char nul[100]; */
  char opt;
  int i;
  printf("Active leyend (Y/N) [Y] :");
  fflush(stdin);
  setvbuf(stdin,"",_IOLBF,0);
  setvbuf(stdout,"",_IOLBF,0);
  
  opt=readc('N');
  if (opt=='N')    (*grafica).leyend.active=0;
  else {
    (*grafica).leyend.active=1;
    printf("Input X position:");
    (*grafica).leyend.x=readf((*grafica).leyend.x);
    /* scanf("%f",&(*grafica).leyend.x); */
    printf("Input Y position: ");
    (*grafica).leyend.y=readf((*grafica).leyend.y);
    /*     scanf("%f",&(*grafica).leyend.y); */
    
    printf("Input Font:");
    (*grafica).leyend.font=readi((*grafica).leyend.font);
/*     scanf("%d",&(*grafica).leyend.font ); */
    printf("Input character height [%f] :",(*grafica).leyend.height);
    (*grafica).leyend.height=readf((*grafica).leyend.height);
/*     scanf("%f",&(*grafica).leyend.height); */
    for(i=0;i<NBUF;i++) {
      if ((*grafica).data[i].status) {
	printf("Input label for buffer %d ",i);
	fflush(stdin);
	setvbuf(stdin,"",_IOLBF,0);
	setvbuf(stdout,"",_IOLBF,0);
	
	reads((*grafica).data[i].label,(*grafica).data[i].label);
/* 	//gets(nul); */
/* 	//strcpy((*grafica).data[i].label,nul); */
      }
    }
  }
}

void ChangeLabels(struct graph *grafica)
{
/*   char nul[100]; */
  printf("Input X label ");
  reads((*grafica).labels.xlabel,(*grafica).labels.xlabel);
/*   //fflush(stdin); */
/*   //gets(nul); */
/*   //strcpy((*grafica).labels.xlabel,nul); */
  printf("Input Y label ");
  reads((*grafica).labels.ylabel,(*grafica).labels.ylabel);
/*   //fflush(stdin); */
/*   //gets(nul); */
/*   //strcpy((*grafica).labels.ylabel,nul); */
  fflush(stdin);
  setvbuf(stdin,"",_IOLBF,0);
  setvbuf(stdout,"",_IOLBF,0);
  
  printf("Input graph title ");
  reads((*grafica).labels.title,(*grafica).labels.title);
/*   //fflush(stdin); */
/*   //gets(nul); */
/*   //strcpy((*grafica).labels.title,nul); */
  fflush(stdin);
  setvbuf(stdin,"",_IOLBF,0);
  setvbuf(stdout,"",_IOLBF,0);
  
  printf("Input character font 1-4 [%d]: ",(*grafica).labels.font);
  (*grafica).leyend.font=readi((*grafica).leyend.font);
/*   scanf("%d",&(*grafica).labels.font); */
  printf("Input character height [%f]: ",(*grafica).labels.height);
  (*grafica).leyend.height=readf((*grafica).leyend.height);
  /* scanf("%f",&(*grafica).labels.height); */
}
void ChangeLimits(struct graph *grafica)
{
  printf("Input xmin: ");
  (*grafica).limits.xmin=readf((*grafica).limits.xmin);
/*   //scanf("%f",&(*grafica).limits.xmin); */
  printf("Input xmax: ");
  (*grafica).limits.xmax=readf((*grafica).limits.xmax);
/*   //scanf("%f",&(*grafica).limits.xmax); */
  printf("Input ymin: ");
  (*grafica).limits.ymin=readf((*grafica).limits.ymin);
/*   //scanf("%f",&(*grafica).limits.ymin); */
  printf("Input ymax: ");
  (*grafica).limits.ymax=readf((*grafica).limits.ymax);
/*   //scanf("%f",&(*grafica).limits.ymax); */
  
}


/* //void Information(struct graph *grafica) */
/* //{ */
/* //  int i; */
/* //  for(i=0;i<NBUF;i++) { */

void ChangeMarkers(struct graph *grafica)
{
  int buffer=1;
  printf("Input buffer 1-%d [1] :",NBUF);
  buffer=readi(buffer);
  buffer--;
  printf("Input number of bins 0-31 :");
  (*grafica).style[buffer].nbin=readi((*grafica).style[buffer].nbin);
  printf("Input color 0-15  :");
  (*grafica).style[buffer].color=readi((*grafica).style[buffer].color);
  printf("Input width 0-15 :");
  (*grafica).style[buffer].width=readi((*grafica).style[buffer].width);
  printf("\n Line options:\n 1 Full line\n 2 Dashed\n 3 Dot-dash-dot-dash\n 4 Dotted\n 5 Dash-dot-dot-dot\n");
  printf("Input line style 0-5 :");
  (*grafica).style[buffer].listy=readi((*grafica).style[buffer].listy);
  printf("\n Filled options: \n 1 Normal\n 2 Filled\n 3 Line drawn\n 4 Dashed 45\n 5 Dashed 135\n");
  printf("Filled histogram option 1-5 :");
  (*grafica).style[buffer].filled=readi((*grafica).style[buffer].filled);
  
}


void SaveOptions(struct graph *graficas)
{
  FILE *fp;
  int i,j;
  int nc,nt;
/*   char nul[2]; */
/*   char file[100]; */
  char ch21[51];
  printf("Name of file  ");
  reads(parfilename,parfilename);
  if((fp=fopen(parfilename,"w")) ==NULL) {
    printf("ERROR: Can't open file\n");
    return;
  }
  fprintf(fp,"COMMENT  Parameter file for HistoTable                                         \n");
  nc=1;
  for(i=0;i<NGRAPH;i++) {
    
    for(j=0;j<NBUF;j++) {
      if(graficas[i].data[j].status) {
	sprintf(ch21,"'%s'",graficas[i].data[j].file);
	fprintf(fp,"FILE%1d_%1d = %-31.31s / File for histogram %d, buffer %d     \n",i,j,ch21,i,j);
	fprintf(fp,"XCOL%1d_%1d =%21d / X column for graph %d, buffer %d                \n",i,j,graficas[i].data[j].xcol,i,j);
	fprintf(fp,"NBIN%1d_%1d =%21d / Number of bins for histogram %d, buffer %d      \n",i,j,graficas[i].style[j].nbin,i,j);
	fprintf(fp,"FILL%1d_%1d =%21d / Filled option for graph %d, buffer %d           \n",i,j,graficas[i].style[j].filled,i,j);
	fprintf(fp,"COLOR%1d_%1d=%21d / Color marker for graph %d, buffer %d            \n",i,j,graficas[i].style[j].color,i,j);
        fprintf(fp,"WIDTH%1d_%1d=%21d / Line width for histogram %d, buffer %d          \n",i,j,graficas[i].style[j].width,i,j);
        fprintf(fp,"LISTY%1d_%1d=%21d / Line style for histogram %d, buffer %d          \n",i,j,graficas[i].style[j].listy,i,j);

	fprintf(fp,"XLOG%1d_%1d =%21c / Apply logs in X for graph %d, buffer %d (Y/N)   \n",i,j,graficas[i].operac[j].xlog,i,j);
	fprintf(fp,"OPT0%1d_%1d =%21c / Aply log for values <=0 for graph %d, buffer %d \n",i,j,graficas[i].operac[j].opt0,i,j);
	nc += 9;
	if(graficas[i].leyend.active) {
	  sprintf(ch21,"'%s'",graficas[i].data[j].label);
	  fprintf(fp,"LEYT%1d_%1d = %-31.31s / Leyend text for graph %d, buffer %d  \n",i,j,ch21,i,j);  
	  nc++;
	}
      }
	
    }
    if(graficas[i].leyend.active) {
      fprintf(fp,"XLEYEN_%1d=%21f / X position for leyend in graph %d              \n",i,graficas[i].leyend.x,i);
      fprintf(fp,"YLEYEN_%1d=%21f / Y position for leyend in graph %d              \n",i,graficas[i].leyend.y,i);
      fprintf(fp,"LEYFON_%1d=%21d / Font for leyend in graph %d                    \n",i,graficas[i].leyend.font,i);
      fprintf(fp,"LEYHEI_%1d=%21f / Height for labels in graph %d                  \n",i,graficas[i].leyend.height,i);
      nc += 4;
      
    }
    sprintf(ch21,"'%s'",graficas[i].labels.xlabel);
    fprintf(fp,"XLABEL_%1d= %-51.51s / X Label graph %d\n",i,ch21,i);
    sprintf(ch21,"'%s'",graficas[i].labels.ylabel);
    fprintf(fp,"YLABEL_%1d= %-51.51s / Y Label graph %d\n",i,ch21,i);
    sprintf(ch21,"'%s'",graficas[i].labels.title);
    fprintf(fp,"TITLE__%1d= %-51.51s / Title   graph %d\n",i,ch21,i);
    fprintf(fp,"LABFON_%1d=%21d / Font for labels in graph %d                    \n",i,graficas[i].labels.font,i);
    fprintf(fp,"LABHEI_%1d=%21f / Height for labels in graph %d                  \n",i,graficas[i].labels.height,i);
    fprintf(fp,"XMIN_%1d  =%21f / Xmin in graph %d                               \n",i,graficas[i].limits.xmin,i);
    fprintf(fp,"XMAX_%1d  =%21f / Xmax in graph %d                               \n",i,graficas[i].limits.xmax,i);
    fprintf(fp,"YMIN_%1d  =%21f / Ymin in graph %d                               \n",i,graficas[i].limits.ymin,i);
    fprintf(fp,"YMAX_%1d  =%21f / Ymax in graph %d                               \n",i,graficas[i].limits.ymax,i);
    nc += 9;
    

  }

  fprintf(fp,"COMMENT  Blank lines to assure the last block of 2880 caracters is complete    \n");
  fprintf(fp,"COMMENT                                                                        \n");
  nc += 2;
  for(nt=nc;nt<(1+(int)(nc/36))*36-1;nt++) {
    fprintf(fp,"COMMENT                                                                        \n");
  }
  fprintf(fp,"END                                                                            \n");
  fclose(fp);
}

void ReadOptions(struct graph *graficas,FILE *filepar, char file[51])
{
  int i,j;
  int nlin;
  char nul3[3],nul1[1];
  char keyf[9]="",key[9]="";
  int status=0;
  char comment[51];
  fitsfile *parfile;
  printf("Using parameter file <<%s>>\n",file);
  if( ffopen2(&parfile,file, READONLY, &status)) fits_report_error(stderr,status);

  for(i=0;i<NGRAPH;i++) {
    sprintf(nul1,"%1d",i);
    strcpy(key,"XLABEL_");
    strcat(key,nul1);
/*     f_kfls(filepar,key,graficas[i].labels.xlabel); */
    ffgky(parfile,TSTRING,key,graficas[i].labels.xlabel,comment,&status);
    if(status) printf(" keyword not found: %s\n",key);
    fits_report_error(stderr,status);
    status=0;

    strcpy(key,"YLABEL_");
    strcat(key,nul1);
/*     f_kfls(filepar,key,graficas[i].labels.ylabel); */
    ffgky(parfile,TSTRING,key,graficas[i].labels.ylabel,comment,&status);
    if(status) printf(" keyword not found: %s\n",key);
    fits_report_error(stderr,status);
    status=0;

    strcpy(key,"TITLE__");
    strcat(key,nul1);
/*     f_kfls(filepar,key,graficas[i].labels.title); */
    ffgky(parfile,TSTRING,key,graficas[i].labels.title,comment,&status);
    if(status) printf(" keyword not found: %s\n",key);
    fits_report_error(stderr,status);
    status=0;    

    strcpy(key,"LABFON_");
    strcat(key,nul1);
/*     f_kfi(filepar,key,&(graficas[i].labels.font)); */
    ffgky(parfile,TINT,key,&(graficas[i].labels.font),comment,&status);
    if(status) printf(" keyword not found: %s\n",key);
    fits_report_error(stderr,status);
    status=0;

    strcpy(key,"LABHEI_");
    strcat(key,nul1);
/*     f_kff(filepar,key,&(graficas[i].labels.height));     */
    ffgky(parfile,TFLOAT,key,&(graficas[i].labels.height),comment,&status);  
    if(status) printf(" keyword not found: %s\n",key);
    fits_report_error(stderr,status);
    status=0;

    strcpy(key,"XMIN_");
    strcat(key,nul1);
/*     f_kff(filepar,key,&(graficas[i].limits.xmin)); */
    ffgky(parfile,TFLOAT,key,&(graficas[i].limits.xmin),comment,&status);  
    if(status) printf(" keyword not found: %s\n",key);
    fits_report_error(stderr,status);
    status=0;

    strcpy(key,"XMAX_");
    strcat(key,nul1);
/*     f_kff(filepar,key,&(graficas[i].limits.xmax)); */
    ffgky(parfile,TFLOAT,key,&(graficas[i].limits.xmax),comment,&status);  
    if(status) printf(" keyword not found: %s\n",key);
    status=0;
    strcpy(key,"YMIN_");
    strcat(key,nul1);
/*     f_kff(filepar,key,&(graficas[i].limits.ymin)); */
    ffgky(parfile,TFLOAT,key,&(graficas[i].limits.ymin),comment,&status);  
    if(status) printf(" keyword not found: %s\n",key);
    status=0;
    strcpy(key,"YMAX_");
    strcat(key,nul1);
/*     f_kff(filepar,key,&(graficas[i].limits.ymax)); */
    ffgky(parfile,TFLOAT,key,&(graficas[i].limits.ymax),comment,&status);  
    if(status) printf(" keyword not found: %s\n",key);
    status=0;

    
    strcpy(key,"XLEYEN_");
    strcat(key,nul1);
    ffgky(parfile,TFLOAT,key,&(graficas[i].leyend.x),comment,&status);
    
    
    if(!status) {
      graficas[i].leyend.active=1;
      strcpy(key,"YLEYEN_");
      strcat(key,nul1);
      ffgky(parfile,TFLOAT,key,&(graficas[i].leyend.y),comment,&status);  
      /*       //printf(" HAY LA YA en y=%f\n",graficas[i].leyend.y); */
	if(status) printf(" keyword not found: %s\n",key);
	status=0;
	strcpy(key,"LEYFON_");
	strcat(key,nul1);
	ffgky(parfile,TINT,key,&(graficas[i].leyend.font),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;
	strcpy(key,"LEYHEI_");
	strcat(key,nul1);
	ffgky(parfile,TFLOAT,key,&(graficas[i].leyend.height),comment,&status);  
	if(status) printf(" keyword not found: %s\n",key);
	status=0; 
    }
    status=0;

    for(j=0;j<NBUF;j++) {
      strcpy(keyf,"FILE");
      strcpy(nul3,"");
      sprintf(nul3,"%1d_%1d",i,j);
      strcat(keyf,nul3);
      status=0;
      if(ffgky(parfile,TSTRING,keyf,graficas[i].data[j].file,comment,&status)!=0) {
	/*       if((c=f_kfls(filepar,keyf,graficas[i].data[j].file))==0) { */
	graficas[i].data[j].status=0;
      }
      else {
	graficas[i].status=1;
	graficas[i].data[j].status=1;
	strcpy(key,"XCOL");
	strcat(key,nul3);
	ffgky(parfile,TINT,key,&(graficas[i].data[j].xcol),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	if(status) {
	  printf("ERROR: No column data for file %s STOP",graficas[i].data[j].file);
	  exit(1);
	}
	status=0;

/* 	c+=f_kfi(filepar,key,&(graficas[i].data[j].xcol)); */
/* 	if(c<2) { */
/* 	  printf("ERROR: No column data for file %s STOP",graficas[i].data[j].file); */
/* 	  exit(1); */
/* 	} */
	strcpy(key,"NBIN");
	strcat(key,nul3);
/* 	f_kfi(filepar,key,&(graficas[i].style[j].nbin)); */
	ffgky(parfile,TINT,key,&(graficas[i].style[j].nbin),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;

	strcpy(key,"COLOR");
	strcat(key,nul3);
	ffgky(parfile,TINT,key,&(graficas[i].style[j].color),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;
/* 	f_kfi(filepar,key,&(graficas[i].style[j].color)); */

	strcpy(key,"FILL");
	strcat(key,nul3);
/* 	f_kfi(filepar,key,&(graficas[i].style[j].filled)); */
	ffgky(parfile,TINT,key,&(graficas[i].style[j].filled),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;

        strcpy(key,"LISTY");
        strcat(key,nul3);
/*         f_kfi(filepar,key,&(graficas[i].style[j].listy)); */
	ffgky(parfile,TINT,key,&(graficas[i].style[j].listy),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;

        strcpy(key,"WIDTH");
        strcat(key,nul3);
/*         f_kfi(filepar,key,&(graficas[i].style[j].width)); */
	ffgky(parfile,TINT,key,&(graficas[i].style[j].width),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;

	strcpy(key,"XLOG");
	strcat(key,nul3);
/* 	f_kfc(filepar,key,&(graficas[i].operac[j].xlog)); */
	ffgky(parfile,TSTRING,key,&(graficas[i].operac[j].xlog),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	fits_report_error(stderr,status);
	status=0;

	strcpy(key,"OPT0");
	strcat(key,nul3);
/* 	f_kfc(filepar,key,&(graficas[i].operac[j].opt0)); */
	ffgky(parfile,TSTRING,key,&(graficas[i].operac[j].opt0),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;

	if(graficas[i].leyend.active) {
	  strcpy(key,"LEYT");
	  strcat(key,nul3);
/* 	  f_kfs(filepar,key,graficas[i].data[j].label); */
	  ffgky(parfile,TSTRING,key,graficas[i].data[j].label,comment,&status);
	  if(status) printf(" keyword not found: %s\n",key);
	}
	status=0;
	
	nlin=FileNLin(graficas[i].data[j].file);
	
	graficas[i].data[j].xdata=malloc(nlin*sizeof(float));
	graficas[i].data[j].xtest=malloc(nlin*sizeof(int));
	
	
	ReadNumcol(graficas[i].data[j].file,graficas[i].data[j].xcol,
		   graficas[i].data[j].xdata,graficas[i].data[j].xtest,
		   &(graficas[i].data[j].xnlin));
	
 
      }
    }
  }
}

void  PerformOper(struct graph *graficas)
{
  char opt;
  static int buffer=1;
  fflush(stdin);
  setvbuf(stdin,"",_IOLBF,0);
  setvbuf(stdout,"",_IOLBF,0);
  
  printf("Input buffer 1-%d:",NBUF);
  buffer=readi(buffer);
  buffer--;
  printf("Do you want to use logaritms in X data Y/N :");
  fflush(stdin);
  setvbuf(stdin,"",_IOLBF,0);
  setvbuf(stdout,"",_IOLBF,0);
  
  opt=readc('N');
  if(opt=='Y') (*graficas).operac[buffer].xlog='Y';
  else (*graficas).operac[buffer].xlog='N';
  printf("Do you want to apply logaritms to values <=0 Y/N :");
  fflush(stdin);
  setvbuf(stdin,"",_IOLBF,0);
  setvbuf(stdout,"",_IOLBF,0);
  
  opt=readc('N');
  if(opt=='Y') (*graficas).operac[buffer].opt0='Y';
  else (*graficas).operac[buffer].opt0='N';

/* //  printf("Enter constant to sum to X data: "); */
/* //  scanf("%f",&(*graficas).operac[buffer].xplus); */
/* //  printf("Enter constant to sum to Y data: "); */
/* //  scanf("%f",&(*graficas).operac[buffer].xplus); */

}


void TrueData(float *xdata,struct dat data,struct opr oper,int *ndata)
{
  int n=0;
  int j;
  for(j=0;j<data.xnlin;j++) {
    if(data.xtest[j]  ) {
      xdata[n]=data.xdata[j]+oper.xplus;
      if(oper.xlog=='Y') {
	if(xdata[n]<=0) {
	  if(oper.opt0=='Y') xdata[n]=0;
	  else continue;
	}
	else xdata[n]=log10(xdata[n]);
      }
      n++;
    }
  }
  *ndata=n;
}

void CursorMode(struct graph grafica)
{
  struct data {
    float *xdat;
    float *ydat;	
    int ndat;
  };
  float x,y;
  char ch='A';
  float rmin,r;
  int i,j;
  int ni,nj,njj;
  struct data data[NBUF];
  while(ch=='A')
    {
      cpgcurs(&x,&y,&ch);
      ni=0;
      nj=0;
      rmin=1e30;      
      i=0;
      while(grafica.data[i].status) {
	
	(data[i].xdat)=malloc((grafica.data[i].xnlin)*sizeof(float));
	(data[i].ydat)=malloc((grafica.data[i].xnlin)*sizeof(float));
	TrueData(data[i].xdat,grafica.data[i],grafica.operac[i],&(data[i].ndat));    
	for(j=0;j<data[i].ndat;j++) {
	  r=(x-data[i].xdat[j])*(x-data[i].xdat[j])+(y-data[i].ydat[j])*(y-data[i].ydat[j]);
	  if(r<rmin) {
	    ni=i;
	    nj=j;
	    rmin=r;
	  }
	}
	i++;

      }
      njj=NumberData(nj,grafica.data[ni],grafica.operac[ni]);
      printf("X cursor: %f\n",x);  
      printf("X data: %f\n",data[ni].xdat[nj]);
      printf("X column: %e\n",grafica.data[ni].xdata[njj]);
    }
}


int NumberData(int nj,struct dat data,struct opr oper)
{
  int n=-1;
  int j;
  float  xtemp;
  for(j=0;j<data.xnlin;j++) {
    if(data.xtest[j]  ) {
      xtemp=data.xdata[j]+oper.xplus;
      
      if(oper.xlog=='Y') {
	if(xtemp<=0) {
	  if(oper.opt0=='Y') ;
	  else continue;
	}
      }
      n++;
    }
    if(nj==n) return(j);
  }
  return(0);
}



void OverPlotDist(struct graph *grafica) {

  printf(" Type of distribution? (0=none, 1=Poisson, 2=Gaussian) ");
  (*grafica).disti.distid=readi((*grafica).disti.distid);
  if((*grafica).disti.distid==1) {
    printf(" Type lmabda parameter (mean) ");
    (*grafica).disti.param1=readf((*grafica).disti.param1);
  }
  else if((*grafica).disti.distid==2) {
    printf(" Type mean parameter ");
    (*grafica).disti.param1=readf((*grafica).disti.param1);
    printf(" Type sigma parameter ");
    (*grafica).disti.param2=readf((*grafica).disti.param2);
  }

}
