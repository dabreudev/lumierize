/* #include "cpgplot.h" */
#include "modulos.h"

#define  NGRAPH 10
#define  NBUF 6
#define NSTR 301

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
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double ratio;
};

struct dat {
  int status;
  char file[200];
  char label[100];
  int xcol;
  int ycol;
  int idcol;
  double *xdata;
  int xnlin;
  int *xtest;
  double *ydata;
  int ynlin;
  int *ytest;
  char *id;
  int idnlin;
  int *idtest;
};

struct opr {
  char xlog;
  char ylog;
  char opt0;
  float xplus;
  float yplus;
};

struct sym {
  int marker;
  float height;
  int color;
  int lw;
  int ls;
/*   int sort; */
};

struct lin {
  int status;
  float a;
  float b;
  int color;
  int width;
};

struct graph {
  int status;
  int pgid;
  struct lab labels;
  struct ley leyend;
  struct lim limits;
  struct dat data[NBUF];
  struct sym symbol[NBUF];
  struct opr operac[NBUF];
  struct lin lines[NBUF];
 };

char parfilename[100];

void LoadData(struct graph *grafica,int idat);
void PlotGraphs(struct graph grafica);
void InicData(struct graph *graficas);
void ChangeLimits(struct graph *grafica);
void ChangeLimitsMouse(struct graph *grafica);
void ChangeLabels(struct graph *grafica);
void ChangeMarkers(struct graph *grafica);
void ChangeLeyend(struct graph *grafica);
void LoadParam(struct graph *graficas,FILE *filepar, char file[50]);
void SaveParam(struct graph *graficas);
void PerformOper(struct graph *graficas);
void DrawLines(struct graph *graficas);
void TrueData(double *xdata,double * ydata,struct dat data,struct opr oper,int *ndata);
void CursorMode(struct graph grafica);
int NumberData(int nj,struct dat data,struct opr oper);
int main(int argc, char **argv)
{
  struct graph graficas[NGRAPH];
  int i,current,current2;
  char option[1]="m";
/*   char opt; */
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
/*     printf(" ANTES\n"); */
    LoadParam(graficas, filepar, argv[1]);
/*      printf(" DESPU\n"); */
  }

/*   printf(" PETO\n"); */

  cpgopen("?");
  cpgask(0);

  PlotGraphs(graficas[0]);

  current=0;

  /* Aqui empieza el menu */
  while(option[0]!='E') {
    printf("\n D   Load Data\n");
    printf(" M   Change markers\n");
    printf(" L   Change labels\n");
    printf(" I   Change limits (keyboard)\n");
    printf(" X   Change limits (mouse)\n");
    printf(" Y   Change leyend\n");
    printf(" N   Draw lines\n");
    printf(" O   Perform operations on data\n");
    printf(" C   Cursor mode\n");
    printf(" A   Change active graph\n");
    printf(" S   Save options\n");
    printf(" E   Exit\n");
    printf("\n Active graph %d\n",current);
    printf("\n\nEnter option: ");


    fflush(NULL);
    setvbuf(stdin,"",_IOLBF,0);
    setvbuf(stdout,"",_IOLBF,0);
    scanf("%s",&option[0]);
/*     opt=getchar(); */
    current2=current;
    switch (option[0]) {
/*       switch (opt) { */
    case 'D' : 
    case 'd' :
      for(i=0;i<NBUF;i++) {
	if(!graficas[current].data[i].status)
	  {
/* 	    printf("Input buffer [%d]:",i); */
/* 	    scanf("%d",&i); */
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
    case 'X' : 
    case 'x' :
      ChangeLimitsMouse(graficas+current);
      break;
    case 'S' : 
    case 's' :
      SaveParam(graficas);
      break;
    case 'E' : 
    case 'e' :
      cpgend();
      exit(0);
      break;
    case 'O' : 
    case 'o' :
/*       printf(" Grafica activge %d  %d\n",current,current2); */
      PerformOper(graficas+current);      
/*       PerformOper(&(graficas[current])); */
/*       printf(" Menu princiapal act%d  %d\n",current,current2); */
      break;
    case 'A' : 
    case 'a' :
      printf("Input graph 1-%d :",NGRAPH);
      current=readi(current+1);
/*       //scanf("%d",&current); */
      current--;
      if(graficas[current].status==0) LoadData(graficas+current,0);

      graficas[current].status=1;
      current2=current;
      break;
    case 'N' : 
    case 'n' :
      DrawLines(graficas+current);
      break;
    case 'C' : 
    case 'c' :
      CursorMode(graficas[current]);
      break;
    case 'Y' : 
    case 'y' :
      ChangeLeyend(graficas+current);
      break;
    }
/*     current=current2; */
/*     printf(" YA ESTA 1111   current %d\n",current); */
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
    graficas[i].labels.height=1.3;
    graficas[i].leyend.active=0;
    graficas[i].limits.xmin=0;
    graficas[i].limits.xmax=0;
    graficas[i].limits.ymin=0;
    graficas[i].limits.ymax=0;
    graficas[i].limits.ratio=0;
    strcpy(graficas[i].labels.xlabel,"\0");
    strcpy(graficas[i].labels.ylabel,"\0");
    strcpy(graficas[i].labels.title,"\0");
    for(j=0;j<NBUF;j++) {
      graficas[i].symbol[j].marker=0;
      graficas[i].symbol[j].height=1.2;
      graficas[i].symbol[j].color=1;
      graficas[i].symbol[j].lw=0;
      graficas[i].symbol[j].ls=1;      
      graficas[i].data[j].status=0;
      graficas[i].data[j].xcol=1;
      graficas[i].data[j].ycol=2;
      graficas[i].data[j].idcol=0;
      graficas[i].operac[j].xlog='N';
      graficas[i].operac[j].ylog='N';
      graficas[i].operac[j].opt0='N';
      graficas[i].operac[j].xplus=0;
      graficas[i].operac[j].yplus=0;
      graficas[i].lines[j].status=0;
      graficas[i].lines[j].width=1;
      graficas[i].lines[j].color=1;

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
  
  printf("Enter file with data: ");
  reads((*grafica).data[idat].file,(*grafica).data[idat].file);
  printf("Enter column with x data: ");
  (*grafica).data[idat].xcol=readi((*grafica).data[idat].xcol);
/*   scanf("%d",&(*grafica).data[idat].xcol);   */
  printf("Enter column with y data: ");
  (*grafica).data[idat].ycol=readi((*grafica).data[idat].ycol);
  printf("Enter column with id data (0=no id): ");
  (*grafica).data[idat].idcol=readi((*grafica).data[idat].idcol);


  
  nlin=FileNLin((*grafica).data[idat].file);
  
  (*grafica).data[idat].xdata=malloc(nlin*sizeof(double));
  (*grafica).data[idat].ydata=malloc(nlin*sizeof(double));
  (*grafica).data[idat].id=malloc(NSTR*nlin*sizeof(char));
  (*grafica).data[idat].xtest=malloc(nlin*sizeof(int));
  (*grafica).data[idat].ytest=malloc(nlin*sizeof(int));
  (*grafica).data[idat].idtest=malloc(nlin*sizeof(int));
  
  ReadDoublecol((*grafica).data[idat].file,(*grafica).data[idat].xcol,
		(*grafica).data[idat].xdata,(*grafica).data[idat].xtest,
		&(*grafica).data[idat].xnlin);
  
  ReadDoublecol((*grafica).data[idat].file,(*grafica).data[idat].ycol,
		(*grafica).data[idat].ydata,(*grafica).data[idat].ytest,
		&(*grafica).data[idat].ynlin);
  if((*grafica).data[idat].idcol!=0) ReadCharcol((*grafica).data[idat].file,(*grafica).data[idat].idcol,
						 (*grafica).data[idat].id,(*grafica).data[idat].idtest,NSTR,
						 &(*grafica).data[idat].idnlin);
  if((*grafica).data[idat].xnlin!=(*grafica).data[idat].ynlin) {
    printf("ERROR: x dat and y data hasn't the same size");
    (*grafica).data[idat].status=0;
  }
  

}


void PlotGraphs(struct graph grafica)
{
  struct data {
    double *xdat;
    double *ydat;	
    int ndat;
  };
  int i,j;
  float ran,xley,yley,yleyr;
  double xmin,xmax,ymin,ymax;
  float nul;
  struct data data[NBUF];
/*   printf(" AQUI si 22222\n"); */
  cpgpage();
  i=0;
  while(grafica.data[i].status && i<NBUF) {
    
    (data[i].xdat)=malloc((grafica.data[i].xnlin)*sizeof(double));
    (data[i].ydat)=malloc((grafica.data[i].xnlin)*sizeof(double));
/*     printf(" Antes de trudata \n"); */
    TrueData(data[i].xdat,data[i].ydat,grafica.data[i],grafica.operac[i],&(data[i].ndat));    
    i++;
  }
  
  if((grafica.limits.xmin == 0. && grafica.limits.xmax == 0.)) {
    xmin=1e30;
    xmax=-1e30;
    i=0;
    
    while(grafica.data[i].status && i<NBUF){
      for(j=0;j<data[i].ndat;j++) {       
	if(xmax<data[i].xdat[j]) xmax=data[i].xdat[j];
	if(xmin>data[i].xdat[j]) xmin=data[i].xdat[j];
      }
      i++;
    }
    ran=xmax-xmin;
    xmax=xmax + ran*.05;
    xmin=xmin - ran*.05;
    
  }
  else {
    xmin=grafica.limits.xmin;
    xmax=grafica.limits.xmax;
  }
  if((grafica.limits.ymin == 0. && grafica.limits.ymax == 0.)) {
    ymin=1e30;
    ymax=-1e30;
    i=0;
    
    while(grafica.data[i].status && i<NBUF){
      for(j=0;j<data[i].ndat;j++) {       
	if(ymax<data[i].ydat[j]) ymax=data[i].ydat[j];
	if(ymin>data[i].ydat[j]) ymin=data[i].ydat[j];
      }
      i++;
    }
    ran=ymax-ymin;
    ymax=ymax + ran*.05;
    ymin=ymin - ran*.05;
  }
  else {
    ymin=grafica.limits.ymin;
    ymax=grafica.limits.ymax;
  }
  
  cpgslw(2);
  cpgscf(grafica.labels.font);
  cpgsch(grafica.labels.height);
  cpgvstd();
  if(grafica.limits.ratio==1) cpgwnad(xmin,xmax,ymin,ymax);
  else  cpgswin(xmin,xmax,ymin,ymax);

  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
  cpglab(grafica.labels.xlabel,grafica.labels.ylabel,grafica.labels.title);
  i=0;
  printf("\n X limits: %f %f\n",xmin,xmax);
  printf(" Y limits: %f %f\n",ymin,ymax);
  xley=grafica.leyend.x+(xmax-xmin)/20;
  while(grafica.data[i].status && i<NBUF) {
    printf(" Buffer %d \n",i);
    
    cpgsci(grafica.symbol[i].color);
    cpgsch(grafica.symbol[i].height);
    cpgpt_d(data[i].ndat,(data[i].xdat),(data[i].ydat),grafica.symbol[i].marker);
    if(grafica.symbol[i].lw!=0) {
      cpgslw(grafica.symbol[i].lw);
      cpgsls(grafica.symbol[i].ls);
      cpgline_d(data[i].ndat,data[i].xdat,data[i].ydat);
    }
    cpgsls(1);
    cpgsci(1);
    cpgslw(2);
/*     printf(" LA LEYENDA ES %d\n",grafica.leyend.active); */
    if(grafica.leyend.active) {
/*       printf(" HACE LA LEYENDA\n"); */
      cpgscf(grafica.leyend.font);
      cpgsch(grafica.leyend.height);
      cpglen(4,"H",&nul,&yleyr);
/*       printf(" LAponae  %f %f\n",grafica.leyend.x,grafica.leyend.y); */
      cpgsci(1);
      yley=grafica.leyend.y-yleyr*1.8*i;
      cpgptxt(xley,yley-yleyr/2,0.0,0.0,grafica.data[i].label);
      cpgsci(grafica.symbol[i].color);
      
      cpgpt1(grafica.leyend.x+(xmax-xmin)/80,yley,grafica.symbol[i].marker);
      if(grafica.symbol[i].lw!=0) {
	cpgslw(grafica.symbol[i].lw);
	cpgsls(grafica.symbol[i].ls);
	cpgmove(grafica.leyend.x-(xmax-xmin)/80,yley);
	cpgdraw(grafica.leyend.x+(xmax-xmin)/30,yley);
      }
      cpgsls(1);

    }
    cpgsci(1);

    printf("   File: %s   X column: %d   Y column: %d\n",grafica.data[i].file,grafica.data[i].xcol,grafica.data[i].ycol);
    printf("   Number of data: %d  Color: %d  Marker %d\n",data[i].ndat,grafica.symbol[i].color,grafica.symbol[i].marker);
    

    if(grafica.lines[i].status && i<NBUF) {
      cpgslw(grafica.lines[i].width);
      cpgsci(grafica.lines[i].color);
      cpgmove(xmin,grafica.lines[i].b+grafica.lines[i].a*xmin);
      cpgdraw(xmax,grafica.lines[i].b+grafica.lines[i].a*xmax);
      cpgsci(1);
      cpgslw(2);
    }
    printf("   Number of data points    : %d\n",data[i].ndat);
    free(data[i].xdat);
    free(data[i].ydat);
    i++;
    
  }
}


void ChangeLabels(struct graph *grafica)
{
/*   //  char nul[100]; */
  printf("Input X label ");
  reads((*grafica).labels.xlabel,(*grafica).labels.xlabel);
/*   //  scanf("%[",(*grafica).labels.xlabel); */
  printf("Input Y label ");
  reads((*grafica).labels.ylabel,(*grafica).labels.ylabel);
/*   //  fgets((*grafica).labels.ylabel,100,stdin); */
/*   //  scanf("%[",(*grafica).labels.ylabel); */
  printf("Input graph title ");
  reads((*grafica).labels.title,(*grafica).labels.title);
/*   //  fgets((*grafica).labels.title,100,stdin); */
/*   //  scanf("%[",(*grafica).labels.title); */
  printf("Input character font 1-4  ");
  (*grafica).labels.font=readi((*grafica).labels.font);
  /* scanf("%d",&(*grafica).labels.font); */

  printf("Input character height  ");
  (*grafica).labels.height=readf((*grafica).labels.height);

/*   //scanf("%f",&(*grafica).labels.height); */
/*   //Asi lo hacia antes: */
/*   //  fflush(stdin); */
/*   //gets(nul); */
/*   //strcpy((*grafica).labels.title,nul); */
}
void ChangeLeyend(struct graph *grafica)
{
/*   char nul[100]; */
  char opt;
  int i;
  printf("Active leyend (Y/N):");
  opt=readc('Y');
/*   //scanf("%s",opt); */
/*   //getline(opt,100); */
  if (opt=='N')    (*grafica).leyend.active=0;
  else {
    (*grafica).leyend.active=1;
    printf("Input X position:");
/*     //getline((*grafica).leyend.x,100); */
    (*grafica).leyend.x=readf((*grafica).leyend.x);
/*     //scanf("%f",&(*grafica).leyend.x); */
    printf("Input Y position:");
/*     //getline((*grafica).leyend.y,100); */
    (*grafica).leyend.y=readf((*grafica).leyend.y);
/*     //scanf("%f",&(*grafica).leyend.y); */
    printf("Input Font:");
    (*grafica).leyend.font=readi((*grafica).leyend.font);
/*     //scanf("%d",&(*grafica).leyend.font ); */
    printf("Input character height:");
    (*grafica).leyend.height=readf((*grafica).leyend.height);
/*     //scanf("%f",&(*grafica).leyend.height); */
    for(i=0;i<NBUF;i++) {
      if ((*grafica).data[i].status) {
	printf("Input label for buffer %d ",i);
	fflush(stdin);
	setvbuf(stdin,"",_IOLBF,0);
	setvbuf(stdout,"",_IOLBF,0);
	reads((*grafica).data[i].label,(*grafica).data[i].label);
/* 	//scanf("%s",(*grafica).data[i].label); */

/* 	//fflush(stdin); */
/* 	//gets(nul); */
/* 	//strcpy((*grafica).data[i].label,nul); */
      }
    }
  }
}
void ChangeLimits(struct graph *grafica)
{
/*   //printf("Input xmin [%f]: ",(*grafica).limits.xmin); */
  printf("Input xmin : ");
  (*grafica).limits.xmin=readf((*grafica).limits.xmin);
/*   //scanf("%f",&(*grafica).limits.xmin); */
/*   //printf("Input xmax [%f]: ",(*grafica).limits.xmax); */
  printf("Input xmax : ");
  (*grafica).limits.xmax=readf((*grafica).limits.xmax);
/*   //scanf("%f",&(*grafica).limits.xmax); */
/*   //printf("Input ymin [%f]: ",(*grafica).limits.ymin); */
  printf("Input ymin : ");
  (*grafica).limits.ymin=readf((*grafica).limits.ymin);
/*   //scanf("%f",&(*grafica).limits.ymin); */
/*   //printf("Input ymax [%f]: ",(*grafica).limits.ymax); */
  printf("Input ymax : ");
  (*grafica).limits.ymax=readf((*grafica).limits.ymax);
/*   //scanf("%f",&(*grafica).limits.ymax); */
/*   //printf("Input aspect ratio (0=not fixed) [%f]: ",(*grafica).limits.ratio); */
  printf("Input aspect ratio (0=not fixed) : ");
  (*grafica).limits.ratio=readf((*grafica).limits.ratio);
/*   //scanf("%f",&(*grafica).limits.ratio); */
}

void ChangeLimitsMouse(struct graph *grafica)
{

  char cnul;

  cpgsci(2);
  printf(" Press bottom left square with mouse...\n");
  cpgcurs_d(&(*grafica).limits.xmin,&(*grafica).limits.ymin,&cnul);
  printf(" Press top right square with mouse...\n");
  cpgband_d(2,1,(*grafica).limits.xmin,(*grafica).limits.ymin,&(*grafica).limits.xmax,&(*grafica).limits.ymax,&cnul);
  cpgsci(1);
}






void ChangeMarkers(struct graph *grafica)
{
  int buffer=1;
/*   //printf("Input buffer 1-%d [1] :",NBUF); */
  printf("Input buffer 1-%d :",NBUF);
  buffer=readi(buffer);
/*   //printf(" Cacita %d\n",buffer); */
/*   //scanf("%d",&buffer); */
  buffer--;
/*   printf("Input marker 0-31 [%d] :",(*grafica).symbol[buffer].marker); */
/*   scanf("%d",&(*grafica).symbol[buffer].marker); */
/*   printf("Input color 0-15 [%d] :",(*grafica).symbol[buffer].color); */
/*   scanf("%d",&(*grafica).symbol[buffer].color); */
/*   printf("Input height [%f] :",(*grafica).symbol[buffer].height); */
/*   scanf("%f",&(*grafica).symbol[buffer].height); */
/*   printf("Input line width for connecting points 1-201 (0=none)  [%d] :",(*grafica).symbol[buffer].lw); */
/*   scanf("%d",&(*grafica).symbol[buffer].lw); */
  printf("Input marker 0-31 :");
  (*grafica).symbol[buffer].marker=readi((*grafica).symbol[buffer].marker);
  printf("Input color 0-15  :");
  (*grafica).symbol[buffer].color=readi((*grafica).symbol[buffer].color);
  printf("Input height  :");
  (*grafica).symbol[buffer].height=readf((*grafica).symbol[buffer].height);
  printf("Input line width for connecting points 1-201 (0=none)  :");
  (*grafica).symbol[buffer].lw=readi((*grafica).symbol[buffer].lw);
  if((*grafica).symbol[buffer].lw!=0) {
/*     printf("Input line style (1 full, 2 dashed, 3 dot-dash-dot-dash, 4 dotted, 5 dash-dot-dot-dot)  [%d] :",(*grafica).symbol[buffer].ls); */
/*     scanf("%d",&(*grafica).symbol[buffer].ls); */
    printf("Input line style (1 full, 2 dashed, 3 dot-dash-dot-dash, 4 dotted, 5 dash-dot-dot-dot) :");
    (*grafica).symbol[buffer].ls=readi((*grafica).symbol[buffer].ls);
  }


}


void SaveParam(struct graph *graficas)
{
  FILE *fp;
  int i,j;
  int nc=0,nt;
/*   char nul[2]; */
/*   char file[100]; */
  char ch21[51];
/*   char opt; */
  printf("Name of file (%s) :",parfilename);
  reads(parfilename,parfilename);
  if((fp=fopen(parfilename,"w")) ==NULL) {
    printf("ERROR: Can't open file\n");
    return;
  }
  fprintf(fp,"COMMENT  Parameter file for GraphTable                                         \n");
  nc++;
  for(i=0;i<NGRAPH;i++) {
    
    for(j=0;j<NBUF;j++) {
      if(graficas[i].data[j].status) {
	sprintf(ch21,"'%s'",graficas[i].data[j].file);
	fprintf(fp,"FILE%1d_%1d = %-51.51s / File G. %d, B. %d\n",i,j,ch21,i,j);
	fprintf(fp,"XCOL%1d_%1d =%21d / X column for graph %d, buffer %d                \n",i,j,graficas[i].data[j].xcol,i,j);
	fprintf(fp,"YCOL%1d_%1d =%21d / Y column  for graph %d, buffer %d               \n",i,j,graficas[i].data[j].ycol,i,j);
	fprintf(fp,"IDCOL%1d_%1d=%21d / ID column  for graph %d, buffer %d              \n",i,j,graficas[i].data[j].idcol,i,j);
	fprintf(fp,"MARK%1d_%1d =%21d / Marker for graph %d, buffer %d                  \n",i,j,graficas[i].symbol[j].marker,i,j);
	fprintf(fp,"HEIG%1d_%1d =%21f / Height marker for graph %d, buffer %d           \n",i,j,graficas[i].symbol[j].height,i,j);
	fprintf(fp,"COLOR%1d_%1d=%21d / Color marker for graph %d, buffer %d            \n",i,j,graficas[i].symbol[j].color,i,j);
	fprintf(fp,"LW%1d_%1d   =%21d / Line width for graph %d, buffer %d              \n",i,j,graficas[i].symbol[j].lw,i,j);
	fprintf(fp,"LS%1d_%1d   =%21d / Line style for graph %d, buffer %d              \n",i,j,graficas[i].symbol[j].ls,i,j);
	fprintf(fp,"XLOG%1d_%1d =%21c / Apply logs in X for graph %d, buffer %d (Y/N)   \n",i,j,graficas[i].operac[j].xlog,i,j);
	fprintf(fp,"YLOG%1d_%1d =%21c / Apply logs in Y for graph %d, buffer %d (Y/N)   \n",i,j,graficas[i].operac[j].ylog,i,j);
	fprintf(fp,"OPT0%1d_%1d =%21c / Aply log for values <=0 for graph %d, buffer %d \n",i,j,graficas[i].operac[j].opt0,i,j);
	nc+=12;
	if(graficas[i].leyend.active) {
	  sprintf(ch21,"'%s'",graficas[i].data[j].label);
	  fprintf(fp,"LEYT%1d_%1d = %-51.51s / Leyend G %d, B %d\n",i,j,ch21,i,j);  
	  nc++;
	}
      }
      if(graficas[i].lines[j].status) {
	fprintf(fp,"LIN_A%1d_%1d=%21f / A data for y=ax+b  graph %d, buffer %d          \n",i,j,graficas[i].lines[j].a,i,j);
	fprintf(fp,"LIN_B%1d_%1d=%21f / B data for y=ax+b  graph %d, buffer %d          \n",i,j,graficas[i].lines[j].b,i,j);
	fprintf(fp,"LIN_C%1d_%1d=%21d / Color line for     graph %d, buffer %d          \n",i,j,graficas[i].lines[j].color,i,j);
	fprintf(fp,"LIN_W%1d_%1d=%21d / Width line for     graph %d, buffer %d          \n",i,j,graficas[i].lines[j].width,i,j);
	nc+=4;
      }
    }
    if(graficas[i].leyend.active) {
      fprintf(fp,"XLEYEN_%1d=%21f / X position for leyend in graph %d              \n",i,graficas[i].leyend.x,i);
      fprintf(fp,"YLEYEN_%1d=%21f / Y position for leyend in graph %d              \n",i,graficas[i].leyend.y,i);
      fprintf(fp,"LEYFON_%1d=%21d / Font for leyend in graph %d                    \n",i,graficas[i].leyend.font,i);
      fprintf(fp,"LEYHEI_%1d=%21f / Height for labels in graph %d                  \n",i,graficas[i].leyend.height,i);
      nc+=4;
      
    }
    
    sprintf(ch21,"'%s'",graficas[i].labels.xlabel);
/*     printf(" garbo el label %d %s\n",i,ch21); */
    fprintf(fp,"XLABEL_%1d= %-51.51s / X Label graph %d\n",i,ch21,i);
    sprintf(ch21,"'%s'",graficas[i].labels.ylabel);
    fprintf(fp,"YLABEL_%1d= %-51.51s / Y Label graph %d\n",i,ch21,i);
    sprintf(ch21,"'%s'",graficas[i].labels.title);
    fprintf(fp,"TITLE__%1d= %-51.51s / X Label graph %d\n",i,ch21,i);
    fprintf(fp,"LABFON_%1d=%21d / Font for labels in graph %d                    \n",i,graficas[i].labels.font,i);
    fprintf(fp,"LABHEI_%1d=%21f / Height for labels in graph %d                  \n",i,graficas[i].labels.height,i);
    fprintf(fp,"XMIN_%1d  =%21f / Xmin in graph %d                               \n",i,graficas[i].limits.xmin,i);
    fprintf(fp,"XMAX_%1d  =%21f / Xmax in graph %d                               \n",i,graficas[i].limits.xmax,i);
    fprintf(fp,"YMIN_%1d  =%21f / Ymin in graph %d                               \n",i,graficas[i].limits.ymin,i);
    fprintf(fp,"YMAX_%1d  =%21f / Ymax in graph %d                               \n",i,graficas[i].limits.ymax,i);
    fprintf(fp,"RATI_%1d  =%21f / Ratio aspect for axes in graph %d              \n",i,graficas[i].limits.ratio,i);
    nc+=10;

  }

  fprintf(fp,"COMMENT  Blank lines to assure the last block of 2880 caracters is complete    \n");
  fprintf(fp,"COMMENT                                                                        \n");
  nc+=2;
  
  for(nt=nc;nt<(1+(int)(nc/36))*36-1;nt++) {
/*     //printf("NC %d %d\n",nt,(1+(int)(nc/36.))*36); */
    fprintf(fp,"COMMENT                                                                        \n");
  }
  fprintf(fp,"END                                                                            \n");
  fclose(fp);
}

void LoadParam(struct graph *graficas, FILE *filepar, char file[50])
{
  int i,j;
  int nlin;
  char nul3[3],nul1[1];
  char keyf[9]="",key[9]="";
  int status=0;
/*   char string[51]; */
  char comment[51];
/*   fitsfile *fits; */
  fitsfile *parfile;
/*   int inte; */
/*   float crpix; */
  printf("Using parameter file <<%s>>\n",file);
  /*   printf(" Y esto antes\n"); */
  if( ffopen2(&parfile,file, READONLY, &status)) fits_report_error(stderr,status);
  /*   if( ffopen(&parfile,file, READONLY, &status)) fits_report_error(stderr,status); */
  /*   printf(" Esto es despues\n"); */
  
  
  /* Seguro que con algo asi tambien se puede abrir. Lo he sacado de ffopen */  
  /*   if (ffopenfile((*fptr)->filename, 0, readwrite, &diskfile,  */
  /* 		 &filesize, status) > 0) */
  /*     { */
  /*       free((*fptr)->filename); */
  /*       free(*fptr); */
  /*       *fptr = 0;    */          /* return null file pointer */ 
  /*       ffpmsg("ffopen failed to find and/or open the following file:"); */
  /*       ffpmsg(&filename[ii]); */
  /*       return(*status = FILE_NOT_OPENED); */
  /*     } */
  /*         (*fptr)->fileptr = diskfile;  */ /* store the file pointer */ 
  /*     } */
  /*     (*fptr)->filesize = filesize;       */   /* physical file size */ 
  /*     (*fptr)->logfilesize = filesize;    */   /* logical file size */ 
  /*     (*fptr)->writemode = readwrite;     */   /* read-write mode    */ 
  /*     (*fptr)->datastart = DATA_UNDEFINED; */  /* unknown start of data */ 
  /*     (*fptr)->curbuf = -1; */  /* undefined current IO buffer */ 
    



  for(i=0;i<NGRAPH;i++) {
/*     printf(" GAPG %d\n",i); */
    status=0;
    
    sprintf(nul1,"%1d",i);
    strcpy(key,"XLABEL_");
    strcat(key,nul1);
/*     printf(" BEFORE\n");  */
    ffgky(parfile,TSTRING,key,graficas[i].labels.xlabel,comment,&status);
    if(status) printf(" keyword not found: %s\n",key);
    fits_report_error(stderr,status);
    status=0;
/*     printf("%s>>  %s \n",key,graficas[i].labels.xlabel);  */
      
      
    /*       //    if(ffgky(parfile,TSTRING,key,graficas[i].labels.xlabel); */
    strcpy(key,"YLABEL_");
    strcat(key,nul1);
    
    ffgky(parfile,TSTRING,key,graficas[i].labels.ylabel,comment,&status);
    /*       //printf("%s>>  %s \n",key,graficas[i].labels.ylabel); */
    if(status) printf(" keyword not found: %s\n",key);
    fits_report_error(stderr,status);
    status=0;
    strcpy(key,"TITLE__");
    strcat(key,nul1);
    ffgky(parfile,TSTRING,key,graficas[i].labels.title,comment,&status);
/*     printf("%s>>  %s \n",key,graficas[i].labels.title);  */
    if(status) printf(" keyword not found: %s\n",key);
    fits_report_error(stderr,status);
    status=0;    
    strcpy(key,"LABFON_");
    strcat(key,nul1);
    
    ffgky(parfile,TINT,key,&(graficas[i].labels.font),comment,&status);
    if(status) printf(" keyword not found: %s\n",key);
    fits_report_error(stderr,status);
    status=0;
    
/*     printf("%s>>  %d \n",key,graficas[i].labels.font);  */
      
      
      
      
    strcpy(key,"LABHEI_");
    strcat(key,nul1);
    
    ffgky(parfile,TFLOAT,key,&(graficas[i].labels.height),comment,&status);  
    if(status) printf(" keyword not found: %s\n",key);
    status=0;
    

    strcpy(key,"XMIN_");
    strcat(key,nul1);
    ffgky(parfile,TDOUBLE,key,&(graficas[i].limits.xmin),comment,&status);  
    if(status) printf(" keyword not found: %s\n",key);
    status=0;
    strcpy(key,"XMAX_");
    strcat(key,nul1);
    ffgky(parfile,TDOUBLE,key,&(graficas[i].limits.xmax),comment,&status);  
    if(status) printf(" keyword not found: %s\n",key);
    status=0;
    strcpy(key,"YMIN_");
    strcat(key,nul1);
    ffgky(parfile,TDOUBLE,key,&(graficas[i].limits.ymin),comment,&status);  
    if(status) printf(" keyword not found: %s\n",key);
    status=0;
    strcpy(key,"YMAX_");
    strcat(key,nul1);
    ffgky(parfile,TDOUBLE,key,&(graficas[i].limits.ymax),comment,&status);  
    if(status) printf(" keyword not found: %s\n",key);
    status=0;
    strcpy(key,"RATI_");
    strcat(key,nul1);
    ffgky(parfile,TDOUBLE,key,&(graficas[i].limits.ratio),comment,&status);  
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
/*     printf(" STATSS %d\n",status);  */
    status=0;
    for(j=0;j<NBUF;j++) {
      strcpy(keyf,"FILE");
      strcpy(nul3,"");
      sprintf(nul3,"%1d_%1d",i,j);
      strcat(keyf,nul3);
/*       //      printf(" Searching keyword  <<%s>>\n",keyf); */
      status=0;
      if(ffgky(parfile,TSTRING,keyf,graficas[i].data[j].file,comment,&status)!=0) {
	graficas[i].data[j].status=0;
	graficas[i].lines[j].status=0;
      }
      else {
/* 	//printf(" Keyword <<%s>>. File %d to read: <<%s>>\n",keyf,j,graficas[i].data[j].file); */
	graficas[i].status=1;
	graficas[i].data[j].status=1;
	strcpy(key,"XCOL");
	strcat(key,nul3);
	ffgky(parfile,TINT,key,&(graficas[i].data[j].xcol),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;
	strcpy(key,"YCOL");
	strcat(key,nul3);
	ffgky(parfile,TINT,key,&(graficas[i].data[j].ycol),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	if(status) {
	  printf("ERROR: No column data for file %s STOP",graficas[i].data[j].file);
	  exit(1);
	}
	strcpy(key,"IDCOL");
	strcat(key,nul3);
	ffgky(parfile,TINT,key,&(graficas[i].data[j].idcol),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);

	status=0;
	strcpy(key,"MARK");
	strcat(key,nul3);
	ffgky(parfile,TINT,key,&(graficas[i].symbol[j].marker),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;

	strcpy(key,"COLOR");
	strcat(key,nul3);
	ffgky(parfile,TINT,key,&(graficas[i].symbol[j].color),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;

	strcpy(key,"LW");
	strcat(key,nul3);
	ffgky(parfile,TINT,key,&(graficas[i].symbol[j].lw),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;

	strcpy(key,"LS");
	strcat(key,nul3);
	ffgky(parfile,TINT,key,&(graficas[i].symbol[j].ls),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;

	strcpy(key,"HEIG");
	strcat(key,nul3);
	ffgky(parfile,TFLOAT,key,&(graficas[i].symbol[j].height),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;

	strcpy(key,"XLOG");
	strcat(key,nul3);
	ffgky(parfile,TSTRING,key,&(graficas[i].operac[j].xlog),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	fits_report_error(stderr,status);
	status=0;
	
/* 	printf("%s>>  <<%c>> \n",key,graficas[i].operac[j].xlog);  */



	strcpy(key,"YLOG");
	strcat(key,nul3);
	ffgky(parfile,TSTRING,key,&(graficas[i].operac[j].ylog),comment,&status);
	strcpy(key,"OPT0");
	strcat(key,nul3);
	ffgky(parfile,TSTRING,key,&(graficas[i].operac[j].opt0),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;
	
	if(graficas[i].leyend.active) {
	  strcpy(key,"LEYT");
	  strcat(key,nul3);
	  ffgky(parfile,TSTRING,key,graficas[i].data[j].label,comment,&status);
	  if(status) printf(" keyword not found: %s\n",key);
	}
	
/* 	//	printf("%s>>  <<%c>> \n",key,graficas[i].operac[j].opt0); */
	
	nlin=FileNLin(graficas[i].data[j].file);
	
	graficas[i].data[j].xdata=malloc(nlin*sizeof(double));
	graficas[i].data[j].ydata=malloc(nlin*sizeof(double));
	graficas[i].data[j].id=malloc(NSTR*nlin*sizeof(char));
	graficas[i].data[j].xtest=malloc(nlin*sizeof(int));
	graficas[i].data[j].ytest=malloc(nlin*sizeof(int));
	graficas[i].data[j].idtest=malloc(nlin*sizeof(int));
	

	printf("Reading column %d from file %s\n",graficas[i].data[j].xcol,graficas[i].data[j].file);
	ReadDoublecol(graficas[i].data[j].file,graficas[i].data[j].xcol,
		   graficas[i].data[j].xdata,graficas[i].data[j].xtest,
		   &(graficas[i].data[j].xnlin));

	printf("Reading column %d from file %s\n",graficas[i].data[j].ycol,graficas[i].data[j].file);
	
	ReadDoublecol(graficas[i].data[j].file,graficas[i].data[j].ycol,
		   graficas[i].data[j].ydata,graficas[i].data[j].ytest,
		      &(graficas[i].data[j].ynlin));
	if(graficas[i].data[j].idcol!=0) ReadCharcol((graficas[i]).data[j].file,(graficas[i]).data[j].idcol,
						       (graficas[i]).data[j].id,(graficas[i]).data[j].idtest,NSTR,
						       &(graficas[i]).data[j].idnlin);
	
	if(graficas[i].data[j].xnlin!=graficas[i].data[j].ynlin) {
	  printf("ERROR: x dat and y data hasn't the same size");
	  graficas[i].data[j].status=0;
	}	
      }
      strcpy(key,"LIN_A");
      strcat(key,nul3);
      if(ffgky(parfile,TFLOAT,key,&(graficas[i].lines[j].a),comment,&status)==0) {
	if(status) printf(" keyword not found: %s\n",key);
	graficas[i].status=1;
	strcpy(key,"LIN_B");
	strcat(key,nul3);
	ffgky(parfile,TFLOAT,key,&(graficas[i].lines[j].b),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;
        strcpy(key,"LIN_C");
        strcat(key,nul3);
        ffgky(parfile,TINT,key,&(graficas[i].lines[j].color),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;

        strcpy(key,"LIN_W");
        strcat(key,nul3);
        ffgky(parfile,TINT,key,&(graficas[i].lines[j].width),comment,&status);
	if(status) printf(" keyword not found: %s\n",key);
	status=0;

	graficas[i].lines[j].status=1;
      }
      
    }
  }
  fits_close_file(parfile,&status);
  

}

void  PerformOper(struct graph *graficas)
{
  char opt;
  int buffer=1;

  fflush(stdin);
  setvbuf(stdin,"",_IOLBF,0);
  setvbuf(stdout,"",_IOLBF,0);
  
/*   printf("Input buffer 1-%d [1] :",NBUF); */
/*   scanf("%d",&buffer); */
  printf("Input buffer 1-%d:",NBUF);
  buffer=readi(buffer);
/*   //printf(" CASSA %d\n",buffer); */
  buffer--;
/*   //printf(" A ver si est abien 333333 %s\n",(*graficas).data[buffer].file); */

  printf("Do you want to use logaritms in X data Y/N :");
  fflush(stdin);
  setvbuf(stdin,"",_IOLBF,0);
  setvbuf(stdout,"",_IOLBF,0);
  
  opt=readc('N');
  if(opt=='Y') (*graficas).operac[buffer].xlog='Y';
  else (*graficas).operac[buffer].xlog='N';
  printf("Do you want to use logaritms in Y data Y/N :");
  fflush(stdin);
  setvbuf(stdin,"",_IOLBF,0);
  setvbuf(stdout,"",_IOLBF,0);
  
  opt=readc('N');
  if(opt=='Y') (*graficas).operac[buffer].ylog='Y';
  else (*graficas).operac[buffer].ylog='N';
  printf("Do you want to apply logaritms to values <=0 Y/N :");
  fflush(stdin);
  setvbuf(stdin,"",_IOLBF,0);
  setvbuf(stdout,"",_IOLBF,0);
  
  opt=readc('N');
  if(opt=='Y') (*graficas).operac[buffer].opt0='Y';
  else (*graficas).operac[buffer].opt0='N';


  printf("Enter constant to sum to X data: ");
  (*graficas).operac[buffer].xplus=readf(0.);
  printf("Enter constant to sum to Y data: ");
  (*graficas).operac[buffer].yplus=readf(0.);
/*   //scanf("%f",&(*graficas).operac[buffer].xplus); */
/*   //printf(" A ver si SIGUE 4444 abien %s\n",(*graficas).data[buffer].file); */

}

void DrawLines(struct graph *graficas)
{
  int buffer=1;
  printf("Input buffer 1-%d:",NBUF);
  buffer=readi(buffer);
/*   //  scanf("%d",&buffer); */
  buffer--;
  (*graficas).lines[buffer].status=1;
  printf("Line y = ax +b\n");
  printf("Input a: ");
  (*graficas).lines[buffer].a=readf((*graficas).lines[buffer].a);
/*   //scanf("%f",&(*graficas).lines[buffer].a); */
  printf("Input b: ");
  (*graficas).lines[buffer].b=readf((*graficas).lines[buffer].b);
/*   //scanf("%f",&(*graficas).lines[buffer].b); */
  
  printf("Input color: ");
  (*graficas).lines[buffer].color=readi((*graficas).lines[buffer].color);
/*   //scanf("%d",&(*graficas).lines[buffer].color); */
  printf("Input width: ");
  (*graficas).lines[buffer].width=readi((*graficas).lines[buffer].width);
/*   //scanf("%d",&(*graficas).lines[buffer].width); */
}

void TrueData(double *xdata,double * ydata,struct dat data,struct opr oper,int *ndata)
{
  int n=0;
  int j;
/*   //printf("Inici a True Dat\n"); */
  for(j=0;j<data.xnlin;j++) {
    if(data.xtest[j] && data.ytest[j] ) {
      xdata[n]=data.xdata[j]+oper.xplus;
      ydata[n]=data.ydata[j]+oper.yplus;
      if(oper.ylog=='Y') {
/* 	printf(" YES\n"); */
	if(ydata[n]<=0) {
	  if(oper.opt0=='Y') ydata[n]=0;
	  else continue;
	}
	else {
/* 	  //printf(" daots : %f log %f\n",ydata[n],log10(ydata[n])); */
	  ydata[n]=log10(ydata[n]);
	}
      }
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
    double *xdat;
    double *ydat;	
    int ndat;
  };
  double x,y;
  char ch='A';
  double rmin,r;
  int i,j;
  int ni,nj,njj;
  struct data data[NBUF];
  while(ch=='A')
    { 
      cpgcurs_d(&x,&y,&ch);
      ni=0;
      nj=0;
      rmin=1e30;      
      i=0;
      while(grafica.data[i].status) {
	
	(data[i].xdat)=malloc((grafica.data[i].xnlin)*sizeof(double));
	(data[i].ydat)=malloc((grafica.data[i].xnlin)*sizeof(double));
	TrueData(data[i].xdat,data[i].ydat,grafica.data[i],grafica.operac[i],&(data[i].ndat));    
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
      printf("Y cursor: %f\n\n",y);  
      printf("X data: %f\n",data[ni].xdat[nj]);
      printf("Y data: %f\n\n",data[ni].ydat[nj]);
      printf("X column: %e\n",grafica.data[ni].xdata[njj]);
      printf("Y column: %e\n",grafica.data[ni].ydata[njj]);
      if(grafica.data[ni].idcol!=0) printf("ID column: %s\n",grafica.data[ni].id+njj*NSTR);
    }
}


int NumberData(int nj,struct dat data,struct opr oper)
{
  int n=-1;
  int j;
  double  xtemp,ytemp;
  for(j=0;j<data.xnlin;j++) {
    if(data.xtest[j] && data.ytest[j] ) {
      xtemp=data.xdata[j]+oper.xplus;
      ytemp=data.ydata[j]+oper.yplus;
      if(oper.ylog=='Y') {
	if(ytemp<=0) {
	  if(oper.opt0=='Y') ;
	  else continue;
	}
      }
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



