#include "modulos.h"
         
   
struct gaussdist { 
  double mean;
  double errmean;
  double sigma;
  double errsigma;
  double covmeansigma;
};   
 
struct histodist {
  int nbin;
  double *p;
  double **covp;
  double *xbin;
  double xmin;
  double xmax; 
}; 
 
struct histodist2D {
/*Creado para tener un histograma de dos dimensiones. Terminado.*/
  int nbinx;
  int nbiny;
  double *p;
  double **covp;
  double *xbin;
  double *ybin;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
};

struct histodist2DFF {
/*Creado para tener un histograma de dos dimensiones en el caso del cálculo de la fracción de fusiones*/
  int nbinx;
  double *p;
  double **covp;
  double *xbin;
  double xmin;
  double xmax;
  double yff;
  double yinter;
};   
              
          
struct sample {
  unsigned int ndata;
  double *xdata;
  double *errdata;
};

struct sample2D {
/*Creado para tener en cuenta las dos variables de entrada.*/
  unsigned int ndata;
  double *xdata;
  double *errxdata;
  double *ydata;
  double *errydata;
};

                                 
 
void ReadSample(); 
void ReadSample2D();     
void Compute_g_g();
void Compute_h_g();
void Compute_hh_gg();
void Compute_ff_gg();
void Generate_g_g();
void Generate_h_g();
void PlotHistogram(int nbin,float *xmin,float *xmax);
 
struct sample Data;
struct sample2D Data2D;
 
int main() { 
/* Modificado para tener en cuenta las dos variables de entrada. Terminado.*/
  
  char opt='E';
  Data.ndata = 0;
  Data2D.ndata = 0;
 
  /* Se añade la opción D para leer los datos 2D y la opción K para obtener el histograma de 2D*/
  /* Se añade la opción F para determinar directamente los valores de la fracción de fusiones */
  do{
    printf("\n R Read sample\n"); 
    printf(" D Read sample for 2D histograms\n"); 
    printf(" G Compute gaussian distribution assuming gaussian errors\n"); 
    printf(" H Compute step distribution assuming gaussian errors\n");
    printf(" K Compute step distribution 2D assuming gaussian errors\n");
    printf(" F Compute the merger fraction of the 2D sample assuming gaussian errors\n");
    printf(" E Exit\n");
    opt=readc(opt); 
    switch (opt) { 
    case 'R':
    case 'r':
      ReadSample();
      break;
    case 'D':
    case 'd':
      ReadSample2D();
      break;
    case 'G':
    case 'g':
      if (Data.ndata == 0) ReadSample();
      Compute_g_g();
      break;
    case 'H':
    case 'h':
      if (Data.ndata == 0) ReadSample();
      Compute_h_g();
      break;
    case 'K':
    case 'k':
      if (Data2D.ndata == 0) ReadSample2D(); 
      Compute_hh_gg();
      break;
    case 'F':
    case 'f':
      if (Data2D.ndata == 0) ReadSample2D(); 
      Compute_ff_gg();
      break;
      
    }
  }while(opt!='E' && opt!='e');
  return(0);
}


void Compute_g_g() {
  
  int i;
  float x;
  float xmin=0,xmax=0;

  static struct gaussdist GD;
  int mlstatus;
  
  static int nbin=20;

  
  printf(" Computing (this may take a while depending on sample size)...\n");
  mlstatus=ML_g_g_d(Data.ndata,Data.xdata,Data.errdata,&(GD.mean),&(GD.sigma),&(GD.errmean),&(GD.errsigma),&(GD.covmeansigma));

   
  if(mlstatus==0)     printf(" Computation exited normaly\n");
  if(mlstatus==2) {
    printf(" WARNING: Computation reached maximum iterations\n");
    printf("          Maybe the solution is not good. Try again\n");    
  }

  
  printf("\n\n Given that the initial distribution was a gaussian\n with parameters mean and sigma, the best estimation is:\n");
  printf("\n Mean  = %.10g +/- %.10g\n",GD.mean,GD.errmean);
  printf("\n Sigma = %.10g +/- %.10g\n",GD.sigma,GD.errsigma);
  
  printf(" Now a plot to show the distribution\n");
  printf(" Input number of bins in histogram: ");
  nbin=readi(nbin);
  printf(" Histogram shows the observed distribution\n");


  printf(" Line shows the estimated original distribution\n");
  cpgopen("?");
  PlotHistogram(nbin,&xmin,&xmax);

  cpgsci(1);
  cpgmove(xmin,0.);
  for(i=0;i<1000;i++) {
    x=xmin+(float)i*(xmax-xmin)/999.;
    cpgdraw(x,Data.ndata*(xmax-xmin)/nbin*gaussian(x,GD.mean,GD.sigma));
  }
  cpgsci(1);
  cpgsls(1);
  cpgclos();
  
}

void Compute_h_g() 
{
  
  int i;
  float x;
  float xmin=0,xmax=0;

  static struct histodist HD;
  int mlstatus;

  static int init=0;

  if(init==0) {
    MinMax_d(Data.ndata,Data.xdata,&(HD.xmin),&(HD.xmax));
    HD.nbin=10;
    init=1;
  }


  printf(" Input minimum x value for the distribution (no allowed values below this number for the original distribution): ");
  HD.xmin=(double)readf(HD.xmin);
  printf(" Input maximum x value for the distribution (no allowed values above this number for the original distribution): ");
  HD.xmax=(double)readf(HD.xmax);
  printf(" Input number of bins for the distribution in x range: ");
  HD.nbin=readi(HD.nbin);
   
  HD.p   =vector_d(HD.nbin);
  HD.xbin=vector_d(HD.nbin+1);
  HD.covp=matrix_d(HD.nbin,HD.nbin);

  for(i=0;i<HD.nbin+1;i++)   HD.xbin[i]=HD.xmin+(float)i*(HD.xmax-HD.xmin)/(HD.nbin);

  printf(" Computing (this may take a while depending on sample size)...\n");
  mlstatus=MLA_h_g_d(Data.ndata,Data.xdata,Data.errdata,HD.nbin,HD.xbin,HD.p,HD.covp);   


  if(mlstatus==0)     printf(" Computation exited normaly\n");
  if(mlstatus==2) {
    printf(" WARNING: Computation reached maximum iterations\n");
    printf("          Maybe the solution is not good. Try again\n");    
  }

  printf("\n\n Given that the initial distribution was a step distribution\n with probability Pk in each interval, the best estimation is:\n");

  for(i=0;i<HD.nbin;i++)     printf(" k %2d  Pk %-17.15f +/- %-17.15f (N=%6d) in interval   x = [%g-%g]\n",i,(HD).p[i],sqrt((HD).covp[i][i]),(int)((HD).p[i]*(float)Data.ndata*((HD).xbin[i+1]-(HD).xbin[i])),(HD).xbin[i],(HD).xbin[i+1]);
  printf("\n The same distribution as above but in histogram format, instead of probabilities is:\n");
  for(i=0;i<HD.nbin;i++)     printf(" k %2d  N %-17.15f +/- %-17.15f in interval   x = [%g-%g]\n",i,(HD).p[i]*(float)Data.ndata*((HD).xbin[i+1]-(HD).xbin[i]),sqrt((HD).covp[i][i])*(float)Data.ndata*((HD).xbin[i+1]-(HD).xbin[i]),(HD).xbin[i],(HD).xbin[i+1]);

  printf(" Now a plot to show the distribution\n");
  printf(" Histogram shows the observed distribution\n");


  
  cpgopen("?");
  xmin=HD.xmin;
  xmax=HD.xmax;
  PlotHistogram(HD.nbin,&xmin,&xmax);

  cpgsci(1);
  cpgmove(xmin,0.);
  for(i=0;i<5000;i++) {
    x=xmin+(float)i*(xmax-xmin)/4999.;
    cpgdraw(x,((float)Data.ndata*(xmax-xmin)/HD.nbin*(HD).p[((int)(HD.nbin*(x-HD.xmin)/(HD.xmax-HD.xmin)))])); 
  }
  cpgsci(1);
  cpgsls(1);
  cpgclos();

  free(HD.p);
  free(HD.xbin);
  free_matrix_d(HD.covp,HD.nbin,HD.nbin);

}

void Compute_hh_gg() 
{
/* Desarrollado para hacer el estudio en dos variables.*/

  int i,j, toplot=0;  
  float xmin=0,xmax=0;
  float ymin=0,ymax=0; 

  static struct histodist2D HD;
  int mlstatus;

  static int init=0; 

  if(init==0) {
    MinMax_d(Data2D.ndata,Data2D.xdata,&(HD.xmin),&(HD.xmax));
    MinMax_d(Data2D.ndata,Data2D.ydata,&(HD.ymin),&(HD.ymax)); 
    HD.nbinx=10;
    HD.nbiny=10;
    init=1;
  }

  printf(" Input minimum x value for the distribution (no allowed values below this number for the original distribution): ");
  HD.xmin=(double)readf(HD.xmin);
  printf(" Input maximum x value for the distribution (no allowed values above this number for the original distribution): ");
  HD.xmax=(double)readf(HD.xmax);
  printf(" Input number of bins in X for the distribution: ");
  HD.nbinx=readi(HD.nbinx);
  
  /*Añadido para obtener los valores de y*/
  printf(" Input minimum y value for the distribution (no allowed values below this number for the original distribution): ");
  HD.ymin=(double)readf(HD.ymin);
  printf(" Input maximum y value for the distribution (no allowed values above this number for the original distribution): ");
  HD.ymax=(double)readf(HD.ymax);
  printf(" Input number of bins in Y for the distribution: ");
  HD.nbiny=readi(HD.nbiny);
  
  HD.p   =vector_d(HD.nbinx*HD.nbiny);
  HD.xbin=vector_d(HD.nbinx+1);
  HD.ybin=vector_d(HD.nbiny+1);
  HD.covp=matrix_d(HD.nbinx*HD.nbiny,HD.nbinx*HD.nbiny);

  for(i=0;i<HD.nbinx+1;i++)
  {
    HD.xbin[i]=HD.xmin+(float)i*(HD.xmax-HD.xmin)/(HD.nbinx); 
  }
  for(i=0;i<HD.nbiny+1;i++) 
  {
    HD.ybin[i]=HD.ymin+(float)i*(HD.ymax-HD.ymin)/(HD.nbiny); 
  }   
       
  printf(" Computing (this may take a while depending on sample size)...\n");
  mlstatus=MLA_hh_gg_d(Data2D.ndata,Data2D.xdata,Data2D.errxdata,Data2D.ydata,Data2D.errydata,HD.nbinx,HD.nbiny,HD.xbin,HD.ybin,HD.p,HD.covp);
  /* Añadido ydata, errydata, xbin e ybin*/ 
    
  if(mlstatus==0)     printf(" Computation exited normaly\n");           
  if(mlstatus==2) {                 
    printf(" WARNING: ComputatioD.yff=(double)readf(HD.yff);n reached maximum iterations\n");   
    printf("          Maybe the solution is not good. Try again\n");            
  } 
 
  printf("\n\n Given that the initial distribution was a step distribution\n with probability Pk in each interval, the best estimation is:\n");
 
  /* Aqui van los resultados. Sale en forma de matriz kl... ¿estará todo ok? :-O */  
  
  for(j=0;j<HD.nbiny;j++)     
  {  
      for (i=0;i<HD.nbinx;i++)            
      {     
        printf(" element %3d%3d  Pkl %-17.15f +/- %-17.15f (N=%6d) in interval   x = [%g-%g]  y = [%g-%g] \n",i,j,(HD).p[j*(HD).nbinx+i],sqrt((HD).covp[j*(HD).nbinx+i][j*(HD).nbinx+i]),(int)((HD).p[j*(HD).nbinx+i]*(float)Data2D.ndata*(HD.xbin[i+1]-HD.xbin[i])*(HD.ybin[j+1]-HD.ybin[j])),HD.xbin[i],HD.xbin[i+1],HD.ybin[j],HD.ybin[j+1]); 
      }
  } 
  printf("\n The same distribution as above but in histogram format, instead of probabilities is:\n");
  for(j=0;j<HD.nbiny;j++) 
  {   
    for (i=0;i<HD.nbinx;i++) 
    {  
       printf(" element kl %2d%2d   N %-17.15f +/- %-17.15f in interval   x = [%g-%g] y = [%g-%g] \n",i,j,(HD).p[j*(HD).nbinx+i]*(float)Data2D.ndata*((HD).xbin[i+1]-(HD).xbin[i])*((HD).ybin[j+1]-(HD).ybin[j]),sqrt((HD).covp[j*(HD).nbinx+i][j*(HD).nbinx+i])*(float)Data2D.ndata*((HD).xbin[i+1]-(HD).xbin[i])*((HD).ybin[j+1]-(HD).ybin[j]),(HD).xbin[i],(HD).xbin[i+1],HD.ybin[j],HD.ybin[j+1]);
    }
  }

  free(HD.p);  
  free(HD.xbin);
  free(HD.ybin);
  free_matrix_d(HD.covp,HD.nbinx*HD.nbiny,HD.nbinx*HD.nbiny);   

}

void Compute_ff_gg() 
{
/* Desarrollado para hacer el estudio en dos variables y determinar directamente la fracción de fusiones*/

  int i,j, toplot=0;  
  float xmin=0,xmax=0;
  float ymin=0,ymax=0; 

  static struct histodist2DFF HD;
  int mlstatus;

  static int init=0; 

  if(init==0) {
    MinMax_d(Data2D.ndata,Data2D.xdata,&(HD.xmin),&(HD.xmax)); 
    HD.nbinx=10;
    init=1;
  }

  printf(" Input minimum x value for the distribution (no allowed values below this number for the original distribution): ");
  HD.xmin=(double)readf(HD.xmin);
  printf(" Input maximum x value for the distribution (no allowed values above this number for the original distribution): ");
  HD.xmax=(double)readf(HD.xmax);
  printf(" Input number of bins in X for the distribution: ");
  HD.nbinx=readi(HD.nbinx);
  
  /*
  printf(" Input minimum y value for the distribution (no allowed values below this number for the original distribution): ");
  HD.ymin=(double)readf(HD.ymin);
  printf(" Input maximum y value for the distribution (no allowed values above this number for the original distribution): ");
  HD.ymax=(double)readf(HD.ymax);
  No me interesa, ya que solo es importante el límite para fusión y el tamaño del intervalo en asimetría*/
  
  printf(" Input the limit in A that define a mayor merger: ");
  HD.yff=(double)readf(HD.yff);
  printf(" Input the size of the two A intervals: ");
  HD.yinter=(double)readf(HD.yinter);
  
  HD.p   =vector_d(HD.nbinx*2);
  HD.xbin=vector_d(HD.nbinx+1); 
  HD.covp=matrix_d(HD.nbinx*2,HD.nbinx*2);

  for(i=0;i<HD.nbinx+1;i++) 
  {
    HD.xbin[i]=HD.xmin+(float)i*(HD.xmax-HD.xmin)/(HD.nbinx); 
  } 
       
  printf(" Computing (this may take a while depending on sample size)...\n");
  mlstatus=MLA_ff_gg_d(Data2D.ndata,Data2D.xdata,Data2D.errxdata,Data2D.ydata,Data2D.errydata,HD.nbinx,HD.xbin,HD.yff,HD.yinter,HD.p,HD.covp);
  /* Añadido ydata, errydata, yff y yinter*/
 
  if(mlstatus==0)     printf(" Computation exited normaly\n");   
  if(mlstatus==2) { 
    printf(" WARNING: Computation reached maximum iterations\n");  
    printf("          Maybe the solution is not good. Try again\n");            
  }
 
  printf("\n\n Given that the initial distribution was a step distribution\n with probability Pk and a merger fraction FFk in each x interval, the best estimation is:\n");
  
  /* Aqui van los resultados. */   
  
 for (i=0;i<HD.nbinx;i++)      
 {     
      printf("element %3d  Pk  %-17.15f +/- %-17.15f (N=%6d) in interval  x = [%g-%g] \n",i,(HD).p[i],sqrt((HD).covp[i][i]),(int)(exp((HD).p[i])*(float)Data2D.ndata*(HD.xbin[i+1]-HD.xbin[i])),HD.xbin[i],HD.xbin[i+1]);
      printf("element %3d  FFk %-17.15f +/- %-17.15f\n", i,(HD).p[HD.nbinx+i],sqrt((HD).covp[HD.nbinx+i][HD.nbinx+i]),HD.xbin[i],HD.xbin[i+1]);       
 }

    
  
  /* No me parece útil
  printf("\n The same distribution as above but in histogram format, instead of probabilities is:\n");
  for(j=0;j<HD.nbiny;j++) 
  {
    for (i=0;i<HD.nbinx;i++) 
    {
       printf(" element kl %2d%2d   N %-17.15f +/- %-17.15f in interval   x = [%g-%g] y = [%g-%g] \n",i,j,(HD).p[j*(HD).nbinx+i]*(float)Data2D.ndata*((HD).xbin[i+1]-(HD).xbin[i])*((HD).ybin[j+1]-(HD).ybin[j]),sqrt((HD).covp[j*(HD).nbinx+i][j*(HD).nbinx+i])*(float)Data2D.ndata*((HD).xbin[i+1]-(HD).xbin[i])*((HD).ybin[j+1]-(HD).ybin[j]),(HD).xbin[i],(HD).xbin[i+1],HD.ybin[j],HD.ybin[j+1]);
    }
  }
  */

  free(HD.p);  
  free(HD.xbin);
  free_matrix_d(HD.covp,HD.nbinx*2,HD.nbinx*2); 

} 

void ReadSample()
{
  static char datafile[100]="";
  static int colx=1,colerrx=2;
  float  *x,*errx;
  int *logx, *logerrx;
  int ndat;
  int i,j=0;

  printf(" Input file name with data: ");
  reads(datafile,datafile);
  printf(" Input column with x data: ");
  colx=readi(colx);
  printf(" Input column with erros in x data: ");
  colerrx=readi(colerrx);
  
  ndat=FileNLin(datafile);
  if(Data2D.ndata !=0)
  {
    free(Data2D.xdata);
    free(Data2D.errxdata);
  }
  Data.xdata     =malloc(ndat*sizeof(double));
  Data.errdata   =malloc(ndat*sizeof(double));
  x              =malloc(ndat*sizeof(float ));
  errx           =malloc(ndat*sizeof(float ));
  logx           =malloc(ndat*sizeof(int   ));
  logerrx        =malloc(ndat*sizeof(int   ));
  ReadNumcol(datafile,colx   ,x   ,logx,&ndat);
  ReadNumcol(datafile,colerrx,errx,logerrx,&ndat);
  
  j=0;
  for(i=0;i<ndat;i++) {
    if(logx[i] && logerrx[i]) {
      (Data).xdata[j]   =   (double)x[i];
      (Data).errdata[j]=(double)errx[i];
      j++;
    }
  }
  Data.ndata=j;

  free(x);
  free(errx);
  free(logx);
  free(logerrx);
}

void ReadSample2D()
{

  static char datafile[100]="";
  static int colx=1,colerrx=2,coly=3,colerry=4;
  float  *x,*errx,*y,*erry;
  int *logx,*logerrx,*logy,*logerry;
  int nrows;
  int i,j=0;

  printf(" Input file name with data: ");
  reads(datafile,datafile);
  printf(" Input column with x data: ");
  colx=readi(colx);
  printf(" Input column with errors in x data: ");
  colerrx=readi(colerrx);
  printf(" Input column with y data: ");
  coly=readi(coly);
  printf(" Input column with errors in y data: ");
  colerry=readi(colerry);
  
  nrows=FileNLin(datafile);
  if(Data2D.ndata !=0)
  {
    free(Data2D.xdata);
    free(Data2D.errxdata);
    free(Data2D.ydata);
    free(Data2D.errydata);
  }
  Data2D.xdata      =malloc(nrows*sizeof(double));
  Data2D.errxdata   =malloc(nrows*sizeof(double));
  Data2D.ydata      =malloc(nrows*sizeof(double));
  Data2D.errydata   =malloc(nrows*sizeof(double));
  x              =malloc(nrows*sizeof(float ));
  errx           =malloc(nrows*sizeof(float ));
  y              =malloc(nrows*sizeof(float ));
  erry           =malloc(nrows*sizeof(float ));
  logx           =malloc(nrows*sizeof(int   ));
  logerrx        =malloc(nrows*sizeof(int   ));
  logy           =malloc(nrows*sizeof(int   ));
  logerry        =malloc(nrows*sizeof(int   ));
  ReadNumcol(datafile,colx   ,x   ,logx,&nrows);
  ReadNumcol(datafile,colerrx,errx,logerrx,&nrows);
  ReadNumcol(datafile,coly   ,y   ,logy,&nrows);
  ReadNumcol(datafile,colerry,erry,logerry,&nrows);
  
  /* Lleno la variable Data2D */
  j=0;
  for(i=0;i<nrows;i++) 
  {
    if (logx[i] && logerrx[i] && logy[i] && logerry[i])
    {
      Data2D.xdata[j]   =   (double)x[i];
      Data2D.errxdata[j]=(double)errx[i];
      Data2D.ydata[j]   =   (double)y[i];
      Data2D.errydata[j]=(double)erry[i];
      j++;
    }
  }
  Data2D.ndata=j;
  free(x);
  free(errx);
  free(y);
  free(erry);
  free(logx);
  free(logerrx);
  free(logy);
  free(logerry);
}

void PlotHistogram(int nbin,float *xmin,float *xmax) 
{
  
  int i;

  int *histo;
  float *fhisto;

  double xmmin,xmmax;
  float  xhmin,xhmax;
  float  yhmin,yhmax;

  float *x;

  static char xlabel[100]="x";
  static char ylabel[100]="N";
  static char tittle[100]="Distribution estimation";
  static int pgcolor=2;
  static int pgstyle=1;



  cpgsch(1.4);
  cpgvstd();
  cpgscf(2);
  
  printf(" Input color for observed histogram ");
  pgcolor=readi(pgcolor);
  printf(" Input line style for observed histogram ");
  pgstyle=readi(pgstyle);
  printf(" Input X label ");
  reads(xlabel,xlabel);
  printf(" Input Y label ");
  reads(ylabel,ylabel);
  printf(" Input graph title ");
  reads(tittle,tittle);


  cpglab(xlabel,ylabel,tittle);
  
  xmmin=*xmin;
  xmmax=*xmax;
  histo=StHisto2_d(Data.ndata,Data.xdata,nbin,&xmmin,&xmmax);
  if(*xmin==0 && *xmax==0) {
    xmmin*=0.95;
    xmmax*=1.05;
  } 
  xhmin=xmmin;
  xhmax=xmmax;
  *xmin=xmmin;
  *xmax=xmmax;
  

  if((fhisto   =malloc(nbin*sizeof(float)))==NULL) { printf("I cannot dimension fhisto    of %d bytes",nbin*sizeof(float));exit(1);} 
  for(i=0;i<nbin;i++)   fhisto[i]=(float)histo[i];

  MinMax(nbin,fhisto,&yhmin,&yhmax);
  
  yhmax*=1.25; 

  cpgswin(xhmin,xhmax,0.,yhmax);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  if((x   =malloc(Data.ndata*sizeof(float)))==NULL) { printf("I cannot dimension x    of %d bytes",Data.ndata*sizeof(float));exit(1);} 
  for(i=0;i<Data.ndata;i++) x[i]=Data.xdata[i];
  cpgsci(pgcolor);
  cpgsls(pgstyle);
  cpghist(Data.ndata,x,xhmin,xhmax,nbin,1);
  cpgsls(1);
  cpgsci(1);
  
}
