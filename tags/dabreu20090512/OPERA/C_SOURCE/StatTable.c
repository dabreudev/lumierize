#include "modulos.h"

/* Parametros del programa */

char tablefile[101];
int xcol;
int sigxcol;
int ycol;
int sigycol;

/* Variables comunes */

double *xdata;
double *ydata;
double *sigxdata;
double *sigydata;
int ndata;

void LoadParam_kbd(void);
void ReadTable(void);
void Main_menu(void);

int  main(void) {

  LoadParam_kbd();
  ReadTable();

  Main_menu();

  return(0);
}


void LoadParam_kbd(void) {


  printf(" Input file with data ");
  reads("",tablefile);
  printf(" Input column with X data ");
  xcol=readi(xcol);
  printf(" Input column with error in X data ");
  sigxcol=readi(sigxcol);
  printf(" Input column with Y data ");
  ycol=readi(ycol);
  printf(" Input column with error in Y data ");
  sigycol=readi(sigycol);


}


void ReadTable(void) {
  
  double *x,*y,*sx,*sy;
  int *ilog;

  int nlines;
  int i;

  nlines=FileNLin(tablefile);
  
  x=vector_d(nlines);
  y=vector_d(nlines);
  sx=vector_d(nlines);
  sy=vector_d(nlines);
  ilog=vector_i(nlines);

  xdata=vector_d(nlines);
  ydata=vector_d(nlines);
  sigxdata=vector_d(nlines);
  sigydata=vector_d(nlines);


  if(xcol!=0)    ReadDoublecol(tablefile,   xcol, x,ilog,&nlines);
  if(ycol!=0)    ReadDoublecol(tablefile,   ycol, y,ilog,&nlines);
  if(sigxcol!=0) ReadDoublecol(tablefile,sigxcol,sx,ilog,&nlines);
  if(sigycol!=0) ReadDoublecol(tablefile,sigycol,sy,ilog,&nlines);

  ndata=0;
  for(i=0;i<nlines;i++) 
  {
    if(ilog[i])
    {
      xdata[ndata]=x[i];
      ydata[ndata]=y[i];
      sigxdata[ndata]=sx[i];
      sigydata[ndata]=sy[i];
      ndata++;
    }
  }
}
    



void Main_menu(void) {

  char opt='E';

  double mean,sigma;
  double errmean,errsigma;
  double covarmeansigma;
  double first,median,third;

  do{
    printf(" S Compute trivial mean and sigma\n");
    printf(" W Compute mean and sigma weighted with errors\n");
    printf(" Q Compute median, first and third quartil\n");
    printf(" D Compute creation gaussian distribution\n");
    printf(" E Exit\n");
    opt=readc(opt); 
    switch (opt) { 
    case 'S':
    case 's':
      mean=StMedia_d(ndata,xdata,&sigma);
      printf(" Mean: %f   Sigma: %f\n",mean,sigma);
      break;
    case 'W':
    case 'w':
      mean=StErrWeightMedia_d(ndata,xdata,sigxdata,&sigma);
      printf(" Mean: %f   Sigma: %f\n",mean,sigma);
      break;
    case 'Q':
    case 'q':
      Quartil_d(ndata,xdata,&first,&median,&third);
      sigma=(third-first)/1.35;
      printf(" Median: %f   First quartil: %f Third quartil: %f\n",median,first,third);
      printf(" Sigma using quartil: %f\n",sigma);
      break;
    case 'D':
    case 'd':
      MLA_g_g_d(ndata,xdata,sigxdata,&mean,&sigma,&errmean,&errsigma,&covarmeansigma);
      printf(" Mean: %f +/- %f  Sigma %f +/- %f\n",mean,errmean,sigma,errsigma);
    }
  }while(opt!='E' && opt!='e');

  exit(1);

}
