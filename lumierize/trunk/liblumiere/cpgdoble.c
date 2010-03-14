#include <cpgplot.h>
#include <stdlib.h>
#include "cpgdoble.h""
#include "alloc.h"


void cpgline_d(int n, const double *xpts, const double *ypts) {

  float *x,*y;
  int i;

  x=vector_f(n);
  y=vector_f(n); 
  for(i=0;i<n;i++) {
    x[i]=xpts[i];
    y[i]=ypts[i];
  }
  cpgline(n,x,y);
  free(x);
  free(y);
}

void cpggray_d(const double *a, int idim, int jdim, int i1, int i2, int j1, int j2, double fg, double bg, const double *tr) {

  float *a_f;
  float fg_f, bg_f;
  float tr_f[6];
  int i;

  a_f=vector_f(idim*jdim);
  for(i=0;i<6;i++) tr_f[i]=tr[i];
  for(i=0;i<idim*jdim;i++) a_f[i]=a[i];
  fg_f=fg;
  bg_f=bg;
  cpggray(a_f, idim,  jdim,  i1, i2, j1, j2, fg_f, bg_f, tr_f);
  free(a_f);
}

void cpgcons_d(const double *a, int idim, int jdim, int i1, int i2, int j1, int j2, const double *c, int nc, const double *tr) {
  
  float *a_f;
  float *c_f;
  float tr_f[6];
  int i;

  a_f=vector_f(idim*jdim);
  c_f=vector_f(nc);
  for(i=0;i<6;i++) tr_f[i]=tr[i];
  for(i=0;i<idim*jdim;i++) a_f[i]=a[i];
  for(i=0;i<nc;i++) c_f[i]=c[i];
  cpgcons(a_f, idim, jdim, i1, i2, j1, j2, c_f, nc, tr_f);
  free(a_f);
  free(c_f);
}

void cpghist_d(int n, const double *data, double datmin, double datmax, int nbin, int pgflag) {

  float *data_f;
  int i;
  
  data_f=vector_f(n);
  for(i=0;i<n;i++) data_f[i]=data[i];
  cpghist(n, data_f, datmin,  datmax, nbin, pgflag);
  free(data_f);
}

void cpgpt_d(int n, const double *xpts, const double *ypts, int symbol) {
  float *x,*y;
  int i;

  x=vector_f(n);
  y=vector_f(n); 
  for(i=0;i<n;i++) {
    x[i]=xpts[i];
    y[i]=ypts[i];
  }
  cpgpt(n,x,y,symbol);
  free(x);
  free(y);
}


void cpgerry_d(int n, const double *x, const double *y1, const double *y2, double t) {
  float *x_f,*y1_f,*y2_f;
  int i;

  x_f=vector_f(n);
  y1_f=vector_f(n); 
  y2_f=vector_f(n); 
  for(i=0;i<n;i++) {
    x_f[i]=x[i];
    y1_f[i]=y1[i];
    y2_f[i]=y2[i];
  }
  cpgerry(n,x_f,y1_f,y2_f,t); 
  free(x_f);
  free(y1_f);
  free(y2_f);
}

int cpgcurs_d(double *x, double *y, char *ch_scalar) {

  float x_f,y_f;
  int status;
  
  status=cpgcurs(&x_f,&y_f,ch_scalar);

  *x=x_f;
  *y=y_f; 

  return(status);
}

int cpgband_d(int mode, int posn, double xref, double yref, double *x, double *y, char *ch_scalar) {

  float x_f,y_f;
  int status;

  status=cpgband(mode,posn,(float)xref,(float)yref,&x_f,&y_f,ch_scalar);

  *x=x_f;
  *y=y_f;
  return(status);
}
void cpgpoly_d(int n, const double *xpts, const double *ypts) {

  float *x,*y;
  int i;

  x=vector_f(n);
  y=vector_f(n); 
  for(i=0;i<n;i++) {
    x[i]=xpts[i];
    y[i]=ypts[i];
  }
  cpgpoly(n,x,y);
  free(x);
  free(y);
}

