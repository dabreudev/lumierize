#include "modulos.h"

static void plotbar(int i, int j, float value, float min, float max);
static void plotrect(int n, float *x, float *y, float val);

void cpghist3d(float *ima, int nx, int ny,  int i1, int i2, int j1, int j2) {


  int i,j;
  float *imaplot;
  int nxp,nyp;
  float min,max;



  nxp=i2-i1+1;
  nyp=j2-j1+1;
  imaplot=vector_f(nxp*nyp);

  cpgswin(0.,(float)(nxp+0.5*nyp),0.,15.);
  
  for(i=i1;i<=i2;i++) {
    for(j=j1;j<=j2;j++) {
      imaplot[i-i1+nxp*(j-j1)]=ima[i+nx*j];
    }
  }

  MinMax(nxp*nyp,imaplot,&min,&max);

  for(j=nyp-1;j>=0;j--) {
    for(i=0;i<nxp;i++) {
      plotbar(i,j,imaplot[i+j*nxp],min,max);
    }
  }
}

static void plotbar(int i, int j, float value, float min, float max) {

  float facx=0.5,facy=0.5;
  float value2;
  float x[4],y[4];

  value2=(value-min)/(max-min)*10.;
  if(value2<0) value2=0;
  if(value2>10) value2=10;
  x[0]=(float)(i+facx*j);
  y[0]=(float)(facy*j);
  x[1]=x[0]+1;
  y[1]=y[0];
  x[2]=x[1];
  y[2]=y[1]+value2;
  x[3]=x[0];
  y[3]=y[2];
  plotrect(4,x,y,value2);

  x[0]=x[1];
  y[0]=y[1];
  x[1]=x[0]+facx;
  y[1]=y[0]+facy;
  x[2]=x[1];
  y[2]=y[1]+value2;
  x[3]=x[0];
  y[3]=y[0]+value2;
  plotrect(4,x,y,value2);


  x[0]=x[3];
  y[0]=y[3];
  x[1]=x[2];
  y[1]=y[2];
  x[2]=x[1]-1;
  y[2]=y[1];
  x[3]=x[0]-1;
  y[3]=y[0];
  plotrect(4,x,y,value2);

}

static void plotrect(int n, float *x, float *y, float val) {
  
  cpgsfs(1);
  cpgsci((int)(val)+1);
  cpgpoly(4,x,y);
  cpgsfs(2);
  cpgsci(0);
  cpgpoly(4,x,y);
  cpgsci(1);
}
