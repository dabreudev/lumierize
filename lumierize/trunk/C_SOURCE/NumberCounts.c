#include "modulos.h"

struct graph {

  float mmin,mmax;
  float ncmin,ncmax;
  char  bandname[51];
  int band;
  int logflag;
  
  int scolor;
  int gcolor;
};

struct ncstarpar1980 {
  float C1;
  float C2;
  float alpha;
  float beta;
  float kapa;
  float eta;
  float delta;
  float lambda;
  float mstar;
  float mdag;

};



struct ncstarpar1984 {
  float mlim;
  float mdim;
  float mbr;
  float cutoff;
  float ns;
  float mc;
  float a;
  float b;
  float d;
  float r0;
  float a0;
  float e;
  float rek;
  float aeband;
  float dm;
  float dma;

};


struct stars {
  float l;
  float b;
  struct ncstarpar1980 nc1980;
};

struct ncstarpar1980 ncsV={200.,400.,-0.2,0.01,-0.26,0.065,2.,1.5,15.,17.5};
struct ncstarpar1980 ncsB={235.,370.,-0.132,0.0,-0.175,0.06,1.5,2.0,17.,18.};
struct graph g;
struct stars st;

int pg1,pg2;

float StarCountsSB1980(float m,float bgal,float lgal, struct ncstarpar1980 p);
float GalCounts_Ferguson2000(float m, int band);

void Init();
void Plot();
void ChangePlot();
void SelectBand();
void ZoomCursor();
void ChangeGalCoord();


int main () {
  char opt='k';



  Init();
  

  while(opt!='E'&& opt!='e') {
    Plot();
    printf(" C Change plot parameters\n"); 
    printf(" G Change galactic coordinates l=%f b=%f\n",st.l*180/M_PI,st.b*180/M_PI); 
    printf(" B Select photometric band\n");
    printf(" Z Zoom with cursor\n");
    printf(" E Exit\n");
    opt=readc('E'); 
    switch (opt) { 
    case 'G':
    case 'g':
      ChangeGalCoord();
      break;
    case 'C':
    case 'c':
      ChangePlot();
      break;
    case 'B':
    case 'b':
      SelectBand();
      break;
    case 'Z':
    case 'z':
      ZoomCursor();
      break;
    }


  }


  return 0;
}

void Init() {


  g.mmin=0;
  g.mmax=30.;
  g.ncmin=1.e-6;
  g.ncmax=1.e7;
  g.logflag=1;
  g.scolor=2;
  g.gcolor=5;
  g.band=2;
  st.l=0./180.*M_PI;
  st.b=90./180.*M_PI;
  memcpy(&(st.nc1980),&ncsV,sizeof(struct ncstarpar1980));


  pg1=cpgopen("?");
  cpgask(0);
  cpgsch(1.7);
  cpgslw(1.5);
  pg2=cpgopen("?");
  cpgask(0);
  cpgsch(1.7);
  cpgslw(1.5);
}

void ChangePlot() {
  static char logflag='Y';
  if(g.logflag==1) logflag='Y';
  else logflag='N';


  printf(" Input minimum magnitude: ");
  g.mmin=readf(g.mmin);
  printf(" Input maximum magnitude: ");
  g.mmax=readf(g.mmax);
  printf(" Input minimum number counts: ");
  g.ncmin=readf(g.ncmin);
  printf(" Input maximum number counts: ");
  g.ncmax=readf(g.ncmax);
  printf("Do you want to apply logaritms in Y axis (Y/N): ");
  logflag=readc(logflag);
  if(logflag=='Y') g.logflag=1;
  else             g.logflag=0;
  
  if(g.logflag && g.ncmin==0) {
    printf(" No valid minimum mag=0 if logs. Repeat\n");
    printf(" Input minimum magnitude: ");
    g.ncmin=readf(g.ncmin);
  }
}

void SelectBand() {
  
  char bandnames[17][51]={"U","B","V","R","I","z","J","H","K","3.2","6.7","15","450","850","2800","8.5Ghz","X rays"};


  printf(" Select one of the following bands (those without numbers are not still available): \n");
  printf(" 1       U\n 2      B\n        V\n        R\n 5      I\n         z\n        J\n        H\n 6      K\n      3.2\n      6.7\n      15\n      450\n      850\n      2800\n      8.5Ghz\n      X rays\n");
  g.band=readi(g.band);
  
  strcpy(g.bandname,bandnames[g.band-1]);
  printf(" banda: %s\n",g.bandname);

}


void Plot() {
  int i;
  int np=100;
  float mag;
  double ntotstar,ntotgal;
  
  cpgslct(pg1);

  cpgpage();
  cpgvstd();
  if(g.logflag)     cpgswin(g.mmin,g.mmax,log10(g.ncmin),log10(g.ncmax));
  else        cpgswin(g.mmin,g.mmax,g.ncmin,g.ncmax);
  if(g.logflag) cpgbox("BCTNS",0,0,"VBCTNSL",0,0);
  else cpgbox("BCTNS",0,0,"BCTNSV",0,0);
  for(i=0;i<np;i++) {
    mag=g.mmin+(g.mmax-g.mmin)/(np-1)*i;
    cpgsci(g.scolor);
    if(g.logflag) cpgpt1(mag,log10(StarCountsSB1980(mag,st.b,st.l,st.nc1980)),3);
    else          cpgpt1(mag,StarCountsSB1980(mag,st.b,st.l,st.nc1980),3);
    cpgsci(g.gcolor);
    if(g.logflag) cpgpt1(mag,log10(GalCounts_Ferguson2000(mag,g.band)),3);
    else          cpgpt1(mag,GalCounts_Ferguson2000(mag,g.band),3);
  }
  cpgsci(1);
  cpglab(g.bandname,"Number / deg2 / mag","Differential counts");

  cpgslct(pg2);
  cpgpage();
  cpgvstd();
  if(g.logflag)     cpgswin(g.mmin,g.mmax,log10(g.ncmin),log10(g.ncmax*2.));
  else        cpgswin(g.mmin,g.mmax,g.ncmin,g.ncmax*2.);
  if(g.logflag) cpgbox("BCTNS",0,0,"VBCTNSL",0,0);
  else cpgbox("BCTNS",0,0,"BCTNSV",0,0);
  ntotstar=0.;ntotgal=0.;
  for(i=0;i<np;i++) {
    mag=g.mmin+(g.mmax-g.mmin)/(np-1)*i;
    cpgsci(g.scolor);
    ntotstar+=(g.mmax-g.mmin)/(np-1)*StarCountsSB1980(mag,st.b,st.l,st.nc1980);
    if(g.logflag) cpgpt1(mag,log10((float)ntotstar),3);
    else          cpgpt1(mag,(float)ntotstar,3);

    cpgsci(g.gcolor);
    ntotgal+=(g.mmax-g.mmin)/(np-1)*GalCounts_Ferguson2000(mag,g.band);
    printf(" M %f NC %f NC_AC %f\n",mag,GalCounts_Ferguson2000(mag,g.band),ntotgal);
    if(g.logflag) cpgpt1(mag,log10((float)ntotgal),3);
    else          cpgpt1(mag,(float)ntotgal,3);
  }
  cpgsci(1);
  cpglab(g.bandname,"Number / deg2 ","Cumulative counts");


}

void ZoomCursor() {
  char cnul;
  cpgsci(2);
  printf(" Press bottom left square with mouse...\n");
  cpgcurs(&(g.mmin),&(g.ncmin),&cnul);
  printf(" Press top right square with mouse...\n");
  cpgband(2,1,g.mmin,g.ncmin,&(g.mmax),&(g.ncmax),&cnul);
  cpgsci(1);

  if(g.logflag) {
    g.ncmin=pow(10.,g.ncmin);
    g.ncmax=pow(10.,g.ncmax);
  }
}



float StarCountsSB1980(float m,float bgal,float lgal, struct ncstarpar1980 p) {
  
  float mu,gamma,sigma;

    

  float coc1,coc2;
  float div1,div2;
  float sor1,sor2;
  float tor1,tor2;

  /*Usando un R-V tipico de R-V=-0.7 para el survey PLAS, que se supone que seran casi todo estrellas*/
/*   p.mstar=15-0.7; */
/*   p.mdag=17.5-0.7; */



  mu=0.0075*(m-12)+0.03;
  gamma=0.04*(12-m)+0.36;
  sigma=1.45-0.20*cos(bgal)*cos(lgal);

  div1=p.C1*pow(10.,p.beta*(m-p.mstar));
  div2=p.C2*pow(10.,p.eta*(m-p.mdag));

  sor1=pow(1+pow(10.,p.alpha*(m-p.mstar)),p.delta);
  tor1=pow(sin(bgal)*(1-mu*cos(lgal)/tan(bgal)),(3-5*gamma));
  sor2=pow(1+pow(10.,p.kapa*(m-p.mdag)),p.lambda);
  tor2=pow(1-cos(bgal)*cos(lgal),sigma);
  coc1=div1/(sor1*tor1);
  coc2=div2/(sor2*tor2);


/*    printf(" m %f l %f b %f div1 %f div2 %f sor1 %f sor2 %f tor1 %f tor2 %f coc1 %f coc2 %f\n",m,lgal,bgal,div1,div2,sor1,sor2,tor1,tor2,coc1,coc2);  */
  
/*   printf(" %f %f %f\n", sin(bgal),(1-mu*cos(lgal)/tan(bgal)),tan(bgal));  */
  return(coc1+coc2);
}



float StarCountsSB1984(float m,float bgal,float lgal, struct ncstarpar1984 p) {
 
 float pi=4*atan(1);
 float sinb,cosb,cosl,cosbl,cos2b;
 float psh,gsh,cfac;

 int imax;
 int nlf,nr;
 float rmax,rmin,r,dr;
 float mabs,absmag,absm;
 float z,x;
 float sh;
 float fms,fg,dnmg,dnmms,dnm,expx,vol;
 float reddening;
 float omega;
 float mg=1,mms=1;
 int indexx;
 
 float lfd[1000];
 float f[1000];

 int i,j;

 float xnmdisk[1000];

 /* Begginning computation of disk component */
 omega=pi/180.*pi/180.; 
 absm=0.165*(1.192-tan(abs(bgal)))/sin(abs(bgal));

 sinb=sin(abs(bgal));
 cosb=cos(bgal);
 cosl=cos(lgal);
 cosbl=cosb*cosl;
 cos2b=cosb*cosb;

 psh=3500;
 gsh=250;

 cfac=1.e-6;

 imax=(int)(p.mlim/p.dm);
 nlf=(int)((p.mdim-p.mbr)/p.dma) + 1;
 for(i=0;i<nlf;i++) {
   mabs=p.mbr + i*p.dma;
   f[i]= 0.44 * exp(1.5e-4 * pow(mabs+8.,3.5));
   if(f[i]>1) f[i]=1.;
 }
 dr=25.0;
 rmax=8000./sinb;
 nr=rmax/dr;
 rmin=0.0;
/*  tot=0.; */
 for(j=0;j<nr;j++) {
   r=j*dr;
   if(r<rmin) continue;
   absmag=absm*(1-exp(-sinb/p.a0*r));
   reddening=absmag/p.aeband;
   z=r*sinb;
   x=sqrt(p.r0*p.r0 + r*r*cos2b - 2*r*p.r0*cosbl);
   expx=exp(- (x-p.r0)/psh);
   vol=r*r*omega*dr*p.dma;
/*    znm=0; */
   for(i=0;i<nlf;i++) {
     mabs= p.mbr + (i-1)*p.dma;
     m=mabs + 5.0*log10(r/10.) + absmag;
     if(m> p.mlim) continue;
     sh=82.5*mabs-97.5;
     if(sh<90.) sh=90.;
     if(sh>325) sh=325;
     fms=f[i];
     fg=1-fms;
     if(mabs<-1.5) fg=0.0;
     dnmg= vol*expx*lfd[i]* fg*exp(-z/gsh);
     dnmms=vol*expx*lfd[i]*fms*exp(-z/sh);
     dnm=dnmg+dnmms;
/*      znm=znm+dnm; */
/*      BVC(1,mabs,dnmg,dnmms,dnm,m,mg,mms,reddening); */
     indexx=(int)(mg/p.dm);
     if(indexx<1) indexx=1;
     if(indexx>imax) indexx=imax;
     xnmdisk[indexx]=xnmdisk[indexx] + dnmg;
     indexx=mms/p.dm;
     if(indexx<1) indexx=1;
     if(indexx>imax) indexx=imax;
     xnmdisk[indexx]=xnmdisk[indexx] + dnmms;
   }
 }

 return (pi); /* No se que es lo que tiene que devolver en realidad */
}


/* void BVC()  { */
/* } */


float GalCounts_Ferguson2000(float m, int band) {
  /* band  i 
      U       1
      B       2
      V   
      R
      I       5
      z
      J
      H
      K       6 
     3.2
     6.7
     15
     450
     850
     2800
     8.5Ghz
     X rays             */

  int nI=53;
  float mI[53]={13.710938,13.906251,14.101563,14.257813,14.414063,14.609376,14.804688,15.039063,15.234376,15.468751,15.703126,15.898438,16.171877,16.445314,16.679689,16.953125,17.187502,17.460939,17.773438,18.046877,18.242189,18.437502,18.671875,18.906250,19.179689,19.414064,19.687502,19.960938,20.234377,20.429689,20.742189,21.015625,21.289064,21.562502,21.835938,22.187500,22.539064,22.851564,23.164062,23.554688,23.828125,24.140625,24.453125,24.765625,25.156252,25.546875,26.132812,26.601562,27.187500,27.695312,28.164062,28.593748,28.750000};
  float lognI[53]={0.415386,0.464800,0.563592,0.629463,0.711796,0.794142,0.892953,0.991780,1.107064,1.205903,1.288287,1.387126,1.502461,1.650732,1.782535,1.947293,2.079115,2.210963,2.375772,2.507642,2.590077,2.672517,2.787920,2.919799,3.051709,3.150665,3.249650,3.365115,3.464115,3.546601,3.645640,3.744661,3.860166,3.975678,4.091199,4.190308,4.289426,4.372047,4.438195,4.537368,4.652942,4.752074,4.867697,4.999812,5.099037,5.231243,5.380118,5.479452,5.562415,5.612341,5.695220,5.761581,5.761726};
 
  int nB=67;
  float mB[67]={14.453125,14.609375,14.843750,15.039062,15.234375,15.468750,15.664062,15.820312,16.015625,16.250000,16.445312,16.640625,16.875000,17.070312,17.304688,17.539062,17.734375,17.929688,18.125000,18.320312,18.515625,18.710938,18.906250,19.101562,19.257812,19.492188,19.687500,19.882812,20.078125,20.078125,20.273438,20.468750,20.664062,20.898438,21.015625,21.171875,21.328125,21.523438,21.718750,21.875000,22.031250,22.226562,22.343750,22.500000,22.617188,22.773438,22.929688,23.164062,23.398438,23.593750,23.750000,23.945312,24.218748,24.570312,24.882812,25.039062,25.273436,25.585938,25.937500,26.250000,26.601561,26.914061,27.265625,27.695312,28.007812,28.320312,28.750000};
  float lognB[67]={0.016458,0.082299,0.181065,0.263375,0.345689,0.444472,0.526796,0.625584,0.707917,0.823186,0.905529,1.004341,1.119629,1.218453,1.333755,1.432598,1.547906,1.646753,1.762073,1.844462,1.943325,2.025724,2.108127,2.190534,2.272934,2.371834,2.454255,2.553153,2.635583,2.635583,2.718018,2.800458,2.882901,2.981839,3.080734,3.163175,3.262094,3.361036,3.443507,3.542439,3.641376,3.723860,3.855740,3.938212,4.020666,4.103145,4.169148,4.251679,4.317736,4.400253,4.466270,4.565277,4.664339,4.730495,4.829597,4.895635,4.994698,5.093821,5.176495,5.259149,5.308864,5.391531,5.457746,5.507536,5.573734,5.623446,5.656761};
  
  int nU=44;
  float mU[44]={17.539061,17.734373,18.007811,18.242188,18.476561,18.476561,18.710936,18.984373,19.257811,19.257811,19.609373,19.921873,20.234373,20.468748,20.820311,21.171873,21.523436,21.914061,22.148436,22.304686,22.656248,22.968748,23.203123,23.437500,23.593750,23.828125,24.062498,24.296873,24.531250,24.804686,25.000000,25.234373,25.546873,25.781248,25.976561,26.249998,26.562498,26.796875,27.031248,27.343748,27.578123,27.812498,28.007811,28.242186};
  float lognU[44]={1.032932,1.148203,1.263482,1.395237,1.527000,1.527000,1.658771,1.790555,1.922350,1.922350,2.070636,2.251875,2.416656,2.564961,2.729778,2.878135,3.042982,3.257287,3.405654,3.521047,3.669473,3.817898,3.933345,4.032319,4.147747,4.263215,4.362207,4.461208,4.560214,4.659246,4.741762,4.791334,4.873920,4.956471,5.039007,5.088617,5.154741,5.187850,5.286913,5.369541,5.402660,5.452271,5.501862,5.534988};
  
  int nK=67;
  float mK[67]={13.750002,13.984377,14.140627,14.296877,14.453127,14.570314,14.687502,14.843752,14.960939,15.117189,15.234377,15.390627,15.507814,15.625003,15.781252,15.898439,16.015627,16.093752,16.289064,16.484377,16.640627,16.796877,17.070314,17.265627,17.460939,17.695314,17.851564,18.046877,18.203127,18.359377,18.515627,18.710939,18.906252,19.062502,19.257814,19.531252,19.726564,19.921877,20.117189,20.273439,20.468752,20.703127,20.859377,21.054689,21.289064,21.562502,21.757814,21.953127,22.148439,22.382814,22.617189,22.929689,23.242189,23.593752,23.945314,24.218750,24.531252,24.843752,25.195314,25.468752,25.703125,26.054689,26.367189,26.640627,26.992189,27.226562,27.500000};
  float lognK[67]={0.748499,0.847339,0.929695,1.012055,1.110878,1.209687,1.275580,1.357953,1.423850,1.489769,1.555670,1.621594,1.687501,1.769871,1.852266,1.918180,1.984096,2.082919,2.198275,2.280710,2.379589,2.494938,2.659763,2.758684,2.874077,2.956569,3.022544,3.105019,3.171001,3.253455,3.335912,3.401932,3.500895,3.599834,3.698807,3.764906,3.847419,3.929936,4.012457,4.094950,4.161007,4.227101,4.293131,4.359198,4.441779,4.507924,4.574003,4.640086,4.689695,4.739343,4.805473,4.888158,4.954372,5.037110,5.103374,5.153086,5.219323,5.285564,5.335368,5.368613,5.418305,5.451635,5.468441,5.501696,5.518546,5.551764,5.568537};


  int n=0;
  float mag[1000];
  float logn[1000];
  float corrmag=0;
 
  switch (band) {
  case 1:
    n=nU;memcpy(mag,mU,n*sizeof(float));memcpy(logn,lognU,n*sizeof(float));corrmag=+0.71;
    break;
  case 2:
    n=nB;memcpy(mag,mB,n*sizeof(float));memcpy(logn,lognB,n*sizeof(float));corrmag=+0.11;
    break;
  case 5:
    n=nI;memcpy(mag,mI,n*sizeof(float));memcpy(logn,lognI,n*sizeof(float));corrmag=-0.48;
    break;
  case 6:
    n=nK;memcpy(mag,mK,n*sizeof(float));memcpy(logn,lognK,n*sizeof(float));corrmag=-1.8;
    break;
    return(-1);
  }
  
  if(m<mag[0]+corrmag) return(0.);
  if(m>mag[n-1]+corrmag) return(0.);

  return(pow(10.,Lagr2(mag,logn,n,m-corrmag)));
 


}

void ChangeGalCoord() {

  printf(" Input new galactic latitude: ");
  st.b=readf(st.b*180/M_PI);
  printf(" Input new galactic longitude: ");
  st.l=readf(st.l*180/M_PI);
  st.b=st.b/180*M_PI;
  st.l=st.l/180*M_PI;
}

