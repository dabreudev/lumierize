#include "modulos.h"
#include "gsl/gsl_math.h"
#define ITMAX 400
#define ITDERMAX 800
#define EPS 3.0e-7
#define FPMIN 1.0e-30

double gaussian(double x,double xmean,double sigma) 
{
  return(exp(-(x-xmean)*(x-xmean)/2./sigma/sigma)/sqrt(2.*M_PI)/sigma);
}

double lngaussian(double x,double xmean,double sigma) 
{
  return(-(x-xmean)*(x-xmean)/2./sigma/sigma - M_LN2 / 2. - M_LNPI / 2. - log(sigma));
}

double poidist(double x, double mean)
{
  double logp;
  logp=(x*log(mean) - mean - gammln((double)x+1.));
/*   printf(" logp %f mean %f x %f  exp %f\n",logp,mean,x,exp(logp)); */
  return(exp(logp));
}


double intgaussian(double x1, double x2, double xmean,double sigma) 
{
  /* Integral de una gaussiana desde x1 hasta x2 */
  float t1,t2;
  t1=(x1-xmean)/(sqrt(2)*sigma);
  t2=(x2-xmean)/(sqrt(2)*sigma);
  return((erf(t2)-erf(t1))/2.);
}
double int2dgaussian(double x1, double x2, double xmean, double y1, double y2, double ymean, double sigma)
{

  /* Integral de una gaussiana en dos dimensiones en el rectangulo
  desde x1 hasta x2 y desde y1 a y2*/
  float t1,t2;
  float s1,s2;
  t1=(x1-xmean)/(sqrt(2)*sigma);
  t2=(x2-xmean)/(sqrt(2)*sigma);
  s1=(y1-ymean)/(sqrt(2)*sigma);
  s2=(y2-ymean)/(sqrt(2)*sigma);
  return((erf(t2)-erf(t1))/2.*(erf(s2)-erf(s1))/2.);
}

double incom(double a,double x)
{
  int i;
  double fac=1/a;
  double s=fac;
  for(i=1;i<ITMAX;i++){
    fac*=x/(a+i);
    s+=fac;
  }
  return pow(x,a)*exp(-x)*s;
}

double lndergamm(double x)
{
  int i;
  double fac;
  double s=log((float)ITMAX);
  for(i=0;i<ITMAX;i++){
    fac=-1./(x+i);
    s+=fac;
  }
  return s;
}

double gammln(double xx)
  /* Sacada de Numerical Recipes Pag 156 (Section 6.1)
     suficientemente testada. Funciona bien.
     */
{
  static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
  double stp=2.5066282746310005;
  double x,y,tmp,ser;
  int j;

  y=x=xx;

  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for(j=0;j<=5;j++) ser += cof[j]/++y;

  return(-tmp+log(stp*ser/x));


}   



double gammq(double a,double x) 
{  
  double gamser,gammcf,gln;
  if( x < 0.0 || a <= 0.0 ) {
    printf(" Invalid arguments in routine gammq x<0 || a <=0. x=%f, a=%f\n",x,a);
    exit(1);
  }
  if(x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return (1.0-gamser);
  }
  else {
    gcf(&gammcf,a,x,&gln);
    return (gammcf);
  }
  
} 

double gammp(double a,double x) 
{   
  double gamser,gammcf,gln;
  if( x < 0.0 || a <= 0.0 ) {
    printf(" Invalid arguments in routine gammp x<0 || a <=0. x=%f, a=%f\n",x,a);
    exit(1);
  }
  if(x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return (gamser);
  }
  else {
    gcf(&gammcf,a,x,&gln);
    return (1.0-gammcf);
  }
  
} 


void gser(double *gamser, double a, double x, double *gln) {
  int n;
  double sum,del,ap;
  *gln=gammln(a);
  if (x <= 0.0) {
    if (x < 0.0)      printf(" x less than 0 in gser\n");
    *gamser=0.0;
    return;
  } 
  else {
    ap=a;
    del=sum=1.0/a;
    for(n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    printf(" Error in gser: ITMAX too small\n");
    return;
  }
}


    
void gcf( double *gammcf, double a, double x, double *gln) 
{
  int i;
  double an,b,c,d,del,h;
  
  *gln=gammln(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for(i=1;i<=ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (i > ITMAX) printf("Error in gcf: a too large, ITMAX too small\n");
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
}


double erfcc(double x)
{
  double t,z,ans;
  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
  return x >= 0.0 ? ans : 2.0-ans;
}



double  Fermi(double x,double mu,double T) {
  if(T==0) {
    if(x<mu) return(1);
    if(x>mu) return(0);
    if(x==mu) return(0.5);
  }
  return(1./(exp((x-mu)/T)+1.));
}
