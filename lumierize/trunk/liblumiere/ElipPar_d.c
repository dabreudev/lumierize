#include "modulos.h"

int ElipPar_d(double a, double b, double c, double d, double f, double e,
	     double *x0, double *y0, double *semax, double *semin,
	     double *t)

{

  double aa,bb,cc,dd,ee,q,xx0,yy0,tmp,z1,z2;


/*  Test: Tengo datos? */
/*  ------------------ */
  if(a == 0 && b == 0) return(0);
 
/*  Calculo el angulo. Test: t=90? */
/*  ------------------------------- */
  tmp=a-b;
  if(fabs(f) > 1E5*fabs(tmp)) {
    *t=M_PI_2;
    if(f*tmp < 0) *t = - *t;
    }
  else {
    *t=f/(a-b);
    *t=atan(*t);
    }
  *t /= 2.0;

/*  Coeficientes de la elipse en el sistema girado */
/*  ----------------------------------------------- */
  aa = a*cos(*t)*cos(*t) + b*sin(*t)*sin(*t) + f*cos(*t)*sin(*t);
  bb = a*sin(*t)*sin(*t) + b*cos(*t)*cos(*t) - f*sin(*t)*cos(*t);
  cc =  c*cos(*t) + d*sin(*t);
  dd = -c*sin(*t) + d*cos(*t);
  ee = e;

/*  Unos test */
/*  --------- */
  if( aa < 0 && bb < 0) {
    aa = -aa;
    bb = -bb;
    cc = -cc;
    dd = -dd;
    ee = -ee;
    }
  if( aa < 0 || bb < 0) return(0);
  if( fabs(cc) > 1E10*fabs(aa) || fabs(dd) > 1e10*fabs(bb)) return(0);

  if(!TestDiv0_d(cc*cc,4.0*aa,&z1,1E-19) ||
     !TestDiv0_d(dd*dd,4.0*bb,&z2,1E-19) ) return(0);
  z1 += z2-ee;
  if(!TestDiv0_d(z1,aa*bb,&q,1E-19)) return(0);
  if( q <= 0) return(0);

  if(!TestDiv0_d(-cc,2.0*aa,&xx0,1E-19)) return(0);
  if(!TestDiv0_d(-dd,2.0*bb,&yy0,1E-19)) return(0);

  *x0 = xx0*cos(*t)-yy0*sin(*t);
  *y0 = xx0*sin(*t)+yy0*cos(*t);


  *semax = sqrt(bb*q);
  *semin = sqrt(aa*q);

  if(*semax < *semin) {
    tmp = *semax;
    *semax = *semin;
    *semin = tmp;
    *t += M_PI_2;
    }

  if(*t > M_PI_2) *t -= M_PI;

  *t *= 180.0/M_PI;
  return(1);
}
