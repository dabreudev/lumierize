#include "modulos.h"  

#define ZMAX_CACHED 20.      /* For cached values in cosmo_init */
#define NUMZ_CACHED 1000000 

#define DEBUG1 1
#define DEBUG2 0

int main() 
{
  int i;
  struct cosmo_param cosmo;
  struct timeval tv1,tv2,tv3;
  struct timezone tz;
  double diff;

  double log10dlum;
  double zstep,res;
  int j;

  /* double c=299792.46;
  double q0;
  double dlum;
  double z = 0.002;

  q0 = cosmo.OM / 2.;
  dlum=c/cosmo.H0*(1-q0*(1-z)-(1-q0)*sqrt(1+2*q0*z))/(q0*q0);
  printf("dlum %g\n",dlum); */

  cosmo_init(&cosmo, 75, 1, 0);

  if(DEBUG2)
  {
    gettimeofday(&tv1,&tz);
    for(i=0;i<100000000;i++)
    {
    zstep=ZMAX_CACHED/NUMZ_CACHED;
    j=rint(2/zstep);
    log10dlum = cosmo.cached.log10dlum[j];
    res=-25+5*log10dlum+25.;
    }
    gettimeofday(&tv2,&tz);
    diff = (double)tv2.tv_sec+(double)tv2.tv_usec/1000000.;
    diff-= (double)tv1.tv_sec+(double)tv1.tv_usec/1000000.;
    printf("diff %g\n",diff);
    gettimeofday(&tv1,&tz);
    for(i=0;i<100000000;i++) mag(2,-25,cosmo);
    gettimeofday(&tv2,&tz);
    diff = (double)tv2.tv_sec+(double)tv2.tv_usec/1000000.;
    diff-= (double)tv1.tv_sec+(double)tv1.tv_usec/1000000.;
    printf("diff %g\n",diff);
    gettimeofday(&tv2,&tz);
    for(i=0;i<100000000;i++) mag_OmegaLambda0(2,-25,cosmo);
    gettimeofday(&tv3,&tz);
    diff = (double)tv3.tv_sec+(double)tv3.tv_usec/1000000.;
    diff-= (double)tv2.tv_sec+(double)tv2.tv_usec/1000000.;
    printf("diff %g\n",diff);
  }

  printf("mag 20 z 2\n");
  printf("Mag %15.10g Magold %15.10g\n",Mag(2,20,cosmo),Mag_OmegaLambda0(2,20,cosmo));
  printf("mag 22 z 3\n");
  printf("Mag %15.10g Magold %15.10g\n",Mag(3,22,cosmo),Mag_OmegaLambda0(3,22,cosmo));
  printf("mag 21 z 1\n");
  printf("Mag %15.10g Magold %15.10g\n",Mag(1,21,cosmo),Mag_OmegaLambda0(1,21,cosmo));

  printf("Mag -25 z 2\n");
  printf("mag %15.10g magold %15.10g\n",mag(2,-25,cosmo),mag_OmegaLambda0(2,-25,cosmo));
  printf("Mag -23 z 3\n");
  printf("mag %15.10g magold %15.10g\n",mag(3,-23,cosmo),mag_OmegaLambda0(3,-23,cosmo));
  printf("Mag -21 z 1\n");
  printf("mag %15.10g magold %15.10g\n",mag(1,-21,cosmo),mag_OmegaLambda0(1,-21,cosmo));
  printf("Z_m m 20 M -21\n");
  printf("Z_m %15.10g Z_mold %15.10g\n",Z_m(20,-21,cosmo),Z_m_OmegaLambda0(20,-21,cosmo));
  printf("Z_m m 21 M -21\n");
  printf("Z_m %15.10g Z_mold %15.10g\n",Z_m(21,-21,cosmo),Z_m_OmegaLambda0(21,-21,cosmo));
  printf("Z_m m 22 M -21\n");
  printf("Z_m %15.10g Z_mold %15.10g\n",Z_m(22,-21,cosmo),Z_m_OmegaLambda0(22,-21,cosmo));

  printf("dVdz z 2\n");
  printf("dVdz %15.10g dVdz %15.10g\n",dVdz(2,cosmo),dVdz_OmegaLambda0(2,cosmo));
  printf("dVdz z 3\n");
  printf("dVdz %15.10g dVdzold %15.10g\n",dVdz(3,cosmo),dVdz_OmegaLambda0(3,cosmo));
  printf("dVdz z 1\n");
  printf("dVdz %15.10g dVdzold %15.10g\n",dVdz(1,cosmo),dVdz_OmegaLambda0(1,cosmo));
  
  printf("Vol z 2\n");
  printf("Vol %15.10g Volold %15.10g\n",Vol(2,cosmo),Vol_OmegaLambda0(2,cosmo));
  printf("Vol z 3\n");
  printf("Vol %15.10g Volold %15.10g\n",Vol(3,cosmo),Vol_OmegaLambda0(3,cosmo));
  printf("Vol z 1\n");
  printf("Vol %15.10g Volold %15.10g\n",Vol(1,cosmo),Vol_OmegaLambda0(1,cosmo));

  printf("Z_m m 28.2824 M -22.3487 -------------> z > 20\n");
  printf("Z_m %15.10g Z_mold %15.10g\n",Z_m(28.2824,-22.3487,cosmo),Z_m_OmegaLambda0(28.2824,-22.3487,cosmo));

  cosmo_free(&cosmo);
  cosmo_init(&cosmo, 70, 0.9, 0.1);
  printf("New cosmology OH 70 OM 0.9 OL 0.1\n");
 
  printf("mag 20 z 2\n");
  printf("Mag %15.10g Magold %15.10g\n",Mag(2,20,cosmo),Mag_OmegaLambda0(2,20,cosmo));
  printf("mag 22 z 3\n");
  printf("Mag %15.10g Magold %15.10g\n",Mag(3,22,cosmo),Mag_OmegaLambda0(3,22,cosmo));
  printf("mag 21 z 1\n");
  printf("Mag %15.10g Magold %15.10g\n",Mag(1,21,cosmo),Mag_OmegaLambda0(1,21,cosmo));

  printf("Mag -25 z 2\n");
  printf("mag %15.10g magold %15.10g\n",mag(2,-25,cosmo),mag_OmegaLambda0(2,-25,cosmo));
  printf("Mag -23 z 3\n");
  printf("mag %15.10g magold %15.10g\n",mag(3,-23,cosmo),mag_OmegaLambda0(3,-23,cosmo));
  printf("Mag -21 z 1\n");
  printf("mag %15.10g magold %15.10g\n",mag(1,-21,cosmo),mag_OmegaLambda0(1,-21,cosmo));
  printf("Z_m m 20 M -21\n");
  printf("Z_m %15.10g Z_mold %15.10g\n",Z_m(20,-21,cosmo),Z_m_OmegaLambda0(20,-21,cosmo));
  printf("Z_m m 21 M -21\n");
  printf("Z_m %15.10g Z_mold %15.10g\n",Z_m(21,-21,cosmo),Z_m_OmegaLambda0(21,-21,cosmo));
  printf("Z_m m 22 M -21\n");
  printf("Z_m %15.10g Z_mold %15.10g\n",Z_m(22,-21,cosmo),Z_m_OmegaLambda0(22,-21,cosmo));

  printf("dVdz z 2\n");
  printf("dVdz %15.10g dVdz %15.10g\n",dVdz(2,cosmo),dVdz_OmegaLambda0(2,cosmo));
  printf("dVdz z 3\n");
  printf("dVdz %15.10g dVdzold %15.10g\n",dVdz(3,cosmo),dVdz_OmegaLambda0(3,cosmo));
  printf("dVdz z 1\n");
  printf("dVdz %15.10g dVdzold %15.10g\n",dVdz(1,cosmo),dVdz_OmegaLambda0(1,cosmo));
  
  printf("Vol z 2\n");
  printf("Vol %15.10g Volold %15.10g\n",Vol(2,cosmo),Vol_OmegaLambda0(2,cosmo));
  printf("Vol z 3\n");
  printf("Vol %15.10g Volold %15.10g\n",Vol(3,cosmo),Vol_OmegaLambda0(3,cosmo));
  printf("Vol z 1\n");
  printf("Vol %15.10g Volold %15.10g\n",Vol(1,cosmo),Vol_OmegaLambda0(1,cosmo));

  cosmo_free(&cosmo);
  cosmo_init(&cosmo, 70, 0.3, 0.7);
  printf("New cosmology OH 70 OM 0.3 OL 0.7\n");
 
  printf("mag 20 z 2\n");
  printf("Mag %15.10g Magold %15.10g\n",Mag(2,20,cosmo),Mag_OmegaLambda0(2,20,cosmo));
  printf("mag 22 z 3\n");
  printf("Mag %15.10g Magold %15.10g\n",Mag(3,22,cosmo),Mag_OmegaLambda0(3,22,cosmo));
  printf("mag 21 z 1\n");
  printf("Mag %15.10g Magold %15.10g\n",Mag(1,21,cosmo),Mag_OmegaLambda0(1,21,cosmo));
  printf("mag 21 z 0.001\n");
  printf("Mag %15.10g Magold %15.10g\n",Mag(0.001,21,cosmo),Mag_OmegaLambda0(0.001,21,cosmo));

  printf("Mag -25 z 2\n");
  printf("mag %15.10g magold %15.10g\n",mag(2,-25,cosmo),mag_OmegaLambda0(2,-25,cosmo));
  printf("Mag -23 z 3\n");
  printf("mag %15.10g magold %15.10g\n",mag(3,-23,cosmo),mag_OmegaLambda0(3,-23,cosmo));
  printf("Mag -21 z 1\n");
  printf("mag %15.10g magold %15.10g\n",mag(1,-21,cosmo),mag_OmegaLambda0(1,-21,cosmo));
  printf("Mag -21 z 0.001\n");
  printf("mag %15.10g magold %15.10g\n",mag(0.001,-21,cosmo),mag_OmegaLambda0(0.001,-21,cosmo));
  printf("Z_m m 20 M -21\n");
  printf("Z_m %15.10g Z_mold %15.10g\n",Z_m(20,-21,cosmo),Z_m_OmegaLambda0(20,-21,cosmo));
  printf("Z_m m 21 M -21\n");
  printf("Z_m %15.10g Z_mold %15.10g\n",Z_m(21,-21,cosmo),Z_m_OmegaLambda0(21,-21,cosmo));
  printf("Z_m m 22 M -21\n");
  printf("Z_m %15.10g Z_mold %15.10g\n",Z_m(22,-21,cosmo),Z_m_OmegaLambda0(22,-21,cosmo));

  printf("dVdz z 2\n");
  printf("dVdz %15.10g dVdz %15.10g\n",dVdz(2,cosmo),dVdz_OmegaLambda0(2,cosmo));
  printf("dVdz z 3\n");
  printf("dVdz %15.10g dVdzold %15.10g\n",dVdz(3,cosmo),dVdz_OmegaLambda0(3,cosmo));
  printf("dVdz z 1\n");
  printf("dVdz %15.10g dVdzold %15.10g\n",dVdz(1,cosmo),dVdz_OmegaLambda0(1,cosmo));
  
  printf("Vol z 2\n");
  printf("Vol %15.10g Volold %15.10g\n",Vol(2,cosmo),Vol_OmegaLambda0(2,cosmo));
  printf("Vol z 3\n");
  printf("Vol %15.10g Volold %15.10g\n",Vol(3,cosmo),Vol_OmegaLambda0(3,cosmo));
  printf("Vol z 1\n");
  printf("Vol %15.10g Volold %15.10g\n",Vol(1,cosmo),Vol_OmegaLambda0(1,cosmo));

  return 0;
}
