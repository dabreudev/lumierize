#ifndef SCHECHTER_H
#define SCHECHTER_H

#include "fermisel.h"
#include "histdist.h"
#include "surveydb.h"
#include "poselfunc.h"
#include "cosmology.h"

#ifdef __cplusplus
extern "C" {
#endif
  

  struct Schlf_M {
    double alfa;
    double erralfa;
    double phistar;
    double errphistar;
    double Mstar;
    double errMstar;

    double covaralfaphistar;
    double covaralfaMstar;
    double covarphistarMstar;
  };

  struct Schlf_L {
    double alfa;
    double erralfa;
    double phistar;
    double errphistar;
    double Lstar;  /*Esta debe ir SIN logaritmos */
    double errLstar;

    double covaralfaphistar;
    double covaralfaLstar;
    double covarphistarLstar;
  };

  struct Steplf_M 
  {
    int nbin;
    double *magni;
    double *errmagni;
    double *lnlf;     /* Estan en log(Phi_i) neperiano */
                      /* Usamos el criterio de que cuando la FL vale 0 lnlf 
                         vale menos infinito (-1/0.) */
    double *errlnlf;  /* Es el error de lf (expresado como antes) */
    double **covarlnlf;
    double *lf;
    double *errlf;
    /* Cuidado!! El bin j va desde magni[j] a magni[j+1], 
       luego magni esta dimensionado
       de nbin+1 y lf de nbin. 
       De esta manera, el intervalo es magni[j+1]-magni[j] */
    /* magni[0] contiene el valor mas negativo, 
       es decir, la magnitud mas brillante.
       Esta por tanto en orden ascendente */
    int *ngalbin; /* number of galaxies in the lf bin */
  };

  struct Steplf_L 
  {
    int nbin; 
    double *lumi;    /* Estas luminosidades estan en log(Lum_i) neperiano */
    double *errlumi; /* Es el error en luminosidades (expresado en log(Lum)) */
    double *lnlf;      /* Estan en log(Phi_i) neperiano */
                       /* Usamos el criterio de que cuando la FL vale 0 lnlf 
                         vale menos infinito (-1/0.) */
    double *errlnlf;   /* Es el error de lf (expresado como antes) */

    double **covarlnlf;
    /* Esta funcion representa la funcion de lum en unidades de 
       (#/Mpc3/W) para el intervalo de luminosidades
       exp(lumi[i]) - exp(lumi[i+1]), con lo que si pintas la funcion
       de luminosidad en luminosidades te sale constante en todo ese 
       intervalo. Otra cosa es que luego  tu definas los intervalos 
       para que sean logaritmicos y luego pintes el eje X de Lumi
       en luminosidades */
  };

  double Schechter_M(double M, struct Schlf_M lf);
  double lnSchechter_M(double M, struct Schlf_M lf);
  double Schechter_L(double L, struct Schlf_L lf);
  double Schechterdev_M(struct Schlf_M lf, double Mlow,double Mup);
  double zSchdev_M(struct Schlf_M lf, double zlow,double zup,double mlow,double mup, struct cosmo_param cosmo);
  double Schechterdev_L(struct Schlf_L lf, double Llow,double Lup);
  double zSchdev_L(struct Schlf_L lf, double zlow,double zup,double ffaint,double fbright, struct cosmo_param cosmo);
  double Int_sch_M(struct Schlf_M lf, double zlow,double zup,double mlim,struct cosmo_param cosmo);
  double Int_sch_L(struct Schlf_L lf, double zlow,double zup,double fluxlim,struct cosmo_param cosmo);
  double Int_sch_M_wC(struct Schlf_M lf, double zlow,double zup, double color_mean, double color_stddev,double mlim,struct cosmo_param cosmo);
  double Int_sch_L_wC(struct Schlf_L lf, double zlow,double zup, double color_mean, double color_stddev,double fluxlim,struct cosmo_param cosmo);
  double Int_sch_f_M(struct Schlf_M lf, double zlow,double zup,struct fermifsel_M fsel, struct cosmo_param cosmo);
  double Int_sch_f_M_wC(struct Schlf_M lf, double zlow,double zup, double color_mean, double color_stddev, struct fermifsel_M fsel, struct cosmo_param cosmo); 
  double Int_sch_f_L(struct Schlf_L lf, double zlow,double zup,struct fermifsel_L fsel,struct cosmo_param cosmo);
  double Int_sch_PO(struct Schlf_L lf, char photband[51], double gamma, double delta, double Kcoc, double zlow,double zup,double mlim, double ewlim, struct Histdist ewd, struct cosmo_param cosmo);
  double Int_sch_f_PO(struct Schlf_L lf, char photband[51], double gamma, double delta, double Kcoc, struct poselfunc fsel, struct SurveyItem si, double ewlim, struct Histdist ewd, struct cosmo_param cosmo);
  double Int_sch_s_f_PO(struct Schlf_L lf, char photband[51], double gamma, double delta, double Kcoc, struct poselfunc fsel, struct SurveyDB sdb , double ewlim, struct Histdist ewd, struct cosmo_param cosmo);
  double Int_sch_e_s_f_PO(struct Schlf_L lf, char photband[51], double gamma, double delta, double Kcoc, double A_ext, double meanebv, double meancocnii, struct poselfunc fsel, struct SurveyDB sdb , double ewlim, struct Histdist ewd, struct cosmo_param cosmo);


  double Sch_rhoz_L(double z, struct Schlf_L lf, double zlow,double zup,double fluxlim,struct cosmo_param cosmo);
  double Sch_rhoz_M(double z, struct Schlf_L lf, double zlow,double zup,double mlim,struct cosmo_param cosmo);

  void PlotSchLF_L_ov( struct Schlf_L lfsch);

#ifdef __cplusplus
}
#endif


#endif
