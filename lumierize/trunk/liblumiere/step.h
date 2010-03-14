#ifndef STEP_H
#define STEP_H

#include "schechter.h"

#ifdef __cplusplus
extern "C" {
#endif

  double Step_M(double M, struct Steplf_M lf);
  double Step_L(double L, struct Steplf_L lf);
  double Steplfdev_M(struct Steplf_M lf, double Mlow,double Mup);
  double zSteplfdev_M(struct Steplf_M lf, double zlow, double zup, double mlow, double mup, struct cosmo_param cosmo);
  double Steplfdev_L(struct Steplf_L lf, double Llow,double Lup);
  double zSteplfdev_L(struct Steplf_L lf, double zlow, double zup, double ffaint, double fbright, struct cosmo_param cosmo);
  double Int_step_M(struct Steplf_M lf, double zlow,double zup,double mlim,struct cosmo_param cosmo);
  double Int_step_L
  (struct Steplf_L lf, double zlow,double zup,double fluxlim,
   struct cosmo_param cosmo);
  void PrintStepLF_L(struct Steplf_L lf);
  void PrintStepLF_M(struct Steplf_M lf);
  void PlotStepLF_L(struct Steplf_L lf);
  void PlotStepLF_L_ov(struct Steplf_L lf);
  void PlotStepLF_M(struct Steplf_M lf);
  void PlotStepSchLF_L(struct Steplf_L lfstep, struct Schlf_L lfsch);
  void PlotStepSchLF_M(struct Steplf_M lfstep, struct Schlf_M lfsch);
  int  FitSch2StepLF_L
  (struct Steplf_L lfstep, struct Schlf_L *lfsch, double *chisq);
  int  FitSch2StepLF_M
  (struct Steplf_M lfstep, struct Schlf_M *lfsch, double *chisq);




#ifdef __cplusplus
}
#endif


#endif
