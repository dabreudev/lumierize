#ifndef GAUSSINT_H
#define GAUSSINT_H

#ifdef __cplusplus
extern "C" {
#endif

  double gaussintleg_d(double (*funk)(double), double x1, double x2, int n);
  void gauleg_d(double x1, double x2, double  x[], double  w[], int n);
  double gaussinther_d(double (*funk)(double), double offset, double scale, int n);
  void gauher_d(double x[], double w[], int n);
  double gaussintlag_d(double (*funk)(double), double scale, double alfa, int n);
  void gaulag_d(double x[], double w[], int n, double alf);


#ifdef __cplusplus
}
#endif


#endif
