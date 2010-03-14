#ifndef RANDOM_H
#define RANDOM_H

#ifdef __cplusplus
extern "C" {
#endif

  /* Funciones de numeros aleatorios */

  float ran1(long *idum);
  float ran2(long *idum);
  int    Pdev(double prob);
  double Poidev(double xm);
  double Expdev(void);
  double Gasdev(void);
  double Bnldev(double pp,int n);
  double Histdev(int k,double *Pk,double *xk);
  double Powdev(double xmin, double xmax,double alfa);
  double Constdev(double xmin, double xmax);
  float ran1(long *idum);
  float ran2(long *idum);

#ifdef __cplusplus
}
#endif


#endif
