#include "modulos.h"


void DeleteRecord_f(float *vec, int n, int irecord) {
  memmove(vec+irecord ,vec+irecord+1,(n-irecord-1)*sizeof(float));
}
void DeleteRecord_i(int *vec, int n, int irecord) {
  memmove(vec+irecord ,vec+irecord+1,(n-irecord-1)*sizeof(int));
}
void DeleteRecord_d(double *vec, int n, int irecord) {
  memmove(vec+irecord ,vec+irecord+1,(n-irecord-1)*sizeof(double));
}
void DeleteRecord_s(char *vec, int n, int nrec, int irecord) {
  memmove(vec+irecord ,vec+irecord+nrec,nrec*(n-irecord-1)*sizeof(char));
}
