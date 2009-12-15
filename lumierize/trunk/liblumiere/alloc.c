#include <stdlib.h>
#include <malloc.h>
#include "alloc.h"

float *******tensor7_f(int n1, int n2, int n3, int n4, int n5, int n6, int n7)  {
  
  int i;
  float *******t;
  
  t=vector_ppppppf(n1);
  for(i=0;i<n1;i++) {
    t[i]=tensor6_f(n2,n3,n4,n5,n6,n7);
  }
  return(t);
}

float ******tensor6_f(int n1, int n2, int n3, int n4, int n5, int n6)  {
  
  int i;
  float ******t;
  
  t=vector_pppppf(n1);
  for(i=0;i<n1;i++) {
    t[i]=tensor5_f(n2,n3,n4,n5,n6);
  }
  return(t);
}

float *****tensor5_f(int n1, int n2, int n3, int n4, int n5)  {
  
  int i;
  float *****t;
  
  t=vector_ppppf(n1);
  for(i=0;i<n1;i++) {
    t[i]=tensor4_f(n2,n3,n4,n5);
  }
  return(t);
}

float ****tensor4_f(int n1, int n2, int n3, int n4)  {
  
  int i;
  float ****t;
  
  t=vector_pppf(n1);
  for(i=0;i<n1;i++) {
    t[i]=tensor_f(n2,n3,n4);
  }
  return(t);
}

float ***tensor_f(int n1, int n2, int n3)  {
  
  int i;
  float ***t;
  
  t=vector_ppf(n1);
  for(i=0;i<n1;i++) {
    t[i]=matrix_f(n2,n3);
  }
  return(t);
}

float **matrix_f(int n1, int n2)  {
  
  int i;
  float **m;
  
  m=vector_pf(n1);
  for(i=0;i<n1;i++) {
    m[i]=vector_f(n2);
  }
  return(m);
}

float *vector_f(int n) {

  float *v;
  
  if((v=(float *) malloc(n*sizeof(float     )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(float));
    exit(1);
  }
  return(v);
}

float **vector_pf(int n) {

  float **v;
  
  if((v=(float **) malloc((unsigned) n*sizeof(float* )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(float *));
    exit(1);
  }
  return(v);
}

float ***vector_ppf(int n) {

  float ***v;
  
  if((v=(float ***) malloc((unsigned) n*sizeof(float** )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(float** ));
    exit(1);
  }
  return(v);
}

float ****vector_pppf(int n) {

  float ****v;
  
  if((v=(float ****) malloc((unsigned) n*sizeof(float*** )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(float*** ));
    exit(1);
  }
  return(v);
}

float *****vector_ppppf(int n) {

  float *****v;
  
  if((v=(float *****) malloc((unsigned) n*sizeof(float**** )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(float**** ));
    exit(1);
  }
  return(v);
}

float ******vector_pppppf(int n) {

  float ******v;
  
  if((v=(float ******) malloc((unsigned) n*sizeof(float***** )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(float***** ));
    exit(1);
  }
  return(v);
}

float *******vector_ppppppf(int n) {

  float *******v;
  
  if((v=(float *******) malloc((unsigned) n*sizeof(float****** )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(float****** ));
    exit(1);
  }
  return(v);
}

void free_tensor_f(float ***tensor, int n1, int n2, int n3)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free_matrix_f(tensor[i],n2,n3);
  }
  free(tensor);
  
}

void free_tensor7_f(float *******tensor7, int n1, int n2, int n3, int n4, int n5, int n6, int n7)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free_tensor6_f(tensor7[i],n2,n3,n4,n5,n6,n7);
  }
  free(tensor7);
  
}

void free_tensor6_f(float ******tensor6, int n1, int n2, int n3, int n4, int n5, int n6)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free_tensor5_f(tensor6[i],n2,n3,n4,n5,n6);
  }
  free(tensor6);
  
}

void free_tensor5_f(float *****tensor5, int n1, int n2, int n3, int n4, int n5)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free_tensor4_f(tensor5[i],n2,n3,n4,n5);
  }
  free(tensor5);
  
}

void free_tensor4_f(float ****tensor4, int n1, int n2, int n3, int n4)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free_tensor_f(tensor4[i],n2,n3,n4);
  }
  free(tensor4);
  
}

void free_matrix_f(float **matrix, int n1, int n2)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free(matrix[i]);
  }
  free(matrix);
  
}

double *******tensor7_d(int n1, int n2, int n3, int n4, int n5, int n6, int n7)  {
  
  int i;
  double *******t;
  
  t=vector_ppppppd(n1);
  for(i=0;i<n1;i++) {
    t[i]=tensor6_d(n2,n3,n4,n5,n6,n7);
  }
  return(t);
}

double ******tensor6_d(int n1, int n2, int n3, int n4, int n5, int n6)  {
  
  int i;
  double ******t;
  
  t=vector_pppppd(n1);
  for(i=0;i<n1;i++) {
    t[i]=tensor5_d(n2,n3,n4,n5,n6);
  }
  return(t);
}

double *****tensor5_d(int n1, int n2, int n3, int n4, int n5)  {
  
  int i;
  double *****t;
  
  t=vector_ppppd(n1);
  for(i=0;i<n1;i++) {
    t[i]=tensor4_d(n2,n3,n4,n5);
  }
  return(t);
}

double ****tensor4_d(int n1, int n2, int n3, int n4)  {
  
  int i;
  double ****t;
  
  t=vector_pppd(n1);
  for(i=0;i<n1;i++) {
    t[i]=tensor_d(n2,n3,n4);
  }
  return(t);
}

double ***tensor_d(int n1, int n2, int n3)  {
  
  int i;
  double ***t;
  
  t=vector_ppd(n1);
  for(i=0;i<n1;i++) {
    t[i]=matrix_d(n2,n3);
  }
  return(t);
}




double **matrix_d(int n1, int n2)  {
  
  int i;
  double **m;
  
  m=vector_pd(n1);
  for(i=0;i<n1;i++) {
    m[i]=vector_d(n2);
  }
  return(m);
}

double *vector_d(int n) {

  double *v;
  
  if((v=(double *) malloc(n*sizeof(double     )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(double));
    exit(1);
  }
  return(v);
}

double **vector_pd(int n) {

  double **v;
  
  if((v=(double **) malloc(n*sizeof(double *    )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(double *));
    exit(1);
  }
  return(v);
}

double ***vector_ppd(int n) {

  double ***v;
  
  if((v=(double ***) malloc(n*sizeof(double **   )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(double **));
    exit(1);
  }
  return(v);
}

double ****vector_pppd(int n) {

  double ****v;
  
  if((v=(double ****) malloc(n*sizeof(double ***   )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(double ***));
    exit(1);
  }
  return(v);
}

double *****vector_ppppd(int n) {

  double *****v;
  
  if((v=(double *****) malloc(n*sizeof(double ****   )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(double ****));
    exit(1);
  }
  return(v);
}

double ******vector_pppppd(int n) {

  double ******v;
  
  if((v=(double ******) malloc(n*sizeof(double *****   )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(double *****));
    exit(1);
  }
  return(v);
}

double *******vector_ppppppd(int n) {

  double *******v;
  
  if((v=(double *******) malloc(n*sizeof(double ******   )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(double ******));
    exit(1);
  }
  return(v);
}

void free_tensor7_d(double *******tensor7, int n1, int n2, int n3, int n4, int n5, int n6, int n7)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free_tensor6_d(tensor7[i],n2,n3,n4,n5,n6,n7);
  }
  free(tensor7);
  
}

void free_tensor6_d(double ******tensor6, int n1, int n2, int n3, int n4, int n5, int n6)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free_tensor5_d(tensor6[i],n2,n3,n4,n5,n6);
  }
  free(tensor6);
  
}

void free_tensor5_d(double *****tensor5, int n1, int n2, int n3, int n4, int n5)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free_tensor4_d(tensor5[i],n2,n3,n4,n5);
  }
  free(tensor5);
  
}

void free_tensor4_d(double ****tensor4, int n1, int n2, int n3, int n4)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free_tensor_d(tensor4[i],n2,n3,n4);
  }
  free(tensor4);
  
}

void free_matrix_d(double **matrix, int n1, int n2)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free(matrix[i]);
  }
  free(matrix);
}

void free_tensor_d(double ***tensor, int n1, int n2, int n3)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free_matrix_d(tensor[i],n2,n3);
  }
  free(tensor);
  
}

int *******tensor7_i(int n1, int n2, int n3, int n4, int n5, int n6, int n7)  {
  
  int i;
  int *******t;
  
  t=vector_ppppppi(n1);
  for(i=0;i<n1;i++) {
    t[i]=tensor6_i(n2,n3,n4,n5,n6,n7);
  }
  return(t);
}

int ******tensor6_i(int n1, int n2, int n3, int n4, int n5, int n6)  {
  
  int i;
  int ******t;
  
  t=vector_pppppi(n1);
  for(i=0;i<n1;i++) {
    t[i]=tensor5_i(n2,n3,n4,n5,n6);
  }
  return(t);
}

int *****tensor5_i(int n1, int n2, int n3, int n4, int n5)  {
  
  int i;
  int *****t;
  
  t=vector_ppppi(n1);
  for(i=0;i<n1;i++) {
    t[i]=tensor4_i(n2,n3,n4,n5);
  }
  return(t);
}

int ****tensor4_i(int n1, int n2, int n3, int n4)  {
  
  int i;
  int ****t;
  
  t=vector_pppi(n1);
  for(i=0;i<n1;i++) {
    t[i]=tensor_i(n2,n3,n4);
  }
  return(t);
}

int ***tensor_i(int n1, int n2, int n3)  {
  
  int i;
  int ***t;
  
  t=vector_ppi(n1);
  for(i=0;i<n1;i++) {
    t[i]=matrix_i(n2,n3);
  }
  return(t);
}

int **matrix_i(int n1, int n2)  {
  
  int i;
  int **m;
  
  m=vector_pi(n1);
  for(i=0;i<n1;i++) {
    m[i]=vector_i(n2);
  }
  return(m);
}

int *vector_i(int n) {

  int *v;
  
  if((v=(int *) malloc(n*sizeof(int     )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(int));
    exit(1);
  }
  return(v);
}

int **vector_pi(int n) {

  int **v;
  
  if((v=(int **) malloc(n*sizeof(int *    )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(int *));
    exit(1);
  }
  return(v);
}

int ***vector_ppi(int n) {

  int ***v;
  
  if((v=(int ***) malloc(n*sizeof(int **   )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(int **));
    exit(1);
  }
  return(v);
}

int ****vector_pppi(int n) {

  int ****v;
  
  if((v=(int ****) malloc(n*sizeof(int ***   )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(int ***));
    exit(1);
  }
  return(v);
}

int *****vector_ppppi(int n) {

  int *****v;
  
  if((v=(int *****) malloc(n*sizeof(int ****   )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(int ****));
    exit(1);
  }
  return(v);
}

int ******vector_pppppi(int n) {

  int ******v;
  
  if((v=(int ******) malloc(n*sizeof(int *****   )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(int *****));
    exit(1);
  }
  return(v);
}

int *******vector_ppppppi(int n) {

  int *******v;
  
  if((v=(int *******) malloc(n*sizeof(int ******   )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(int ******));
    exit(1);
  }
  return(v);
}

void free_tensor7_i(int *******tensor7, int n1, int n2, int n3, int n4, int n5, int n6, int n7)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free_tensor6_i(tensor7[i],n2,n3,n4,n5,n6,n7);
  }
  free(tensor7);
  
}

void free_tensor6_i(int ******tensor6, int n1, int n2, int n3, int n4, int n5, int n6)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free_tensor5_i(tensor6[i],n2,n3,n4,n5,n6);
  }
  free(tensor6);
  
}

void free_tensor5_i(int *****tensor5, int n1, int n2, int n3, int n4, int n5)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free_tensor4_i(tensor5[i],n2,n3,n4,n5);
  }
  free(tensor5);
  
}

void free_tensor4_i(int ****tensor4, int n1, int n2, int n3, int n4)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free_tensor_i(tensor4[i],n2,n3,n4);
  }
  free(tensor4);
  
}

void free_matrix_i(int **matrix, int n1, int n2)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free(matrix[i]);
  }
  free(matrix);
}

void free_tensor_i(int ***tensor, int n1, int n2, int n3)  {
  
  int i;
  
  for(i=0;i<n1;i++) {
    free_matrix_i(tensor[i],n2,n3); 
  }
  free(tensor);
  
}


char **vector_s(int n, int nchar) {

  int i;
  char **v;
  
  v=vector_pps(n);
  for(i=0;i<n;i++) {
    v[i]=alloc_s(nchar);
  }
  return(v);
}

char **vector_pps(int n) {
 
  char **v;
  
  if((v=(char **) malloc(n*sizeof(char *  )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",n*sizeof(char * ));
    exit(1);
  }
  return(v);
}

char *alloc_s(int nchar) {

  char *v;
  
  if((v=(char *) malloc(nchar*sizeof(char  )))==NULL) {
    printf("I cannot dimension v of %d bytes \n",nchar*sizeof(char ));
    exit(1);
  }
  return(v);
}
  

void free_vector_s(char **v, int n, int nchar) {

  int i;
  
  for(i=0;i<n;i++) {
    free(v[i]);
  }
  free(v);
}
