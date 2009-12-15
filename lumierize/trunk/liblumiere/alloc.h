#ifndef ALLOC_H
#define ALLOC_H

#ifdef __cplusplus
extern "C" {
#endif

	/* Funciones de allocateo */
	float *******tensor7_f(int n1, int n2, int n3, int n4, int n5, int n6, int n7);
	float ******tensor6_f(int n1, int n2, int n3, int n4, int n5, int n6);
	float *****tensor5_f(int n1, int n2, int n3, int n4, int n5);
	float ****tensor4_f(int n1, int n2, int n3, int n4);
	float ***tensor_f(int n1, int n2, int n3);
	float **matrix_f(int n1, int n2);
	float *vector_f(int n);
	float **vector_pf(int n);
	float ***vector_ppf(int n);
	float ****vector_pppf(int n);
	float *****vector_ppppf(int n);
	float ******vector_pppppf(int n);
	float *******vector_ppppppf(int n);
	double *******tensor7_d(int n1, int n2, int n3, int n4, int n5, int n6, int n7);
	double ******tensor6_d(int n1, int n2, int n3, int n4, int n5, int n6);
	double *****tensor5_d(int n1, int n2, int n3, int n4, int n5);
	double ****tensor4_d(int n1, int n2, int n3, int n4);
	double ***tensor_d(int n1, int n2, int n3);
	double **matrix_d(int n1, int n2);
	double *vector_d(int n);
	double **vector_pd(int n);
	double ***vector_ppd(int n);
	double ****vector_pppd(int n);
	double *****vector_ppppd(int n);
	double ******vector_pppppd(int n);
	double *******vector_ppppppd(int n);
	int *******tensor7_i(int n1, int n2, int n3, int n4, int n5, int n6, int n7);
	int ******tensor6_i(int n1, int n2, int n3, int n4, int n5, int n6);
	int *****tensor5_i(int n1, int n2, int n3, int n4, int n5);
	int ****tensor4_i(int n1, int n2, int n3, int n4);
	int ***tensor_i(int n1, int n2, int n3);
	int **matrix_i(int n1, int n2);
	int *vector_i(int n);
	int **vector_pi(int n);
	int ***vector_ppi(int n);
	int ****vector_pppi(int n);
	int *****vector_ppppi(int n);
	int ******vector_pppppi(int n);
	int *******vector_ppppppi(int n);
	char **vector_s(int n, int nchar);
	char **vector_pps(int n);
	char *alloc_s(int nchar);
	void free_tensor7_f(float  *******tensor7, int n1, int n2, int n3, int n4, int n5, int n6, int n7);
	void free_tensor7_d(double *******tensor7, int n1, int n2, int n3, int n4, int n5, int n6, int n7);
	void free_tensor7_i(int    *******tensor7, int n1, int n2, int n3, int n4, int n5, int n6, int n7);
	void free_tensor6_f(float  ******tensor6, int n1, int n2, int n3, int n4, int n5, int n6);
	void free_tensor6_d(double ******tensor6, int n1, int n2, int n3, int n4, int n5, int n6);
	void free_tensor6_i(int    ******tensor6, int n1, int n2, int n3, int n4, int n5, int n6);
	void free_tensor5_f(float  *****tensor5, int n1, int n2, int n3, int n4, int n5);
	void free_tensor5_d(double *****tensor5, int n1, int n2, int n3, int n4, int n5);
	void free_tensor5_i(int    *****tensor5, int n1, int n2, int n3, int n4, int n5);
	void free_tensor4_f(float  ****tensor4, int n1, int n2,int n3, int n4);
	void free_tensor4_d(double ****tensor4, int n1, int n2,int n3, int n4);
	void free_tensor4_i(int    ****tensor4, int n1, int n2,int n3, int n4);
	void free_tensor_f(float  ***tensor, int n1, int n2, int n3);
	void free_tensor_d(double ***tensor, int n1, int n2, int n3);
	void free_tensor_i(int    ***tensor, int n1, int n2, int n3);
	void free_matrix_f(float  **mat, int n1, int n2);
	void free_matrix_d(double **mat, int n1, int n2);
	void free_matrix_i(int    **mat, int n1, int n2);
	void free_vector_s(char **v, int n, int nchar);
	void DeleteRecord_f(float *vec, int n, int irecord);
	void DeleteRecord_i(int *vec, int n, int irecord);
	void DeleteRecord_d(double *vec, int n, int irecord);
	void DeleteRecord_s(char *vec, int n, int nrec, int irecord);

#ifdef __cplusplus
}
#endif


#endif
