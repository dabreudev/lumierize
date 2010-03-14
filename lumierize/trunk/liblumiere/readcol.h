#ifndef READCOL_H
#define READCOL_H

#ifdef __cplusplus
extern "C" {
#endif

  int ReadDoublecol(char file[],int col,double *vector,int *lvec,int *nlin);
  int ReadNumcol(char file[],int col,float *vector,int *lvec,int *nlin);
  void LeeWord(char a[],int nw,char word[]);

#ifdef __cplusplus
}
#endif


#endif
