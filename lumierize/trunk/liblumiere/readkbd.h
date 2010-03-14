#ifndef READKBD_H
#define READKBD_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

  void reads(const char* def,char *outstr);
  char readc (char c);
  float readf(float n);
  double readd(double n);
  int readi(int n);
  void kbdpause(void);
  int getstrline(char s[], int lim);
  int fgetline(FILE *stream, char s[], int lim);


#ifdef __cplusplus
}
#endif


#endif
