#ifndef SURVEY_H
#define SURVEY_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

  struct SurveyItem {

    /*Astrometric information */
    double alfac, deltac;
    float xdim, ydim;
    float rot;
    float epoch;
    float stdx, stdy;

    /*Flux calibration variables. mag = a + b log10(flux)*/
    float a,b;
    float erra,errb,covab;
    float rms;

    /*Limiting magnitude*/
    float mmaximum;
    float mmode;
    float mfermicut;
    float deltafermi;
    float gammapowerlaw;
    float errmfermicut;
    float errdeltafermi;
    float errgammapowerlaw;

    /*Seeing variables */
    float seeing;
    float stdseeing;

    /*Image variables */
    char image[101];

    /*File cross variables */
    char calfile[101];
    int nobj;

    /* Sky britghness variables */
    float sky;
    float sigsky;

    /* Spectra variables */
    char specfile[101];
    char respfile[101];

   /* Transparency */
   float transparency;   /* con respecto a la teorica de 1, por ejemplo, o por lo menos algo proporcional a la transparencia  */

   /* Instrumental setup */
   char instsetup[51];

  };


  struct SurveyDB {

    struct SurveyItem *si;
    int nitems;

  };

  void ReadSurDB(char *dbfile, struct SurveyDB  *sdb);
  float Surveyrad(struct SurveyDB sdb);
  int  whithinimage(double ra, double dec, struct SurveyItem sitem);

#ifdef __cplusplus
}
#endif


#endif
