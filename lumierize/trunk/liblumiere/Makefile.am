## Process this file with automake to produce Makefile.in

##   This file is part of the ESO Common Pipeline Library
##   Copyright (C) 2001-2006 European Southern Observatory

##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.

##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.

##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

include $(top_srcdir)/config/common.am

AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~

AM_CPPFLAGS = $(GSL_CFLAGS)

#SUBDIRS = tests

include_HEADERS =  

noinst_HEADERS =  

lib_LTLIBRARIES = liblumiere.la

liblumiere_la_SOURCES = alloc.c Amoe_Ini_d.c Amoe_Init.d Amoeba_d.c Amoeba.c \
            Bnldev.c Constdev.c cosmology.c cpgdoble.c Expdev.c Ecu2Plac.c \
            FileNLin.c functions.c gaussint_d.c \
            gaussj.c gaussj_d.c getline.c \
            Gasdev.c hessian.c Histdev.c Histfunc.c \
            Lagr2_d.c mlamoebaerr.c mlamoebaerr_d.c minmax.c MinMax_d.c MCElipN_d.c \
            MLA_STY_L.c MLA_STY_p_M_wC.c MLA_SWML_L.c\
            MLA_SWML_M.c MLA_g_g_d.c MLA_STY_M.c \
            MLA_STY_s_p_f_PO.c MLA_SWML_p_L.c MLA_g_g_f_d.c \
            MLA_STY_p_f_L.c MLA_SWML_p_M.c \
            MLA_STY_gm_p_M.c MLA_STY_p_f_M.c \
            MLA_STY_gmz_p_f_M_wC.c MLA_STY_gc_p_M_wC.c\
            ML_g_g.c MLA_h_g_d.c MLA_STY_g_p_L.c \
            MLA_STY_p_L.c ML_g_g_corr.c MLA_h_g_f_d.c MLA_ff_gg_d.c \
            ML_g_g_corr_g.c MLA_hh_gg_d.c MLA_STY_gz_p_L.c MLA_STY_p_M.c \
            ML_g_g_d.c \
            mrqmin.c mrqmin_d.c pgLimits_d.c \
            Powdev.c Quartil_d.c random.c ReadWCScol.c\
            ReadKbd.c ReadDoublecol.c ReadNumcol.c ReadCharcol.c \
            schechter.c schechterdev.c step.c StMedia_d.c StModa_d.c StSuma1_d.c\
            StErrWeightMedia_d.c StHisto2_d.c StHisto_d.c StHisto2D_d.c\
            StHisto2DFF_d.c survey.c VVmax_M.c VVmax_L.c
liblumiere_la_LDFLAGS =
liblumiere_la_LIBADD =  -lm $(GSL_LIBS) -lcpgplot -lpgplot
liblumiere_la_DEPENDENCIES = 
