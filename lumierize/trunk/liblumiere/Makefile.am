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

#SUBDIRS = tests

include_HEADERS =  amoeba.h

noinst_HEADERS =  

lib_LTLIBRARIES = liblumiere.la

liblumiere_la_SOURCES =  alloc.c Amoe_Ini_d.c Amoe_Init.d Amoeba_d.c Amoeba.c
liblumiere_la_LDFLAGS =
liblumiere_la_LIBADD =  -lm
liblumiere_la_DEPENDENCIES = 
