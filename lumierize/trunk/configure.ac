# Process this file with autoconf to produce a configure script.

AC_INIT([Luminosity function tool], [0.1], [cquike@gmail.com], [lumierize])

AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_AUX_DIR([admin])

AC_CANONICAL_BUILD

AM_INIT_AUTOMAKE

AM_MAINTAINER_MODE

AC_PREFIX_DEFAULT([/usr/local/])

# Checks for programs.
AC_PROG_CXX
AM_PROG_LIBTOOL
AC_PROG_INSTALL

#Set languaje to C++
AC_LANG_CPLUSPLUS

AC_ENABLE_SHARED(yes)

AC_SUBST(LIBTOOL_DEPS)

# Checks for header files.
AC_HEADER_STDC

AC_CONFIG_FILES([Makefile
                 liblumiere/Makefile
                 src/Makefile])

#Check for libraries
PKG_CHECK_MODULES([GSL], [gsl],
                  AC_DEFINE([HAVE_LIBGSL], [1], [1 if GSL present]),
                  AC_MSG_ERROR([No GSL available]))


AC_OUTPUT
