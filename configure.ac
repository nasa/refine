dnl Process this file with autoconf to produce a configure script.
AC_INIT([src/main.c])
AM_INIT_AUTOMAKE([refine],[0.3])

AC_ARG_WITH(SDK,
	[  --with-SDK                enable GEOLAB SDK support],
	[with_SDK=$withval],
	[with_SDK=$GEO_SDK])

AC_ARG_WITH(CAPRI,
	[  --with-CAPRI              enable CAPRI support],
	[with_CAPRI=$withval],
	[with_CAPRI=$CAPRI_ROOT])

dnl Checks for std C stuff
AC_PROG_CC
AC_HEADER_STDC

dnl Checks for external libs
AC_CHECK_FILE($with_SDK/CADGeom/CADGeom.h,AC_DEFINE([HAVE_CADGEOM]))
AC_CHECK_FILE($with_SDK/CADGeom/libCADGeom-Linux-CAPrI.a)

AC_CHECK_FILE($with_CAPRI/include/capri.h,AC_DEFINE([HAVE_CAPRI]))

AC_SUBST(with_SDK)

dnl output
AC_OUTPUT(Makefile
	  src/Makefile)
