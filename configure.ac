dnl Process this file with autoconf to produce a configure script.
AC_INIT([src/main.c])
AM_INIT_AUTOMAKE([refine],[0.3])

AC_ARG_WITH(SDK,
	[  --with-SDK[=ARG]        use GEOLAB SDK package [ARG=\$GEO_SDK]],
	[with_SDK=$withval],
	[with_SDK=$GEO_SDK])

AC_ARG_WITH(CAPRI,
	[  --with-CAPRI[=ARG]      use CAPRI package [ARG=\$CAPRI_ROOT]],
	[with_CAPRI=$withval],
	[with_CAPRI=$CAPRI_ROOT])

dnl Checks for std C stuff
AC_PROG_CC
AC_HEADER_STDC


dnl Check for CAPRI header

have_CAPRI='no'
if test "$with_CAPRI" != 'no' 
then
    AC_CHECK_FILE([$with_CAPRI/include/capri.h],
                  [have_CAPRI='yes'],[have_CAPRI='no'])
fi

if test "$have_CAPRI" != 'no'
then
   AC_PATH_X
   CAPRI=$with_CAPRI
   AC_DEFINE([HAVE_CAPRI])
   CAPRILib="-L$CAPRI/Native/LINUX -lcapri -lX11"
   if test "$x_libraries" != "" ; then
       CAPRILib="$CAPRILib -L$x_libraries -lX11"
   fi
   AC_SUBST([CAPRILib])
fi 

dnl Check for SDK

have_SDK='no'
if test "$with_SDK" != 'no' 
then
   AC_CHECK_FILE([$with_SDK/MeatLib/Common.h],
                 [have_SDK='yes'],[have_SDK='no'])
fi

if test "$have_SDK" != 'no'
then
   SDK=$with_SDK
   AC_DEFINE([HAVE_SDK])
   CPPFLAGS="$CPPFLAGS -I$SDK"
   SDK_ARCH=`uname -s`
   CADGeomLib="-L$SDK/CADGeom -lCADGeom-$SDK_ARCH-CAPrI"
   AC_SUBST([CADGeomLib])
fi 



AC_SUBST(with_SDK)

dnl output
AC_OUTPUT(Makefile
	  src/Makefile)
