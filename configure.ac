dnl Process this file with autoconf to produce a configure script.
AC_INIT([src/main.c])
AM_INIT_AUTOMAKE([refine],[0.6])

AC_ARG_WITH(SDK,
	[  --with-SDK[=ARG]          use GEOLAB SDK package [ARG=GEO_SDK]],
	[with_SDK=$withval],
	[with_SDK=$GEO_SDK])

AC_ARG_WITH(CAPRI,
	[  --with-CAPRI[=ARG]        use CAPRI package [ARG=CAPRI_ROOT]],
	[with_CAPRI=$withval],
	[with_CAPRI=$CAPRI_ROOT])

AC_ARG_WITH(KERNEL,
	[  --with-KERNEL[=ARG]       use CAPRI CAD Kernel [ARG=Native]],
	[with_KERNEL=$withval],
	[with_KERNEL=Native])

AC_ARG_WITH(FELISA,
	[  --with-FELISA[=ARG]       use FELISA lib [ARG=FELISA_LIB]],
	[with_FELISA=$withval],
	[with_FELISA=$FELISA_LIB])

dnl Define RanLib

AC_PROG_RANLIB

dnl Checks for std C stuff
AC_PROG_CC
AC_HEADER_STDC


dnl Check for SDK

have_SDK='no'
if test "$with_SDK" != 'no' 
then
   AC_CHECK_FILE([$with_SDK/MeatLib/Common.h],
                 [have_SDK='yes'],[have_SDK='no'])
fi

SDK_ARCH=`uname -s`
if test "$SDK_ARCH" = 'IRIX64'
then  
  ABI='-n32'
else
  ABI=''  
fi
SDK_ARCH="$SDK_ARCH$ABI"

if test "$have_SDK" != 'no'
then
   SDK=$with_SDK
   AC_DEFINE([HAVE_SDK])
   CPPFLAGS="$CPPFLAGS -I$SDK"
   CADGeom="-L$SDK/CADGeom -lCADGeom-$SDK_ARCH-CAPrI"
   Goolache="-L$SDK/Goolache -lGoolache-$SDK_ARCH"
   MeatLib="-L$SDK/MeatLib -lMeat-$SDK_ARCH"
   AC_SUBST([CADGeom])
   AC_SUBST([Goolache])
   AC_SUBST([MeatLib])
fi 

AC_SUBST(with_SDK)

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
   KERNEL=$with_KERNEL
   if test `uname` = 'Linux'
   then  
     CAPRI_ARCH='LINUX'
   else
     CAPRI_ARCH='SGIn32'
   fi
   AC_DEFINE([HAVE_CAPRI])
   CAPRILib="-L$CAPRI/$KERNEL/$CAPRI_ARCH -lcapri"
   if test "$KERNEL" == 'Parasolid' ; then
     if test `uname` = 'Linux'
     then  
       FG_LIB=$PARASOLID/fg.lib
       AC_CHECK_FILE([$FG_LIB],[have_FG_LIB='yes'],[have_FG_LIB='no'])
       CAPRILib="$CAPRILib $FG_LIB -L$PARASOLID -lpskernel"
     else
       FG_LIB=$CAPRI/$KERNEL/$CAPRI_ARCH/fg.lib
       AC_CHECK_FILE([$FG_LIB],[have_FG_LIB='yes'],[have_FG_LIB='no'])
       CAPRILib="$CAPRILib $FG_LIB -L$PARASOLID/uglibs -lpskernel"
     fi
   fi
   if test "$KERNEL" == 'OpenCASCADE' ; then
     CAPRILib="$CAPRILib -L$CASROOT/Linux/lib -lFWOSPlugin -lPTKernel -lTKBool -lTKCDF -lTKernel -lTKFeat -lTKGeomAlgo -lTKMath -lTKOffset -lTKPrim -lTKPShape -lTKService -lTKTopAlgo -lTKBRep -lTKG2d -lTKG3d -lTKGeomBase -lTKShHealing"
   fi
   if test "$x_libraries" != "" ; then
     CAPRILib="$CAPRILib -L$x_libraries -lX11"
   fi
   AC_SUBST([CAPRILib])
fi 


dnl Check for FELISA Lib

have_FELISA='no'
if test "$with_FELISA" != 'no' 
then
    AC_CHECK_FILE([$with_FELISA/MG_util.h],
                  [have_FELISA='yes'],[have_FELISA='no'])
fi

if test "$have_FELISA" != 'no'
then
   FELISA=$with_FELISA
   AC_DEFINE([HAVE_FELISA])
   FELISALib="-L$FELISA -lFELISA-$SDK_ARCH"
   AC_SUBST([FELISALib])
fi 

dnl output
AC_OUTPUT(Makefile
	  src/Makefile
	  test/Makefile)
