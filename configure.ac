# -*- Autoconf -*-

AC_PREREQ(2.53)
AC_INIT([src/grid.c])
AM_INIT_AUTOMAKE([refine],[0.7.3])

AC_ARG_WITH(SDK,
	[  --with-SDK[=ARG]          use GEOLAB SDK package [ARG=GEO_SDK]],
	[with_SDK=$withval],         [with_SDK=$GEO_SDK])

AC_ARG_WITH(CAPRI,
	[  --with-CAPRI[=ARG]        use CAPRI package [ARG=CAPRI_ROOT]],
	[with_CAPRI=$withval],       [with_CAPRI=$CAPRI_ROOT])

AC_ARG_WITH(PARASOLID,
	[  --with-PARASOLID[=ARG]    use PARASOLID kernel [ARG=PARASOLID]],
	[with_PARASOLID=$withval],   [with_PARASOLID=$PARASOLID])

AC_ARG_WITH(OpenCASCADE,
	[  --with-OpenCASCADE[=ARG]  use OpenCASCADE kernel [ARG=CASROOT]],
	[with_OpenCASCADE=$withval], [with_OpenCASCADE=$CASROOT])

AC_PROG_RANLIB

AC_PROG_CC
AC_HEADER_STDC

arch=`uname`

dnl special arch libs

if test "$arch" = 'Linux'
then  
  ARCHlibrary=
  dynamicFlag="-rdynamic"
else
  ARCHlibrary="-lgen"
  dynamicFlag=
fi
AC_SUBST([ARCHlibrary])
AC_SUBST([dynamicFlag])

dnl Check for SDK

have_SDK='no'
if test "$with_SDK" != 'no' 
then
  AC_CHECK_FILE([$with_SDK/include/MeatLib/Common.h],
                 [have_SDK='yes'],[have_SDK='no'])
fi

if test "$have_SDK" != 'no'
then
  AC_DEFINE([HAVE_SDK])
  AM_CONDITIONAL(BUILD_SDK_DEPENDENTS,true)
  SDKheader=-I$with_SDK/include
  SDKlibrary=-L$with_SDK/lib
else
  AM_CONDITIONAL(BUILD_SDK_DEPENDENTS,false)
  SDKheader=
  SDKlibrary=
fi
AC_SUBST([SDKheader])
AC_SUBST([SDKlibrary])

dnl Check for X header

AC_PATH_X
if test "$x_includes" != ''
then
  Xheader=-I$x_includes
  Xlibrary=-L$x_libraries
else
  Xheader=
  Xlibrary=
fi
AC_SUBST([Xheader])
AC_SUBST([Xlibrary])

dnl Check for CAPRI header

have_CAPRI='no'
if test "$with_CAPRI" != 'no'
then
  CAPRIinclude=$with_CAPRI/include
  AC_CHECK_FILE([$CAPRIinclude/capri.h],
                [have_CAPRI='yes'],[have_CAPRI='no'])
fi

GridExKernels=
if test "$have_CAPRI" != 'no'
then
  AC_DEFINE([HAVE_CAPRI])
  CAPRIheader=-I$CAPRIinclude
  CAPRI=$with_CAPRI
  if test "$arch" = 'Linux'
  then  
    CAPRIarch='LINUX'
  else
    CAPRIarch='SGIn32'
  fi
  AM_CONDITIONAL(BUILD_NATIVE_KERNEL,true)
  CAPRINativelibrary="-L$CAPRI/Native/$CAPRIarch -lcapri"
else
  CAPRIheader=
  AM_CONDITIONAL(BUILD_NATIVE_KERNEL,false)
  CAPRINativelibrary=
fi
AC_SUBST([CAPRIheader])
AC_SUBST([CAPRINativelibrary])

dnl Check for Parasolid

have_PARASOLID='no'
have_FG_LIB='no'
if test "$with_PARASOLID" != 'no' 
then
  if test "$arch" == 'Linux'
  then  
    AC_CHECK_FILE([$with_PARASOLID/libpskernel.so],
                  [have_PARASOLID='yes'],[have_PARASOLID='no'])
    FG_LIB=$with_PARASOLID/fg.lib
    AC_CHECK_FILE([$FG_LIB],[have_FG_LIB='yes'],[have_FG_LIB='no'])
  else
    AC_CHECK_FILE([$with_PARASOLID/uglibs/libpskernel.so],
                  [have_PARASOLID='yes'],[have_PARASOLID='no'])
    FG_LIB=$CAPRI/Parasolid/$CAPRI_ARCH/fg.lib
    AC_CHECK_FILE([$FG_LIB],[have_FG_LIB='yes'],[have_FG_LIB='no'])
  fi
fi

if test "$have_FG_LIB" != 'no'
then
  AM_CONDITIONAL(BUILD_PARASOLID_KERNEL,true)
  if test "$arch" == 'Linux'
  then  
    CAPRIParasolidlibrary="-L$CAPRI/Parasolid/$CAPRIarch -lcapri $FG_LIB -L$with_PARASOLID -lpskernel"
  else
    CAPRIParasolidlibrary="-L$CAPRI/Parasolid/$CAPRIarch -lcapri $FG_LIB -L$with_PARASOLID/uglibs -lpskernel"
  fi
else
  AM_CONDITIONAL(BUILD_PARASOLID_KERNEL,false)
  CAPRIParasolidlibrary=
fi 
AC_SUBST([CAPRIParasolidlibrary])

have_CAS='no'
if test "$with_OpenCASCADE" != 'no'
then
  CAS_LIB=$with_OpenCASCADE/lin/lib
  AC_CHECK_FILE([$CAS_LIB/libPTKernel.so],[have_CAS='yes'],[have_CAS='no'])
fi

if test "$have_CAS" != 'no'
then
  AM_CONDITIONAL(BUILD_OPENCASCADE_KERNEL,true)
  CAPRIOpenCASCADElibrary="-L$CAPRI/OpenCASCADE/$CAPRIarch -lcapri -L$CAS_LIB -lFWOSPlugin -lPTKernel -lTKBool -lTKCDF -lTKernel -lTKFeat -lTKGeomAlgo -lTKMath -lTKOffset -lTKPrim -lTKPShape -lTKService -lTKTopAlgo -lTKBRep -lTKG2d -lTKG3d -lTKGeomBase -lTKShHealing"
else
  AM_CONDITIONAL(BUILD_OPENCASCADE_KERNEL,false)
  CAPRIOpenCASCADElibrary=
fi
AC_SUBST([CAPRIOpenCASCADElibrary])

AC_SUBST([GridExKernels])

AC_OUTPUT(Makefile src/Makefile)
