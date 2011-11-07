# -*- Autoconf -*-

AC_PREREQ(2.53)

AC_INIT([refine], [1.8.2], [fun3d-support@lists.nasa.gov])
AC_CONFIG_SRCDIR([src/refine_defs.h])
AM_INIT_AUTOMAKE

AC_ARG_ENABLE(exec,
	[  --enable-exec[=ARG]       build stand alone executables [ARG=yes]],
	[enable_exec=$enableval],    [enable_exec="yes"])

AC_ARG_ENABLE(CAPRI2,
	[  --enable-CAPRI2[=ARG]     use the CAPRI 2.0 API [ARG=no]],
	[enable_CAPRI2=$enableval],  [enable_CAPRI2="no"])

AC_ARG_WITH(SDK,
	[  --with-SDK[=ARG]          use GEOLAB SDK package [ARG=no]],
	[with_SDK=$withval],         [with_SDK="no"])

AC_ARG_WITH(CAPRI,
	[  --with-CAPRI[=ARG]        use CAPRI package [ARG=no]],
	[with_CAPRI=$withval],       [with_CAPRI="no"])

AC_ARG_WITH(PARASOLID,
	[  --with-PARASOLID[=ARG]    use PARASOLID kernel [ARG=no]],
	[with_PARASOLID=$withval],   [with_PARASOLID="no"])

AC_ARG_WITH(OpenCASCADE,
	[  --with-OpenCASCADE[=ARG]  use OpenCASCADE kernel [ARG=no]],
	[with_OpenCASCADE=$withval], [with_OpenCASCADE="no"])

AC_ARG_WITH(PROE,
	[  --with-PROE[=ARG]         use Pro-Engineer kernel [ARG=no]],
	[with_PROE=$withval],        [with_PROE="no"])

AC_PROG_RANLIB

AC_PROG_CC
AC_HEADER_STDC
AM_PROG_CC_C_O

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

dnl Enable stand alone executable build

if test "$enable_exec" != 'no'
then
  AM_CONDITIONAL(BUILD_EXECUTABLES,true)
else
  AM_CONDITIONAL(BUILD_EXECUTABLES,false)
fi

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
  if test "$enable_CAPRI2" != 'no'
  then
    AC_DEFINE([HAVE_CAPRI2],[1],[CAPrI 2.0 is available])
    AM_CONDITIONAL(BUILD_CAPRI2,true)
    CAPRIlibrary="-L$with_CAPRI/lib -ldcapri -lcapriDyn"
  else
    AM_CONDITIONAL(BUILD_CAPRI2,false)
    CAPRIlibrary=
  fi
else
  AM_CONDITIONAL(BUILD_CAPRI2,false)
fi
AC_SUBST([CAPRIlibrary])

if test "$have_CAPRI" != 'no'
then
  AC_DEFINE([HAVE_CAPRI])
  CAPRIheader=-I$CAPRIinclude
  CAPRI=$with_CAPRI
  AM_CONDITIONAL(BUILD_NATIVE_KERNEL,true)
  CAPRINativelibrary="-L$with_CAPRI/Native/$CAPRI_ARCH -lcapri"
  AC_SUBST([CAPRIheader])
  AC_SUBST([CAPRINativelibrary])
else
  AM_CONDITIONAL(BUILD_NATIVE_KERNEL,false)
fi

dnl Check for Parasolid

if test "$have_CAPRI" != 'no' && test "$enable_CAPRI2" != 'no'
then
  AM_CONDITIONAL(BUILD_NATIVE_KERNEL,false)
  AM_CONDITIONAL(BUILD_PARASOLID_KERNEL,false)
  AM_CONDITIONAL(BUILD_OPENCASCADE_KERNEL,false)
  AM_CONDITIONAL(BUILD_PROE_KERNEL,false)
else
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
      FG_LIB=$with_CAPRI/Parasolid/$CAPRI_ARCH/fg.lib
      AC_CHECK_FILE([$FG_LIB],[have_FG_LIB='yes'],[have_FG_LIB='no'])
    fi
  fi

  if test "$have_FG_LIB" != 'no'
  then
    AM_CONDITIONAL(BUILD_PARASOLID_KERNEL,true)
    if test "$arch" == 'Linux'
    then  
      CAPRIParasolidlibrary="-L$with_CAPRI/Parasolid/$CAPRI_ARCH -lcapri $FG_LIB -L$with_PARASOLID -lpskernel"
    else
      CAPRIParasolidlibrary="-L$with_CAPRI/Parasolid/$CAPRI_ARCH -lcapri $FG_LIB -L$with_PARASOLID/uglibs -lpskernel"
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
    CAPRIOpenCASCADElibrary="-L$with_CAPRI/OpenCASCADE/$CAPRI_ARCH -lcapri -L$CAS_LIB -lFWOSPlugin -lPTKernel -lTKBool -lTKCDF -lTKernel -lTKFeat -lTKGeomAlgo -lTKMath -lTKOffset -lTKPrim -lTKPShape -lTKService -lTKTopAlgo -lTKBRep -lTKG2d -lTKG3d -lTKGeomBase -lTKShHealing"
  else
    AM_CONDITIONAL(BUILD_OPENCASCADE_KERNEL,false)
    CAPRIOpenCASCADElibrary=
  fi
  AC_SUBST([CAPRIOpenCASCADElibrary])

  have_PROETool='no'
  have_PROEDev='no'
  if test "$with_PROE" != 'no'
  then
    PROE_BASE=$with_PROE
    PROE_TOOL=$PROE_BASE/protoolkit/i486_linux/obj/protk_dll.a
    AC_CHECK_FILE([$PROE_TOOL],[have_PROETool='yes'],[have_PROETool='no'])
    PROE_DEV=$PROE_BASE/prodevelop/i486_linux/obj/prodev_dll.a
    AC_CHECK_FILE([$PROE_DEV],[have_PROEDev='yes'],[have_PROEDev='no'])
  fi

  if test "$have_PROETool" != 'no' &&  test "$have_PROEDev" != 'no' 
  then
    AM_CONDITIONAL(BUILD_PROE_KERNEL,true)
    CAPRIProElibrary="-L$with_CAPRI/ProE/$CAPRI_ARCH -lcapri $PROE_TOOL $PROE_DEV"
  else
    AM_CONDITIONAL(BUILD_PROE_KERNEL,false)
    CAPRIProElibrary=
  fi
  AC_SUBST([CAPRIProElibrary])
fi

AC_OUTPUT( \
	   Makefile \
	   two/Makefile \
	   src/Makefile \
	   util/Makefile \
	   )
