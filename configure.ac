# -*- Autoconf -*-

AC_PREREQ(2.53)

dnl define([git_describe], esyscmd([sh -c "git describe --tag | tr -d '\n'"]))dnl
dnl AC_INIT([refine], [git_describe], [fun3d-support@lists.nasa.gov])

AC_INIT([refine], [1.8.11], [fun3d-support@lists.nasa.gov])

AC_CONFIG_HEADERS(config.h)
AC_CONFIG_SRCDIR([src/refine_defs.h])
AM_INIT_AUTOMAKE

AC_ARG_ENABLE(exec,
	[  --enable-exec[=ARG]       build stand alone executables [ARG=yes]],
	[enable_exec=$enableval],    [enable_exec="yes"])

AC_ARG_WITH(SDK,
	[  --with-SDK[=ARG]          use GEOLAB SDK package [ARG=no]],
	[with_SDK=$withval],         [with_SDK="no"])

AC_ARG_WITH(CAPRI,
	[  --with-CAPRI[=ARG]        use CAPRI package [ARG=no]],
	[with_CAPRI=$withval],       [with_CAPRI="no"])

AC_PROG_RANLIB

AC_PROG_CC
AC_HEADER_STDC
AM_PROG_CC_C_O

AC_FC_WRAPPERS

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

dnl turn off mpi stuff in stand alone refine (set by fun3d)
AM_CONDITIONAL(BUILD_MPI,false)
AC_SUBST([MPIINC])
AC_SUBST([mpi_ldadd])

dnl Check for SDK

have_SDK='no'
if test "$with_SDK" != 'no' 
then
  AC_CHECK_FILE([$with_SDK/include/MeatLib/Common.h],
                 [have_SDK='yes'],[have_SDK='no'])
fi

if test "$have_SDK" != 'no'
then
  AC_DEFINE([HAVE_SDK],[1],[Geometry Laboratory SDK available])
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

  if test "$have_CAPRI" != 'no'
  then
    AC_DEFINE([HAVE_CAPRI],[3],[CAPrI 3.0 API])
    CAPRIheader=-I$CAPRIinclude
    CAPRIlibrary="-L$with_CAPRI/lib -lcapriDyn -ldcapri"
    AC_SUBST([CAPRIlibrary])
    AC_SUBST([CAPRIheader])
  else
    AC_MSG_ERROR([capri.h not found in $with_CAPRI/include])
  fi
fi

AC_ARG_WITH(zoltan,
	[  --with-zoltan[=ARG]       use Zoltan partitioner [ARG=no]],
	[with_zoltan=$withval],        [with_zoltan="no"])

if test "$with_zoltan" != 'no'
then
  AC_CHECK_FILE([$with_zoltan/include/zoltan.h],
                [zoltan_h_path=$with_zoltan/include],[zoltan_h_path='no'])
  AC_CHECK_FILE([$with_zoltan/lib/libzoltan.a],
                [zoltan_lib_path=$with_zoltan/lib],[zoltan_lib_path='no'])

  if test "$zoltan_h_path" != 'no'
  then
    if test "$zoltan_lib_path" != 'no'
    then
      AC_DEFINE([HAVE_ZOLTAN],[1],[Zoltan is available])
      zoltan_include="-I$zoltan_h_path"
      zoltan_ldadd="-L$zoltan_lib_path -lzoltan"
    else
      AC_MSG_ERROR([libzoltan.a not found in $with_zoltan/lib])
    fi
  else
    AC_MSG_ERROR([zoltan.h not found in $with_zoltan/include])
  fi
  AC_SUBST([zoltan_include])
  AC_SUBST([zoltan_ldadd])
fi

AC_ARG_WITH(parmetis,
	[  --with-parmetis[=ARG]       use ParMETIS partitioner [ARG=no]],
	[with_parmetis=$withval],        [with_parmetis="no"])

if test "$with_parmetis" != 'no'
then
  AC_CHECK_FILE([$with_parmetis/include/parmetis.h],
                [parmetis_h_path=$with_parmetis/include],
                [parmetis_h_path='no'])
  AC_CHECK_FILE([$with_parmetis/include/metis.h],
                [metis_h_path=$with_parmetis/include],
                [metis_h_path='no'])
  AC_CHECK_FILE([$with_parmetis/lib/libparmetis.a],
                [libparmetis_a_path=$with_parmetis/lib],
		[libparmetis_a_path='no'])
  AC_CHECK_FILE([$with_parmetis/lib/libmetis.a],
                [libmetis_a_path=$with_parmetis/lib],
		[libmetis_a_path='no'])

  if test "$parmetis_h_path" == 'no'
  then
    AC_MSG_ERROR([parmetis.h not found in $with_parmetis/include])
  fi

  if test "$metis_h_path" == 'no'
  then
    AC_MSG_ERROR([metis.h not found in $with_parmetis/include])
  fi

  if test "$libparmetis_a_path" == 'no'
  then
    AC_MSG_ERROR([libparmetis.a not found in $with_parmetis/lib])
  fi

  if test "$libmetis_a_path" == 'no'
  then
    AC_MSG_ERROR([libmetis.a not found in $with_parmetis/lib])
  fi

  parmetis_include="-I$parmetis_h_path"
  parmetis_ldadd="-L$libparmetis_a_path -lparmetis -lmetis"

  AC_SUBST([parmetis_include])
  AC_SUBST([parmetis_ldadd])
  AC_DEFINE([HAVE_PARMETIS],[1],[ParMETIS is available])
  AM_CONDITIONAL(BUILD_PARMETIS_SUPPORT,true)
else
  AM_CONDITIONAL(BUILD_PARMETIS_SUPPORT,false)
fi

AC_OUTPUT( \
	   Makefile \
	   two/Makefile \
	   src/Makefile \
	   util/Makefile \
	   )
