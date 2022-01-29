# -*- Autoconf -*-

AC_PREREQ(2.53)

define([git_describe], esyscmd([sh -c "git describe --tag | tr -d '\n'"]))dnl
AC_INIT([refine], [git_describe], [fun3d-support@lists.nasa.gov])

dnl AC_INIT([refine], [1.9.4], [fun3d-support@lists.nasa.gov])

AC_CONFIG_HEADERS(config.h)
AC_CONFIG_SRCDIR([src/ref_defs.h])
AM_INIT_AUTOMAKE([tar-pax])

AC_ARG_ENABLE(exec,
	[  --enable-exec[=ARG]       build stand alone executables [ARG=yes]],
	[enable_exec=$enableval],    [enable_exec="yes"])

AC_PROG_RANLIB

AC_PROG_CC
AM_PROG_CC_C_O

dnl compiler rpath
AC_ARG_WITH(compiler-rpath,
	[  --with-compiler-rpath[=ARG] use rpath for compiler [ARG=no]],
	[with_compiler_rpath=$withval],        [with_compiler_rpath="no"])
if test "$with_compiler_rpath" != 'no'
then
  compiler_rpath_ldadd="-Wl,-rpath,${with_compiler_rpath} -L${with_compiler_rpath}"
  AC_SUBST([compiler_rpath_ldadd])
fi

dnl mpi autostuff
AC_ARG_WITH(mpi,
        [[  --with-mpi[=ARG]        Path to MPI library installation [ARG=no]]],
        [with_mpi="$withval"],      [with_mpi="no"] )
AC_ARG_WITH(mpicc,
        [[  --with-mpicc[=ARG]      MPI C compiler wrapper absolute path [with_mpi/bin/mpicc]]],
        [with_mpicc="$withval"],    [with_mpicc="with_mpi/bin/mpicc"] )
AC_ARG_WITH(mpiinc,
        [[  --with-mpiinc[=ARG]     MPI include absolute path for mpif.h [with_mpi/include]]],
        [with_mpiinc="$withval"],   [with_mpiinc="with_mpi/include"] )
if test "$with_mpi" != 'no'
then
  AC_MSG_NOTICE([checking MPI installation])
  if test "$with_mpiinc" != 'no'
  then
    if test "$with_mpiinc" = 'with_mpi/include'
    then
      mpi_include="-I${with_mpi}/include"
    else
      mpi_include="-I${with_mpiinc}"
    fi
  else
    mpi_include=""
  fi
  AC_SUBST([mpi_include])

  if test "$with_mpicc" != 'no'
  then
    if test "$with_mpicc" = 'with_mpi/bin/mpicc'
    then
      MPICC="${with_mpi}/bin/mpicc"
    else
      MPICC="${with_mpicc}"
    fi
    mpi_ldadd="-L${with_mpi}/lib -lmpi"
  else
    MPICC="${CC}"
    mpi_ldadd="-L${with_mpi}/lib -lmpi"
  fi
  AC_SUBST([MPICC])
  AC_SUBST([mpi_ldadd])

  AM_CONDITIONAL(BUILD_MPI,true)
else
  AM_CONDITIONAL(BUILD_MPI,false)
fi

AC_ARG_WITH(zoltan,
	[  --with-zoltan[=ARG]       use Zoltan partitioner [ARG=no]],
	[with_zoltan=$withval],        [with_zoltan="no"])

if test "$with_zoltan" != 'no'
then
  AC_CHECK_FILE([$with_zoltan/include/zoltan.h],
                [zoltan_h_path=$with_zoltan/include],[zoltan_h_path='no'])
  AC_CHECK_FILE([$with_zoltan/lib/libzoltan.a],
                [zoltan_a_path=$with_zoltan/lib],[zoltan_a_path='no'])
  AC_CHECK_FILE([$with_zoltan/lib/libzoltan.so],
                [zoltan_so_path=$with_zoltan/lib],[zoltan_so_path='no'])
  AC_CHECK_FILE([$with_zoltan/lib/libzoltan.dylib],
                [zoltan_dylib_path=$with_zoltan/lib],[zoltan_dylib_path='no'])
  AC_CHECK_FILE([$with_zoltan/lib/libzoltan.dll],
                [zoltan_libdll_path=$with_zoltan/lib],[zoltan_libdll_path='no'])
  AC_CHECK_FILE([$with_zoltan/bin/zoltan.dll],
                [zoltan_dll_path=$with_zoltan/bin],[zoltan_dll_path='no'])

  if test "$zoltan_h_path" != 'no'
  then
    AC_DEFINE([HAVE_ZOLTAN],[1],[Zoltan is available])
    zoltan_include="-I$zoltan_h_path"
  else
    AC_MSG_ERROR([zoltan.h not found in $with_zoltan/include])
  fi
  zoltan_ldadd=""
  if test "$zoltan_a_path" != 'no'
  then
      zoltan_ldadd="-L${zoltan_a_path} -lzoltan"
  fi
  if test "$zoltan_so_path" != 'no'
  then
      zoltan_ldadd="-Wl,-rpath,${zoltan_so_path} -L${zoltan_so_path} -lzoltan"
  fi
  if test "$zoltan_dylib_path" != 'no'
  then
      zoltan_ldadd="-Wl,-rpath,${zoltan_dylib_path} -L${zoltan_dylib_path} -lzoltan"
  fi
  if test "$zoltan_libdll_path" != 'no'
  then
      zoltan_ldadd="-L${zoltan_libdll_path} -lzoltan"
  fi
  if test "$zoltan_dll_path" != 'no'
  then
      zoltan_ldadd="-L${zoltan_dll_path} -lzoltan"
  fi
  if test "empty$zoltan_ldadd" == 'empty'
  then
    AC_MSG_ERROR([libzoltan.{a,so,dylib,dll} not found in $with_zoltan/lib])
  fi

  AC_SUBST([zoltan_include])
  AC_SUBST([zoltan_ldadd])
fi

AC_ARG_WITH(parmetis,
	[  --with-parmetis[=ARG]     use ParMETIS partitioner [ARG=no]],
	[with_parmetis=$withval],        [with_parmetis="no"])

AC_ARG_WITH(metis,
	[  --with-metis[=ARG]        use METIS partitioner [ARG=ParMETIS path]],
	[with_metis=$withval],        [with_metis="$with_parmetis"])

AC_ARG_WITH(metislibdir,
	[  --with-metislibdir[=ARG]  METIS lib dir suffix [ARG=lib]],
	[with_metislibdir=$withval],        [with_metislibdir="lib"])

if test "$with_parmetis" != 'no'
then
  AC_CHECK_FILE([$with_parmetis/include/parmetis.h],
                [parmetis_h_path=$with_parmetis/include],
                [parmetis_h_path='no'])

  AC_CHECK_FILE([$with_parmetis/lib/libparmetis.a],
                [parmetis_a_path=$with_parmetis/lib],
		[parmetis_a_path='no'])
  AC_CHECK_FILE([$with_parmetis/lib/libparmetis.so],
                [parmetis_so_path=$with_parmetis/lib],
		[parmetis_so_path='no'])
  AC_CHECK_FILE([$with_parmetis/lib/libparmetis.dylib],
                [parmetis_dylib_path=$with_parmetis/lib],
		[parmetis_dylib_path='no'])
  AC_CHECK_FILE([$with_parmetis/lib/libparmetis.dll],
                [parmetis_libdll_path=$with_parmetis/lib],
		[parmetis_libdll_path='no'])
  AC_CHECK_FILE([$with_parmetis/bin/parmetis.dll],
                [parmetis_dll_path=$with_parmetis/bin],
		[parmetis_dll_path='no'])

  AC_CHECK_FILE([$with_metis/include/metis.h],
                [metis_h_path=$with_metis/include],
                [metis_h_path='no'])

  AC_CHECK_FILE([$with_metis/$with_metislibdir/libmetis.a],
                [metis_a_path=$with_metis/$with_metislibdir],
		[metis_a_path='no'])
  AC_CHECK_FILE([$with_metis/$with_metislibdir/libmetis.so],
                [metis_so_path=$with_metis/$with_metislibdir],
		[metis_so_path='no'])
  AC_CHECK_FILE([$with_metis/$with_metislibdir/libmetis.dylib],
                [metis_dylib_path=$with_metis/$with_metislibdir],
		[metis_dylib_path='no'])
  AC_CHECK_FILE([$with_metis/$with_metislibdir/libmetis.dll],
                [metis_libdll_path=$with_metis/$with_metislibdir],
		[metis_libdll_path='no'])
  AC_CHECK_FILE([$with_metis/bin/metis.dll],
                [metis_dll_path=$with_metis/bin],
		[metis_dll_path='no'])

  if test "$metis_h_path" != 'no'
  then
    metis_include="-I$metis_h_path"
  else
    AC_MSG_ERROR([metis.h not found in $with_metis/include])
  fi

  metis_ldadd=""
    if test "$metis_a_path" != 'no'
  then
      metis_ldadd="-L${metis_a_path} -lmetis"
  fi
  if test "$metis_so_path" != 'no'
  then
      metis_ldadd="-Wl,-rpath,${metis_so_path} -L${metis_so_path} -lmetis"
  fi
  if test "$metis_dylib_path" != 'no'
  then
      metis_ldadd="-Wl,-rpath,${metis_dylib_path} -L${metis_dylib_path} -lmetis"
  fi
  if test "$metis_libdll_path" != 'no'
  then
      metis_ldadd="-L${metis_libdll_path} -lmetis"
  fi
  if test "$metis_dll_path" != 'no'
  then
      metis_ldadd="-L${metis_dll_path} -lmetis"
  fi
  if test "empty$metis_ldadd" == 'empty'
  then
    AC_MSG_ERROR([libmetis.{a,so,dylib,dll} not found in $with_metis/lib])
  fi

  if test "$parmetis_h_path" != 'no'
  then
    parmetis_include="-I$parmetis_h_path ${metis_include}"
  else
    AC_MSG_ERROR([parmetis.h not found in $with_parmetis/include])
  fi

  parmetis_ldadd=""
    if test "$parmetis_a_path" != 'no'
  then
      parmetis_ldadd="-L${parmetis_a_path} -lparmetis ${metis_ldadd}"
  fi
  if test "$parmetis_so_path" != 'no'
  then
      parmetis_ldadd="-Wl,-rpath,${parmetis_so_path} -L${parmetis_so_path} -lparmetis ${metis_ldadd}"
  fi
  if test "$parmetis_dylib_path" != 'no'
  then
      parmetis_ldadd="-Wl,-rpath,${parmetis_dylib_path} -L${parmetis_dylib_path} -lparmetis ${metis_ldadd}"
  fi
  if test "$parmetis_libdll_path" != 'no'
  then
      parmetis_ldadd="-L${parmetis_libdll_path} -lparmetis ${metis_ldadd}"
  fi
  if test "$parmetis_dll_path" != 'no'
  then
      parmetis_ldadd="-L${parmetis_dll_path} -lparmetis ${metis_ldadd}"
  fi
  if test "empty$parmetis_ldadd" == 'empty'
  then
    AC_MSG_ERROR([libparmetis.{a,so,dylib,dll} not found in $with_parmetis/lib])
  fi

  AC_SUBST([parmetis_include])
  AC_SUBST([parmetis_ldadd])
  AC_DEFINE([HAVE_PARMETIS],[1],[ParMETIS is available])
  AM_CONDITIONAL(BUILD_PARMETIS_SUPPORT,true)
else
  AM_CONDITIONAL(BUILD_PARMETIS_SUPPORT,false)
fi

AC_ARG_ENABLE(lite,
        [[  --enable-lite build with EGADSlite [no]]],
        [egadslite=$enableval],[egadslite="no"])

AC_ARG_WITH(EGADS,
	[  --with-EGADS[=ARG]        use EGADS package [ARG=no]],
	[with_egads=$withval],     [with_egads="no"])

if test "$with_egads" != 'no'
then
  if test "$egadslite" != 'no'
  then
    egads_lib_root="egadslite"
    AM_CONDITIONAL(BUILD_EGADS_LITE,true)
  else
    egads_lib_root="egads"
    AM_CONDITIONAL(BUILD_EGADS_LITE,false)
  fi

  AC_CHECK_FILE([$with_egads/include/egads.h],
                [egads_h_path=$with_egads/include],[egads_h_path='no'])
  if test "$egads_h_path" != 'no'
  then
    egads_include="-I${egads_h_path}"
  else
    AC_MSG_ERROR([egads.h not found in ${with_egads}/include])
  fi

  AC_CHECK_FILE([$with_egads/lib/lib${egads_lib_root}.a],
                [egads_a_path=$with_egads/lib],[egads_a_path='no'])
  AC_CHECK_FILE([$with_egads/lib/lib${egads_lib_root}.so],
                [egads_so_path=$with_egads/lib],[egads_so_path='no'])
  AC_CHECK_FILE([$with_egads/lib/lib${egads_lib_root}.dylib],
                [egads_dylib_path=$with_egads/lib],[egads_dylib_path='no'])
  AC_CHECK_FILE([$with_egads/lib/lib${egads_lib_root}.dll],
                [egads_libdll_path=$with_egads/lib],[egads_libdll_path='no'])
  AC_CHECK_FILE([$with_egads/bin/${egads_lib_root}.dll],
                [egads_dll_path=$with_egads/bin],[egads_dll_path='no'])

  egads_ldadd=""
  if test "$egads_a_path" != 'no'
  then
      egads_ldadd="-L${egads_a_path} -l${egads_lib_root}"
  fi
  if test "$egads_so_path" != 'no'
  then
      egads_ldadd="-Wl,-rpath,${egads_so_path} -L${egads_so_path} -l${egads_lib_root}"
  fi
  if test "$egads_dylib_path" != 'no'
  then
      egads_ldadd="-Wl,-rpath,${egads_dylib_path} -L${egads_dylib_path} -l${egads_lib_root}"
  fi
  if test "$egads_libdll_path" != 'no'
  then
      egads_ldadd="-L${egads_libdll_path} -l${egads_lib_root}"
  fi
  if test "$egads_dll_path" != 'no'
  then
      egads_ldadd="-L${egads_dll_path} -l${egads_lib_root}"
  fi
  if test "empty$egads_ldadd" == 'empty'
  then
    AC_MSG_ERROR([lib${egads_lib_root}.{a,so,dylib} not found in $with_egads/lib])
  fi

save_egads_effective_CPPFLAGS="${CPPFLAGS}"
save_egads_effective_LIBS="${LIBS}"
CPPFLAGS="${egads_include} ${CPPFLAGS}"
LIBS="${egads_ldadd} ${LIBS}"
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
  [[#include "egads.h"]],
  [[ego tess = NULL; double angle = 10.0; ego effective;
    EG_initEBody(tess, angle, &effective);]])],
               [egads_effective='available'],
               [egads_effective='missing'])
CPPFLAGS="${save_egads_effective_CPPFLAGS}"
LIBS="${save_egads_effective_LIBS}"
  if test "$egads_effective" == 'available'
  then
    AC_DEFINE([HAVE_EGADS_EFFECTIVE],[1],[EGADS implements Effective Geometry])
  fi

  AC_SUBST([egads_include])
  AC_SUBST([egads_ldadd])

  AC_CHECK_FILE([$with_egads/lib/libegadslite.a],
                [egadslite_a_path=$with_egads/lib],[egadslite_a_path='no'])
  AC_CHECK_FILE([$with_egads/lib/libegadslite.so],
                [egadslite_so_path=$with_egads/lib],[egadslite_so_path='no'])
  AC_CHECK_FILE([$with_egads/lib/libegadslite.dylib],
                [egadslite_dylib_path=$with_egads/lib],[egadslite_dylib_path='no'])
  AC_CHECK_FILE([$with_egads/lib/libegadslite.dll],
                [egadslite_libdll_path=$with_egads/lib],[egadslite_libdll_path='no'])
  AC_CHECK_FILE([$with_egads/bin/egadslite.dll],
                [egadslite_dll_path=$with_egads/bin],[egadslite_dll_path='no'])

  egadslite_ldadd=""
  if test "$egadslite_a_path" != 'no'
  then
      egadslite_ldadd="-L${egadslite_a_path} -legadslite"
  fi
  if test "$egadslite_so_path" != 'no'
  then
      egadslite_ldadd="-Wl,-rpath,${egadslite_so_path} -L${egadslite_so_path} -legadslite"
  fi
  if test "$egadslite_dylib_path" != 'no'
  then
      egadslite_ldadd="-Wl,-rpath,${egadslite_dylib_path} -L${egadslite_dylib_path} -legadslite"
  fi
  if test "$egadslite_libdll_path" != 'no'
  then
      egadslite_ldadd="-L${egadslite_libdll_path} -legadslite"
  fi
  if test "$egadslite_dll_path" != 'no'
  then
      egadslite_ldadd="-L${egadslite_dll_path} -legadslite"
  fi
  if test "empty$egadslite_ldadd" == 'empty'
  then
    AC_MSG_ERROR([libegadslite.{a,so,dylib} not found in $with_egads/lib])
  fi

  AC_SUBST([egadslite_include])
  AC_SUBST([egadslite_ldadd])
  AM_CONDITIONAL(BUILD_EGADS,true)
else
  AM_CONDITIONAL(BUILD_EGADS,false)
  AM_CONDITIONAL(BUILD_EGADS_LITE,false)
fi

AC_ARG_WITH(OpenCASCADE,
	[  --with-OpenCASCADE[=ARG]  use OpenCASCADE CAD Kernel [ARG=no]],
	[with_OpenCASCADE=$withval],        [with_OpenCASCADE="no"])

if test "$with_OpenCASCADE" != 'no'
then
  AC_CHECK_FILE([$with_OpenCASCADE/lib/libTKernel.a],
                [opencascade_a_path=$with_OpenCASCADE/lib],[opencascade_a_path='no'])
  AC_CHECK_FILE([$with_OpenCASCADE/lib/libTKernel.so],
                [opencascade_so_path=$with_OpenCASCADE/lib],[opencascade_so_path='no'])
  AC_CHECK_FILE([$with_OpenCASCADE/lib/libTKernel.dylib],
                [opencascade_dylib_path=$with_OpenCASCADE/lib],[opencascade_dylib_path='no'])
  AC_CHECK_FILE([$with_OpenCASCADE/lib/libTKernel.dll],
                [opencascade_libdll_path=$with_OpenCASCADE/lib],[opencascade_libdll_path='no'])
  AC_CHECK_FILE([$with_OpenCASCADE/bin/TKernel.dll],
                [opencascade_dll_path=$with_OpenCASCADE/bin],[opencascade_dll_path='no'])

  opencascade_ldadd=""
  if test "$opencascade_a_path" != 'no'
  then
      opencascade_ldadd="-L${opencascade_a_path}"
  fi
  if test "$opencascade_so_path" != 'no'
  then
      opencascade_ldadd="-Wl,-rpath,${opencascade_so_path} -L${opencascade_so_path}"
  fi
  if test "$opencascade_dylib_path" != 'no'
  then
      opencascade_ldadd="-Wl,-rpath,${opencascade_dylib_path} -L${opencascade_dylib_path}"
  fi
  if test "$opencascade_libdll_path" != 'no'
  then
      opencascade_ldadd="-L${opencascade_libdll_path}"
  fi
  if test "$opencascade_dll_path" != 'no'
  then
      opencascade_ldadd="-L${opencascade_dll_path}"
  fi
  if test "empty$opencascade_ldadd" == 'empty'
  then
    AC_MSG_ERROR([libTKernel.{a,so,dylib,dll} for OpenCASCADE not found in $with_OpenCASCADE/lib])
  fi

  AC_SUBST([opencascade_ldadd])
fi

AC_ARG_WITH(MeshLink,
	[  --with-MeshLink[=ARG]  use MeshLink [ARG=no]],
	[with_MeshLink=$withval],        [with_MeshLink="no"])

if test "$with_MeshLink" != 'no'
then
  AC_CHECK_FILE([$with_MeshLink/src/meshlink/h/Types.h],
                [meshlink_include_path=$with_MeshLink/src/meshlink/h],
		[meshlink_inlcude_path='no'])

  meshlink_platform="none"
  if test "$meshlink_platform" == 'none'
  then
    AC_CHECK_FILE([$with_MeshLink/src/meshlink/dist/macosx/bin/meshlink.a],
                  [meshlink_platform='macosx'],
		  [meshlink_platform='none'])
  fi
  if test "$meshlink_platform" == 'none'
  then
    AC_CHECK_FILE([$with_MeshLink/src/meshlink/dist/linux_x86_64/bin/meshlink.a],
                  [meshlink_platform='linux_x86_64'],
		  [meshlink_platform='none'])
  fi
  if test "empty$meshlink_platform" == 'emptynone'
  then
    AC_MSG_ERROR([unable to find meshlink.a platform])
  fi

  AC_CHECK_FILE([$with_MeshLink/src/meshlink/dist/${meshlink_platform}/bin/meshlink.a],
                [meshlink_a_path=$with_MeshLink/src/meshlink/dist/${meshlink_platform}/bin/meshlink.a],
		[meshlink_a_path='no'])

  AC_CHECK_FILE([$with_MeshLink/src/mlparser_xerces/MeshLinkParser_xerces_c.h],
                [xerces_include_path=$with_MeshLink/src/mlparser_xerces],
		[xerces_inlcude_path='no'])

  AC_CHECK_FILE([$with_MeshLink/src/mlparser_xerces/dist/${meshlink_platform}/bin/mlparser_xerces.a],
                [xerces_a_path=$with_MeshLink/src/mlparser_xerces/dist/${meshlink_platform}/bin/mlparser_xerces.a],
		[xerces_a_path='no'])

  AC_CHECK_FILE([$with_MeshLink/external/xerces/${meshlink_platform}/lib/libxerces-c.dylib],
                [xerces_lib_path=$with_MeshLink/external/xerces/${meshlink_platform}/lib],
		[xerces_lib_path='no'])
  if test "$xerces_lib_path" == 'no'
  then
    AC_CHECK_FILE([$with_MeshLink/external/xerces/${meshlink_platform}/lib/libxerces-c.so],
                  [xerces_lib_path=$with_MeshLink/external/xerces/${meshlink_platform}/lib],
		  [xerces_lib_path='no'])
  fi

  AC_CHECK_FILE([$with_MeshLink/src/mlkernel_geode/GeomKernel_Geode_c.h],
                [geode_include_path=$with_MeshLink/src/mlkernel_geode],
		[geode_inlcude_path='no'])

  AC_CHECK_FILE([$with_MeshLink/src/mlkernel_geode/dist/${meshlink_platform}/bin/mlkernel_geode.a],
                [geode_a_path=$with_MeshLink/src/mlkernel_geode/dist/${meshlink_platform}/bin/mlkernel_geode.a],
		[geode_a_path='no'])

  AC_CHECK_FILE([$with_MeshLink/external/gelite/${meshlink_platform}/lib/libgeom.a],
                [geode_lib_path=$with_MeshLink/external/gelite/${meshlink_platform}/lib],
		[geode_lib_path='no'])

  meshlink_include=""
  if test "$meshlink_include_path" != 'no'
  then
      meshlink_include="-I${meshlink_include_path}"
  fi
  if test "$xerces_include_path" != 'no'
  then
      meshlink_include="${meshlink_include} -I${xerces_include_path}"
  fi
  if test "$geode_include_path" != 'no'
  then
      meshlink_include="${meshlink_include} -I${geode_include_path}"
  fi
  
  meshlink_ldadd=""
  if test "$meshlink_a_path" != 'no'
  then
      meshlink_ldadd="${meshlink_a_path}"
  fi

  if test "$xerces_a_path" != 'no'
  then
      meshlink_ldadd="${meshlink_ldadd} ${xerces_a_path}"
  fi
  if test "$xerces_lib_path" != 'no'
  then
      meshlink_ldadd="${meshlink_ldadd} -Wl,-rpath,${xerces_lib_path} -L${xerces_lib_path} -lxerces-c"
  fi

  if test "$geode_a_path" != 'no'
  then
      meshlink_ldadd="${meshlink_ldadd} ${geode_a_path}"
  fi
  if test "$geode_lib_path" != 'no'
  then
      meshlink_ldadd="${meshlink_ldadd} -L${geode_lib_path} -lnmb -lgeom -lvmath -lstdc++"
  fi
  
  AC_SUBST([meshlink_include])
  AC_SUBST([meshlink_ldadd])
  AM_CONDITIONAL(BUILD_MESHLINK,true)
else
  AM_CONDITIONAL(BUILD_MESHLINK,false)
fi
  
AC_CONFIG_FILES( \
	   Makefile \
	   src/Makefile \
	   util/Makefile
	   )
AC_CONFIG_FILES( \
	   acceptance/Makefile
	   )

AC_OUTPUT
