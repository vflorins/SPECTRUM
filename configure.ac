#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# Building from zero:
# autoreconf
# automake --add-missing
# ./configure CXXFLAGS="-Ofast" --with-mpi=openmpi --with-trajectory=GUIDING --with-rkmethod=0 --with-server=CARTESIAN
# ./configure CXXFLAGS="-Ofast" --with-mpi=openmpi --with-trajectory=LORENTZ --with-rkmethod=29 --with-server=SELF

# Init commands
AC_PREREQ([2.69])
AC_INIT([spectrum_batl], [1.0])
AM_INIT_AUTOMAKE([-Wall foreign subdir-objects no-exeext])

AC_CONFIG_SRCDIR([src/trajectory_base.cc])
AC_CONFIG_HEADERS([config.h])

# Configure compilers
AC_LANG([C++])
AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
CXXFLAGS="$CXXFLAGS -std=c++17"

# Check for MPI
AC_ARG_WITH([mpi], [AS_HELP_STRING([--with-mpi=MPI], [use MPI=mpich|openmpi])], [], [])
PKG_CONFIG_PATH="/usr/lib64/pkgconfig"
MODULE_NAME=""
AS_IF([test "x$with_mpi" == "x"],
      [AC_MSG_ERROR([This package cannot be built without MPI])],
      [test $with_mpi == "mpich"],
      [path="/usr/lib64/mpich/lib/pkgconfig"
       AS_IF([test -d $path],
             [PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$path],
             [])
       MODULE_NAME="mpich"],
      [test $with_mpi == "openmpi"],
      [for path in "/usr/lib64/openmpi/lib/pkgconfig" "/opt/openmpi/lib/pkgconfig/";
       do
             AS_IF([test -d $path],
                   [PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$path]
                   [])
       done
       MODULE_NAME="ompi"
       CXXFLAGS="$CXXFLAGS -DOMPI_SKIP_MPICXX"],
      [AC_MSG_ERROR([MPI must be either mpich or openmpi])])

# If the path is empty or contains no .pc file, the following will return an error
export PKG_CONFIG_PATH
PKG_CHECK_MODULES([MPI], [$MODULE_NAME])
AS_IF([test "x$MODULE_NAME" != "x"],
      [AC_DEFINE([USE_MPI], [1], [Using MPI])],
      [])

# Check for GSL
MODULE_NAME="gsl"
PKG_CHECK_MODULES([GSL], [$MODULE_NAME])

# Check for SILO
AC_SEARCH_LIBS([DBPutQuadmesh], [silo siloh5], [AC_DEFINE([USE_SILO], [1], [Using SILO])], [AC_MSG_ERROR([SILO library was not found])])

# Check for SLURM
AC_SEARCH_LIBS([slurm_get_rem_time], [slurm], [AC_DEFINE([USE_SLURM], [1], [Using SLURM])], [AC_MSG_ERROR([SLURM library was not found])])

AC_DEFINE([TRAJ_FIELDLINE], [199], [Field line tracer (-,-,-)])
AC_DEFINE([TRAJ_LORENTZ], [200], [Newton-Lorentz model (px,py,pz)])
AC_DEFINE([TRAJ_GUIDING], [201], [Guiding center model (p_para,-,p_perp)])
AC_DEFINE([TRAJ_GUIDING_SCATT], [202], [Guiding center model with PA scattering (p_para,-,p_perp)])
AC_DEFINE([TRAJ_GUIDING_DIFF], [203], [Guiding center model with perp. diffusion (p_para,-,p_perp)])
AC_DEFINE([TRAJ_GUIDING_DIFF_SCATT], [204], [Guiding center model with perp. diffusion and PA scattering (p_para,-,p_perp)])
AC_DEFINE([TRAJ_FOCUSED], [205], [Focused transport model (p,mu,-)])
AC_DEFINE([TRAJ_PARKER], [206], [Isotropic model (p,-,-)])

# Set up the trajectory type
AC_ARG_WITH([trajectory], [AS_HELP_STRING([--with-trajectory=TRAJECTORY], [use TRAJECTORY=FIELDLINE|LORENTZ|GUIDING|GUIDING_SCATT|GUIDING_DIFF|GUIDING_DIFF_SCATT|FOCUSED|PARKER])], [], [])
AS_IF([test "x$with_trajectory" == "x"],
      [AC_MSG_ERROR([A value of TRAJECTORY is required])],
      [test $with_trajectory == "FIELDLINE" || test $with_trajectory == "LORENTZ" || test $with_trajectory == "GUIDING" || test $with_trajectory == "GUIDING_SCATT" || test $with_trajectory == "GUIDING_DIFF" || test $with_trajectory == "GUIDING_DIFF_SCATT" || test $with_trajectory == "FOCUSED" || test $with_trajectory == "PARKER"],
      [AC_DEFINE_UNQUOTED([TRAJ_TYPE], [TRAJ_$with_trajectory], [Choice of the trajectory integrator])],
      [AC_MSG_ERROR([Invalid TRAJECTORY value])])
AC_MSG_NOTICE([Using "$with_trajectory" as the trajectory model])

# Set up the RK integrator
AC_ARG_WITH([rkmethod], [AS_HELP_STRING([--with-rkmethod=RKMETHOD], [use RKMETHOD=0..29])], [], [])
AS_IF([test "x$with_rkmethod" == "x"],
      [AC_MSG_ERROR([A value of RKMETHOD is required])],
      [test $with_rkmethod -ge 0 && test $with_rkmethod -le 29],
      [AC_DEFINE_UNQUOTED([RK_INTEGRATOR_TYPE], [$with_rkmethod], [Choice of the Runge-Kutta method])],
      [AC_MSG_ERROR([RKMETHOD must be between 0 and 29])])
AC_MSG_NOTICE([Using "$with_rkmethod" as the Runge-Kutta merhod])

AC_DEFINE([SERVER_SELF], [299], [No server])
AC_DEFINE([SERVER_CARTESIAN], [300], [Cartesian server with uniform grid])
AC_DEFINE([SERVER_BATL], [301], [Cartesian server with AMR (part of BATS-R-US)])

# Set up the server type
AC_ARG_WITH([server], [AS_HELP_STRING([--with-server=SERVER], [use SERVER=SELF|CARTESIAN|BATL])], [], [])
AS_IF([test "x$with_server" == "x"],
      [AC_MSG_ERROR([A value of SERVER is required])],
      [test $with_server == "SELF" || test $with_server == "CARTESIAN" || test $with_server == "BATL"],
      [AC_DEFINE_UNQUOTED([SERVER_TYPE], [SERVER_$with_server], [Choice of the background server])],
      [AC_MSG_ERROR([Invalid SERVER value])])
AC_MSG_NOTICE([Using "$with_server" as the background server])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_CHECK_FUNCS([memset pow sqrt])

# CUDA option
AC_ARG_ENABLE([cuda], [AS_HELP_STRING([--enable-cuda], [enable CUDA [default=no]])], [], [])
AS_IF([test "x$enable_cuda" == "xyes"],
      [AC_PROG_CC[nvcc]
       AC_PROG_CXX[nvcc]
       AC_DEFINE([USE_CUDA], [1], [Using CUDA])
       CXXFLAGS="$CXXFLAGS --x=cu -dc --expt-relaxed-constexpr"
       AC_MSG_NOTICE([CUDA is enabled])],
      [AC_MSG_NOTICE([CUDA is disabled])])

# Files to generate
AC_CONFIG_FILES([Makefile benchmarks/Makefile tests/Makefile runs/Makefile])

AC_OUTPUT

