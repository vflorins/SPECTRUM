#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# Building from zero:
# autoreconf
# automake --add-missing
# ./configure CXXFLAGS="-O0 -g" --with-mpi=openmpi --enable-silo --with-execution=PARALLEL --with-trajectory=FIELDLINE --with-time_flow=FORWARD --with-rkmethod=0 --with-server=BATL --with-server_interp_order=1 --with-server_num_gcs=1

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
CXXFLAGS="$CXXFLAGS -std=c++20"
FCFLAGS="$FCFLAGS"
LDFLAGS="$LDFLAGS -L/usr/local/slurm/lib/slurm"

# MPI option
AC_ARG_WITH([mpi], [AS_HELP_STRING([--with-mpi=MPI], [use MPI=mpich|openmpi])], [], [])
MODULE_NAME=""
AS_IF([test "x$with_mpi" == "x"],
      [AC_MSG_ERROR([This package cannot be built without MPI])],
      [test $with_mpi == "mpich"],
      [MODULE_NAME="mpich"],
      [test $with_mpi == "openmpi"],
      [MODULE_NAME="ompi"
       CXXFLAGS="$CXXFLAGS -DOMPI_SKIP_MPICXX"],
      [AC_MSG_ERROR([MPI must be either mpich or openmpi])])

full_path=$(pkg-config --variable=prefix "$MODULE_NAME")
AS_IF([test "x$full_path" == "x"],
      [AC_MSG_ERROR([Path to MPI pkgconfig was not found])])
PKG_CHECK_MODULES([MPI], [$MODULE_NAME])
AS_IF([test "x$MODULE_NAME" != "x"],
      [AC_DEFINE([USE_MPI], [1], [Using MPI])],
      [])

# Check for GSL
MODULE_NAME="gsl"
full_path=$(pkg-config --variable=prefix "$MODULE_NAME")
AS_IF([test "x$full_path" == "x"],
      [AC_MSG_ERROR([Path to GSL pkgconfig was not found])])
PKG_CHECK_MODULES([GSL], [$MODULE_NAME])

# SILO option
AC_ARG_ENABLE([silo], [AS_HELP_STRING([--enable-silo], [enable SILO [default=no]])], [], [])
AS_IF([test "x$enable_silo" == "xyes"],
      [AC_SEARCH_LIBS([DBPutQuadmesh], [silo siloh5], [AC_DEFINE([USE_SILO], [1], [Using SILO])], [AC_MSG_WARN([SILO library was not found.])])],
      [AC_MSG_NOTICE([SILO is disabled])])

# Check for SLURM (not required)
AC_SEARCH_LIBS([slurm_get_rem_time], [slurm slurmfull], [AC_DEFINE([USE_SLURM], [1], [Using SLURM])], [AC_MSG_WARN([SLURM library was not found.])])

# Define simulation execution types
AC_DEFINE([EXEC_SERIAL], [0], [Serial execution])
AC_DEFINE([EXEC_PARALLEL], [1], [Parallel execution])

# Set up simulation execution type
AC_ARG_WITH([execution], [AS_HELP_STRING([--with-execution=EXECUTION], [use EXECUTION=SERIAL|PARALLEL])], [], [])
AS_IF([test "x$with_execution" == "x"],
      [AC_MSG_ERROR([A value of EXECUTION is required])],
      [test $with_execution == "SERIAL" || test $with_execution == "PARALLEL"],
      [AC_DEFINE_UNQUOTED([EXEC_TYPE], [EXEC_$with_execution], [Choice of simulation execution type])],
      [AC_MSG_ERROR([Invalid EXECUTION value])])
AC_MSG_NOTICE([Using "$with_execution" as the simulation execution type])

# Define trajectory types
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

# Define trajectory time flow direction types
AC_DEFINE([TRAJ_TIME_FLOW_FORWARD], [0], [Forward-in-time simulation])
AC_DEFINE([TRAJ_TIME_FLOW_BACKWARD], [1], [Backward-in-time simulation])

# Set up trajectory time flow direction
AC_ARG_WITH([time_flow], [AS_HELP_STRING([--with-time_flow=TIME_FLOW], [use TIME_FLOW=FORWARD|BACKWARD])], [], [])
AS_IF([test "x$with_time_flow" == "x"],
      [AC_MSG_ERROR([A value of TIME_FLOW is required])],
      [test $with_time_flow == "FORWARD" || test $with_time_flow == "BACKWARD"],
      [AC_DEFINE_UNQUOTED([TRAJ_TIME_FLOW], [TRAJ_TIME_FLOW_$with_time_flow], [Choice of trajectory time flow direction])],
      [AC_MSG_ERROR([Invalid TIME_FLOW value])])
AC_MSG_NOTICE([Using "$with_time_flow" as the trajectory time flow direction])

# Set up the RK integrator
AC_ARG_WITH([rkmethod], [AS_HELP_STRING([--with-rkmethod=RKMETHOD], [use RKMETHOD=0..29])], [], [])
AS_IF([test "x$with_rkmethod" == "x"],
      [AC_MSG_ERROR([A value of RKMETHOD is required])],
      [test $with_rkmethod -ge 0 && test $with_rkmethod -le 29],
      [AC_DEFINE_UNQUOTED([RK_INTEGRATOR_TYPE], [$with_rkmethod], [Choice of the Runge-Kutta method])],
      [AC_MSG_ERROR([RKMETHOD must be between 0 and 29])])
AC_MSG_NOTICE([Using "$with_rkmethod" as the Runge-Kutta method])

# Define server types
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
AS_IF([test $with_server != "SELF" && test $with_execution == "SERIAL"],
      [AC_MSG_ERROR([SERIAL execution can only used with a SELF server type])],
      [AC_MSG_NOTICE([Using "$with_server" as the background server])])

# Set up server interpolation order
AC_ARG_WITH([server_interp_order], [AS_HELP_STRING([--with-server_interp_order=SERVER_INTERP_ORDER], [use SERVER_INTERP_ORDER=-1|0|1])], [], [])
AS_IF([test "x$with_server_interp_order" == "x" && test $with_server != "SELF"],
      [AC_MSG_ERROR([A value of SERVER_INTERP_ORDER is required for all server types other than SELF])],
      [test "x$with_server_interp_order" == "x" || test $with_server_interp_order == "-1" || test $with_server_interp_order == "0" || test $with_server_interp_order == "1"],
      [AC_DEFINE_UNQUOTED([SERVER_INTERP_ORDER], [$with_server_interp_order], [Choice of background server interpolation order])],
      [AC_MSG_ERROR([Invalid SERVER_INTERP_ORDER value])])
AC_MSG_NOTICE([Using "$with_server_interp_order" as the background server interpolation order])

# Set up server number of ghost cells per side of blocks
AC_ARG_WITH([server_num_gcs], [AS_HELP_STRING([--with-server_num_gcs=SERVER_NUM_GHOST_CELLS], [use SERVER_NUM_GHOST_CELLS=0|1])], [], [])
AS_IF([test "x$with_server_num_gcs" == "x" && test $with_server != "SELF"],
      [AC_MSG_ERROR([A value of SERVER_NUM_GHOST_CELLS is required for all server types other than SELF])],
      [test "x$with_server_num_gcs" == "x" || test $with_server_num_gcs == "0" || test $with_server_num_gcs == "1"],
      [AC_DEFINE_UNQUOTED([SERVER_NUM_GHOST_CELLS], [$with_server_num_gcs], [Choice of background server number of ghost cells per side of block])],
      [AC_MSG_ERROR([Invalid SERVER_NUM_GHOST_CELLS value])])
AC_MSG_NOTICE([Using "$with_server_num_gcs" as the background server number of ghost cells per side of block])

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
AC_CONFIG_FILES([Makefile runs/Makefile])

AC_OUTPUT

