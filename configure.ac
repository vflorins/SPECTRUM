#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# Building from zero:
# autoreconf
# automake --add-missing

# Use for serial testing
# ./configure --enable-debug --with-execution=SERIAL --with-trajectory=FIELDLINE --with-time_flow=FORWARD --with-rkmethod=0 --with-server=SELF

# Use for parallel production run
# ./configure --with-mpi=openmpi --enable-silo --enable-gsl --with-execution=PARALLEL --with-trajectory=FIELDLINE --with-time_flow=FORWARD --with-rkmethod=0 --with-server=SELF --with-server_interp_order=1 --with-server_num_gcs=1 --with-cuda=13.0

# Init commands
AC_PREREQ([2.69])
AC_INIT([spectrum_batl], [1.0])
AM_INIT_AUTOMAKE([-Wall foreign subdir-objects no-exeext])

AC_CONFIG_SRCDIR([src/trajectory_base.cc])
AC_CONFIG_HEADERS([config.h])

# Check for programs.
AC_LANG([C++])
AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# Need to reset these because autoconf populates them with values that are not always appropriate (-O2, etc.)
CXXFLAGS="-std=c++20"
FCFLAGS=""
AC_CHECK_TOOL([LD], [ld])

# Debug option
AC_ARG_ENABLE([debug], [AS_HELP_STRING([--enable-debug], [enable debugging [default=no]])], [], [])
AS_IF([test "x$enable_debug" == "xyes"], [
   CXXFLAGS="$CXXFLAGS -g -O0"
   FFLAGS="$FCFLAGS -g -O0"
   AC_DEFINE([GEO_DEBUG], [1], [Enable debugging])
   AC_MSG_NOTICE([Debugging is enabled])
],
[
   CXXFLAGS="$CXXFLAGS -O3"
   FCFLAGS="$FCFLAGS -O3"
   AC_MSG_NOTICE([Debugging is disabled])
])

# MPI configuration
HAVE_MPI=0
MPI_CFLAGS=""
AC_ARG_WITH([mpi], [AS_HELP_STRING([--with-mpi=MPI], [use MPI=mpich|openmpi])], [], [])
MODULE_NAME=""
AS_IF([test "x$with_mpi" == "x"], [
   AC_MSG_NOTICE([Building without MPI])
],
[
   AS_IF([test "$with_mpi" == "mpich"], [
      MODULE_NAME="mpich"
   ], 
   [
      AS_IF([test "$with_mpi" == "openmpi"], [
         MODULE_NAME="ompi"
      ],
      [
         AC_MSG_ERROR([MPI must be either mpich or openmpi])
      ])
   ])

   PKG_CHECK_MODULES([MPI], [$MODULE_NAME])
   AS_IF([test "x$MODULE_NAME" != "x"], [
      AC_DEFINE([USE_MPI], [1], [Using MPI])
      MPI_CFLAGS="$MPI_CFLAGS -DOMPI_SKIP_MPICXX"
      AC_SUBST([MPI_CFLAGS])
      AC_MSG_NOTICE([MPI is enabled])
      HAVE_MPI=1
   ],
   [
      AC_MSG_ERROR([MPI was not found])
   ])
])
AM_CONDITIONAL([ENABLE_MPI], [test "x$HAVE_MPI" = x1])

# CUDA configuration
HAVE_CUDA=0
CUDA_CFLAGS=""
AC_ARG_WITH([cuda], [AS_HELP_STRING([--with-cuda=<version>], [Specify CUDA version])], [], [])
MODULE_NAME=""
AS_IF([test "x$with_cuda" == "x"], [
   AC_MSG_NOTICE([Building without CUDA])
],
[
   MODULE_NAME="cudart-$with_cuda"
   PKG_CHECK_MODULES([CUDA], [$MODULE_NAME], [
      cuda_prefix=$(pkg-config --variable=cudaroot "$MODULE_NAME")
      AC_MSG_CHECKING([for nvcc in $cuda_prefix/bin])
      AS_IF([test -x "$cuda_prefix/bin/nvcc"], [
         AC_MSG_RESULT([found])
         AC_DEFINE([USE_CUDA], [1], [Using CUDA])
         HAVE_CUDA=1

# CUDA requires replacing C++ compiler ($CXX) and adding switches to fix common errors.
         CUDA_CFLAGS="-I$cuda_prefix/include -x=cu --expt-relaxed-constexpr"
         AC_SUBST([CUDA_CFLAGS])
         AC_SUBST([NVCC], [$cuda_prefix/bin/nvcc])
         CXX=$NVCC
         AC_MSG_NOTICE([CUDA is enabled])
      ],
      [
         AC_MSG_RESULT([not found])
      ])
   ],
   [
      AC_MSG_ERROR([CUDA is disabled])
   ])
])
AM_CONDITIONAL([ENABLE_CUDA], [test "x$HAVE_CUDA" = x1])

# CUDA does not link properly, revert to g++. The makefiles should override per target $LINK to remove the compile-only switches.
AC_CHECK_PROG([CXXLD], [g++], [g++], [false])
AC_SUBST([CXXLD])

# GSL configuration
HAVE_GSL=0
AC_ARG_ENABLE([gsl], [AS_HELP_STRING([--enable-gsl], [enable GSL [default=no]])], [], [])
MODULE_NAME=""
AS_IF([test "x$enable_gsl" == "x"], [
   AC_MSG_NOTICE([Building without GSL])
],
[
   MODULE_NAME="gsl"
   PKG_CHECK_MODULES([GSL], [$MODULE_NAME])
   AS_IF([test "x$MODULE_NAME" != "x"], [
      AC_DEFINE([USE_GSL], [1], [Using GSL])
      AC_MSG_NOTICE([GSL is enabled])
      HAVE_GSL=1
   ],
   [
      AC_MSG_NOTICE([GSL was not found])
   ])
])
AM_CONDITIONAL([ENABLE_GSL], [test "x$HAVE_GSL" = x1])

# Check for SILO
HAVE_SILO=0
AC_ARG_ENABLE([silo], [AS_HELP_STRING([--enable-silo], [enable SILO [default=no]])], [], [])
AS_IF([test "x$enable_silo" == "x"], [
   AC_MSG_NOTICE([Building without SILO])
],
[
   AC_SEARCH_LIBS([DBPutQuadmesh], [silo siloh5], [
      AC_DEFINE([USE_SILO], [1], [Using SILO])
      AC_MSG_NOTICE([SILO is enabled])
      HAVE_SILO=1
   ],
   [
      AC_MSG_WARN([SILO was not found])
   ])
])
AM_CONDITIONAL([ENABLE_SILO], [test "x$HAVE_SILO" = x1])

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Everything below and up to the dashed line should go into configuration files
#-----------------------------------------------------------------------------------------------------------------------------------------------------

# Check for SLURM (not required)
#AC_SEARCH_LIBS([slurm_get_rem_time], [slurm slurmfull], [AC_DEFINE([USE_SLURM], [1], [Using SLURM])], [AC_MSG_WARN([SLURM library was not found.])])

# Define simulation execution types
AC_DEFINE([EXEC_SERIAL], [0], [Serial execution])
AC_DEFINE([EXEC_PARALLEL], [1], [Parallel execution])

# Set up simulation execution type
AC_ARG_WITH([execution], [AS_HELP_STRING([--with-execution=EXECUTION], [use EXECUTION=SERIAL|PARALLEL])], [], [])
AS_IF([test "x$with_execution" == "x"], [
   AC_MSG_ERROR([A value of EXECUTION is required])
],
[
   AS_IF([test "$with_execution" == "SERIAL" -o "$with_execution" == "PARALLEL"], [
      AC_DEFINE_UNQUOTED([EXEC_TYPE], [EXEC_$with_execution], [Choice of simulation execution type])
   ],
   [
      AC_MSG_ERROR([Invalid EXECUTION value])
   ])
])
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
AC_DEFINE([TRAJ_PARKER_SOURCE], [207], [Isotropic model (p,-,-) with source term])

# Set up the trajectory type
AC_ARG_WITH([trajectory], [AS_HELP_STRING([--with-trajectory=TRAJECTORY], [use TRAJECTORY=FIELDLINE|LORENTZ|GUIDING|GUIDING_SCATT|GUIDING_DIFF|GUIDING_DIFF_SCATT|FOCUSED|PARKER|PARKER_SOURCE])], [], [])
AS_IF([test "x$with_trajectory" == "x"], [
   AC_MSG_ERROR([A value of TRAJECTORY is required])
],
[
   AS_IF([test "$with_trajectory" == "FIELDLINE" -o "$with_trajectory" == "LORENTZ" -o "$with_trajectory" == "GUIDING" -o "$with_trajectory" == "GUIDING_SCATT" \
            -o "$with_trajectory" == "GUIDING_DIFF" -o "$with_trajectory" == "GUIDING_DIFF_SCATT" -o "$with_trajectory" == "FOCUSED" -o "$with_trajectory" == "PARKER" -o "$with_trajectory" == "PARKER_SOURCE"], [
      AC_DEFINE_UNQUOTED([TRAJ_TYPE], [TRAJ_$with_trajectory], [Choice of the trajectory physics])
   ],
   [
      AC_MSG_ERROR([Invalid TRAJECTORY value])
   ])
])
AC_MSG_NOTICE([Using "$with_trajectory" as the trajectory model])

# Define trajectory time flow direction types
AC_DEFINE([TRAJ_TIME_FLOW_FORWARD], [0], [Forward-in-time simulation])
AC_DEFINE([TRAJ_TIME_FLOW_BACKWARD], [1], [Backward-in-time simulation])

# Set up trajectory time flow direction
AC_ARG_WITH([time_flow], [AS_HELP_STRING([--with-time_flow=TIME_FLOW], [use TIME_FLOW=FORWARD|BACKWARD])], [], [])
AS_IF([test "x$with_time_flow" == "x"], [
   AC_MSG_ERROR([A value of TIME_FLOW is required])
],
[
   AS_IF([test "$with_time_flow" == "FORWARD" -o "$with_time_flow" == "BACKWARD"], [
      AC_DEFINE_UNQUOTED([TRAJ_TIME_FLOW], [TRAJ_TIME_FLOW_$with_time_flow], [Choice of trajectory time flow direction])
   ],
   [
      AC_MSG_ERROR([Invalid TIME_FLOW value])
   ])
])
AC_MSG_NOTICE([Using "$with_time_flow" as the trajectory time flow direction])

# Set up the RK integrator
AC_ARG_WITH([rkmethod], [AS_HELP_STRING([--with-rkmethod=RKMETHOD], [use RKMETHOD=0..29])], [], [])
AS_IF([test "x$with_rkmethod" == "x"], [
   AC_MSG_ERROR([A value of RKMETHOD is required])
],
[
   AS_IF([test "$with_rkmethod" -ge 0 -a "$with_rkmethod" -le 29], [
      AC_DEFINE_UNQUOTED([RK_INTEGRATOR_TYPE], [$with_rkmethod], [Choice of the Runge-Kutta method])
   ],
   [
      AC_MSG_ERROR([RKMETHOD must be between 0 and 29])
   ])
])
AC_MSG_NOTICE([Using "$with_rkmethod" as the Runge-Kutta method])

# Define server types
AC_DEFINE([SERVER_SELF], [299], [No server])
AC_DEFINE([SERVER_CARTESIAN], [300], [Cartesian server with uniform grid])
AC_DEFINE([SERVER_BATL], [301], [Cartesian server with AMR (part of BATS-R-US)])

# Set up the server type
AC_ARG_WITH([server], [AS_HELP_STRING([--with-server=SERVER], [use SERVER=SELF|CARTESIAN|BATL])], [], [])
AS_IF([test "x$with_server" == "x"], [
   AC_MSG_ERROR([A value of SERVER is required])
],
[
   AS_IF([test "$with_server" == "SELF" -o "$with_server" == "CARTESIAN" -o "$with_server" == "BATL"], [
      AC_DEFINE_UNQUOTED([SERVER_TYPE], [SERVER_$with_server], [Choice of the background server])
   ],
   [
      AC_MSG_ERROR([Invalid SERVER value])
   ])
])
AS_IF([test "$with_server" != "SELF" -a "$with_execution" == "SERIAL"], [
   AC_MSG_ERROR([SERIAL execution can only used with a SELF server type])
],
[
   AC_MSG_NOTICE([Using "$with_server" as the background server])
])

# Set up server interpolation order
AC_ARG_WITH([server_interp_order], [AS_HELP_STRING([--with-server_interp_order=SERVER_INTERP_ORDER], [use SERVER_INTERP_ORDER=-1|0|1])], [], [])
AS_IF([test "x$with_server_interp_order" == "x" -a $with_server != "SELF"], [
   AC_MSG_ERROR([A value of SERVER_INTERP_ORDER is required for all server types other than SELF])
],
[
   AS_IF([test "x$with_server_interp_order" == "x" -o "$with_server_interp_order" == "-1" -o "$with_server_interp_order" == "0" -o "$with_server_interp_order" == "1"], [
      AC_DEFINE_UNQUOTED([SERVER_INTERP_ORDER], [$with_server_interp_order], [Choice of background server interpolation order])
   ],
   [
      AC_MSG_ERROR([Invalid SERVER_INTERP_ORDER value])
   ])
])
AC_MSG_NOTICE([Using "$with_server_interp_order" as the background server interpolation order])

# Set up server number of ghost cells per side of blocks
AC_ARG_WITH([server_num_gcs], [AS_HELP_STRING([--with-server_num_gcs=SERVER_NUM_GHOST_CELLS], [use SERVER_NUM_GHOST_CELLS=0|1])], [], [])
AS_IF([test "x$with_server_num_gcs" == "x" -a $with_server != "SELF"], [
   AC_MSG_ERROR([A value of SERVER_NUM_GHOST_CELLS is required for all server types other than SELF])
],
[
   AS_IF([test "x$with_server_num_gcs" == "x" -o "$with_server_num_gcs" == "0" -o "$with_server_num_gcs" == "1"], [
      AC_DEFINE_UNQUOTED([SERVER_NUM_GHOST_CELLS], [$with_server_num_gcs], [Choice of background server number of ghost cells per side of block])
   ],
   [
      AC_MSG_ERROR([Invalid SERVER_NUM_GHOST_CELLS value])
   ])
])
AC_MSG_NOTICE([Using "$with_server_num_gcs" as the background server number of ghost cells per side of block])

#-----------------------------------------------------------------------------------------------------------------------------------------------------

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T

# Checks for library functions.
AC_CHECK_FUNCS([floor memset pow sqrt])

# Files to generate
AC_CONFIG_FILES([Makefile runs/Makefile tests/Makefile benchmarks/Makefile])

AC_OUTPUT

