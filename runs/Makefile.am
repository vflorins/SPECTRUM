## Process this file with automake to produce Makefile.in

AM_CXXFLAGS = $(MPI_CFLAGS)
AM_LDFLAGS = $(MPI_LIBS)

bin_PROGRAMS = main_test_units

SPBL_COMMON_DIR = ../common
SPBL_SOURCE_DIR = ../src

main_test_units_SOURCES = main_test_units.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/physics.cc

main_test_units_LDADD = $(MPI_LIBS) $(GSL_LIBS)