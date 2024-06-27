## Process this file with automake to produce Makefile.in

AM_CXXFLAGS = $(MPI_CFLAGS)
AM_LDFLAGS = $(MPI_LIBS)

bin_PROGRAMS = main_test_units \
               main_tesstest \
               main_drawsector \
               main_onegridblock \
               main_test_exchange_site

SPBL_COMMON_DIR = ../common
SPBL_SOURCE_DIR = ../src
SPBL_GEODESIC_DIR = ../geodesic
SPBL_FLUID_DIR = ../fluid
SPBL_GEOMETRY_DIR = ../geometry


main_test_units_SOURCES = main_test_units.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/physics.cc

main_test_units_LDADD = $(MPI_LIBS) $(GSL_LIBS)


main_tesstest_SOURCES = main_tesstest.cc \
   $(SPBL_GEODESIC_DIR)/drawable_tesselation.cc \
   $(SPBL_GEODESIC_DIR)/drawable_tesselation.hh \
   $(SPBL_GEODESIC_DIR)/spherical_tesselation.cc \
   $(SPBL_GEODESIC_DIR)/spherical_tesselation.hh \
   $(SPBL_GEODESIC_DIR)/polyhedron.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh

main_tesstest_LDADD = 


main_drawsector_SOURCES = main_drawsector.cc \
   $(SPBL_GEODESIC_DIR)/drawable_sector.cc \
   $(SPBL_GEODESIC_DIR)/drawable_sector.hh \
   $(SPBL_GEODESIC_DIR)/geodesic_sector.cc \
   $(SPBL_GEODESIC_DIR)/geodesic_sector.hh \
   $(SPBL_GEODESIC_DIR)/polygonal_addressing.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh

main_drawsector_LDADD = 


main_onegridblock_SOURCES = main_onegridblock.cc \
   $(SPBL_GEODESIC_DIR)/grid_block.cc \
   $(SPBL_GEODESIC_DIR)/grid_block.hh \
   $(SPBL_GEODESIC_DIR)/geodesic_sector.cc \
   $(SPBL_GEODESIC_DIR)/geodesic_sector.hh \
   $(SPBL_GEODESIC_DIR)/polygonal_addressing.hh \
   $(SPBL_GEODESIC_DIR)/spherical_slab.hh \
   $(SPBL_GEODESIC_DIR)/traversable_tesselation.cc \
   $(SPBL_GEODESIC_DIR)/traversable_tesselation.hh \
   $(SPBL_GEODESIC_DIR)/requestable_tesselation.cc \
   $(SPBL_GEODESIC_DIR)/requestable_tesselation.hh \
   $(SPBL_GEODESIC_DIR)/spherical_tesselation.cc \
   $(SPBL_GEODESIC_DIR)/spherical_tesselation.hh \
   $(SPBL_GEODESIC_DIR)/polyhedron.hh \
   $(SPBL_GEOMETRY_DIR)/distance_map.cc \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh

main_onegridblock_LDADD = $(SILO_LIBS)


main_test_exchange_site_SOURCES = main_test_exchange_site.cc \
   $(SPBL_COMMON_DIR)/exchange_site.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

main_test_exchange_site_LDADD =
