## Process this file with automake to produce Makefile.in

AM_CXXFLAGS = $(MPI_CFLAGS) -DGEO_DEBUG -Wall -Wno-comment
#AM_CXXFLAGS = $(MPI_CFLAGS)
#AM_LDFLAGS = $(MPI_LIBS)
#BATL_HOME = /home/vladimir/Science/Global3D/Michigan/BATL-master

AM_FCFLAGS = -J$(BATL_HOME)/share/include

bin_PROGRAMS = main_test_mpi_config \
               main_test_data_container \
               main_test_physics \
               main_test_background_bochum \
               main_test_drawable_tesselation \
               main_test_drawable_sector \
               main_test_exchange_site \
               main_test_shared_site \
               main_test_distance_map \
               main_test_polygonal_addressing \
               main_test_geodesic_sector \
               main_test_grid_block \
               main_test_buffered_block \
               main_test_domain_partition

SPBL_COMMON_DIR = ../common
SPBL_SOURCE_DIR = ../src
SPBL_GEODESIC_DIR = ../geodesic
SPBL_FLUID_DIR = ../fluid
SPBL_GEOMETRY_DIR = ../geometry

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_mpi_config_SOURCES = main_test_mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh

main_test_mpi_config_LDADD = $(MPI_LIBS)

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_data_container_SOURCES = main_test_data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh   

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_physics_SOURCES = main_test_physics.cc \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_background_bochum_SOURCES = main_test_background_bochum.cc \
   $(SPBL_SOURCE_DIR)/background_vlism_bochum.cc \
   $(SPBL_SOURCE_DIR)/background_vlism_bochum.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh

main_test_background_bochum_LDADD = $(GSL_LIBS)

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_drawable_tesselation_SOURCES = main_test_drawable_tesselation.cc \
   $(SPBL_GEODESIC_DIR)/drawable_tesselation.cc \
   $(SPBL_GEODESIC_DIR)/drawable_tesselation.hh \
   $(SPBL_GEODESIC_DIR)/spherical_tesselation.cc \
   $(SPBL_GEODESIC_DIR)/spherical_tesselation.hh \
   $(SPBL_GEODESIC_DIR)/polyhedron.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_drawable_sector_SOURCES = main_test_drawable_sector.cc \
   $(SPBL_GEODESIC_DIR)/drawable_sector.cc \
   $(SPBL_GEODESIC_DIR)/drawable_sector.hh \
   $(SPBL_GEODESIC_DIR)/geodesic_sector.cc \
   $(SPBL_GEODESIC_DIR)/geodesic_sector.hh \
   $(SPBL_GEODESIC_DIR)/polygonal_addressing.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_exchange_site_SOURCES = main_test_exchange_site.cc \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/exchange_site.hh \
   $(SPBL_COMMON_DIR)/communication_site.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

main_test_exchange_site_LDADD = $(MPI_LIBS)

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_shared_site_SOURCES = main_test_shared_site.cc \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/shared_site.hh \
   $(SPBL_COMMON_DIR)/communication_site.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

main_test_shared_site_LDADD = $(MPI_LIBS)

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_distance_map_SOURCES = main_test_distance_map.cc \
   $(SPBL_GEOMETRY_DIR)/distance_map.cc \
   $(SPBL_GEOMETRY_DIR)/distance_map.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh \
   $(SPBL_COMMON_DIR)/random.hh

main_test_distance_map_LDADD = $(GSL_LIBS)

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_polygonal_addressing_SOURCES = main_test_polygonal_addressing.cc \
   $(SPBL_GEODESIC_DIR)/polygonal_addressing.hh

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_geodesic_sector_SOURCES = main_test_geodesic_sector.cc \
   $(SPBL_GEODESIC_DIR)/geodesic_sector.cc \
   $(SPBL_GEODESIC_DIR)/geodesic_sector.hh \
   $(SPBL_GEODESIC_DIR)/polygonal_addressing.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_grid_block_SOURCES = main_test_grid_block.cc \
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
   $(SPBL_GEOMETRY_DIR)/distance_map.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh

main_test_grid_block_LDADD = $(GSL_LIBS)

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_buffered_block_SOURCES = main_test_buffered_block.cc \
   $(SPBL_GEODESIC_DIR)/buffered_block.cc \
   $(SPBL_GEODESIC_DIR)/buffered_block.hh \
   $(SPBL_GEODESIC_DIR)/stenciled_block.cc \
   $(SPBL_GEODESIC_DIR)/stenciled_block.hh \
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
   $(SPBL_GEOMETRY_DIR)/distance_map.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/communication_site.hh \
   $(SPBL_COMMON_DIR)/exchange_site.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh

main_test_buffered_block_LDADD = $(MPI_LIBS) $(GSL_LIBS)

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_domain_partition_SOURCES = main_test_domain_partition.cc \
   $(SPBL_GEODESIC_DIR)/domain_partition.cc \
   $(SPBL_GEODESIC_DIR)/domain_partition.hh \
   $(SPBL_GEODESIC_DIR)/buffered_block.cc \
   $(SPBL_GEODESIC_DIR)/buffered_block.hh \
   $(SPBL_GEODESIC_DIR)/stenciled_block.cc \
   $(SPBL_GEODESIC_DIR)/stenciled_block.hh \
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
   $(SPBL_GEOMETRY_DIR)/distance_map.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/communication_site.hh \
   $(SPBL_COMMON_DIR)/exchange_site.hh \
   $(SPBL_COMMON_DIR)/shared_site.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh

main_test_domain_partition_LDADD = $(MPI_LIBS) $(GSL_LIBS)

