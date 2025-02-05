## Process this file with automake to produce Makefile.in

AM_CXXFLAGS = $(MPI_CFLAGS)
AM_LDFLAGS = $(MPI_LIBS)

bin_PROGRAMS = gcr_modulation_example postprocess_gcr_modulation_results

SPBL_COMMON_DIR = ../common
SPBL_SOURCE_DIR = ../src

gcr_modulation_example_SOURCES = gcr_modulation_example.cc \
   $(SPBL_SOURCE_DIR)/simulation.cc \
   $(SPBL_SOURCE_DIR)/simulation.hh \
   $(SPBL_SOURCE_DIR)/trajectory_parker.cc \
   $(SPBL_SOURCE_DIR)/trajectory_parker.hh \
   $(SPBL_SOURCE_DIR)/trajectory_base.cc \
   $(SPBL_SOURCE_DIR)/trajectory_base.hh \
   $(SPBL_SOURCE_DIR)/background_solarwind.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/boundary_momentum.cc \
   $(SPBL_SOURCE_DIR)/boundary_momentum.hh \
   $(SPBL_SOURCE_DIR)/boundary_space.cc \
   $(SPBL_SOURCE_DIR)/boundary_space.hh \
   $(SPBL_SOURCE_DIR)/boundary_time.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.hh \
   $(SPBL_SOURCE_DIR)/boundary_base.cc \
   $(SPBL_SOURCE_DIR)/boundary_base.hh \
   $(SPBL_SOURCE_DIR)/initial_time.cc \
   $(SPBL_SOURCE_DIR)/initial_time.hh \
   $(SPBL_SOURCE_DIR)/initial_space.cc \
   $(SPBL_SOURCE_DIR)/initial_space.hh \
   $(SPBL_SOURCE_DIR)/initial_momentum.cc \
   $(SPBL_SOURCE_DIR)/initial_momentum.hh \
   $(SPBL_SOURCE_DIR)/initial_base.cc \
   $(SPBL_SOURCE_DIR)/initial_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_other.cc \
   $(SPBL_SOURCE_DIR)/diffusion_other.hh \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_SOURCE_DIR)/distribution_other.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.hh \
   $(SPBL_SOURCE_DIR)/distribution_templated.cc \
   $(SPBL_SOURCE_DIR)/distribution_templated.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_SOURCE_DIR)/traj_config.hh \
   $(SPBL_SOURCE_DIR)/server_config.hh \
   $(SPBL_COMMON_DIR)/workload_manager.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

gcr_modulation_example_LDADD = $(MPI_LIBS) $(GSL_LIBS)

postprocess_gcr_modulation_results_SOURCES = postprocess_gcr_modulation_results.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

postprocess_gcr_modulation_results_LDADD = $(MPI_LIBS) $(GSL_LIBS)