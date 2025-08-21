## Process this file with automake to produce Makefile.in

## **** spdata-fields update ****:
## The build is the same for all the benchmarks.
## This is done on purpose in case this should change while the spdata update is in progress.

#AM_CXXFLAGS = $(MPI_CFLAGS) -DGEO_DEBUG
AM_CXXFLAGS = $(MPI_CFLAGS)
AM_LDFLAGS = $(MPI_LIBS)
AM_FCFLAGS = -J$(BATL_INCL_DIR)

bin_PROGRAMS = main_test_dipole_visualization \
               main_test_turb_waves \
               main_test_parker_spiral \
               main_test_init_cond_records \
               main_test_dipole_periods \
               main_test_pa_distro_isotrop \
               main_test_pa_scatt \
               main_test_perp_diff main_test_full_diff \
               main_test_modulation_cartesian_parker \
               main_postprocess_modulation_cartesian_parker \
               main_generate_cartesian_solarwind_background

SPBL_COMMON_DIR = ../common
SPBL_SOURCE_DIR = ../src

SPBL_COMMON_SOURCES = \
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

main_test_dipole_visualization_SOURCES = main_test_dipole_visualization.cc $(SPBL_COMMON_SOURCES)
main_test_dipole_visualization_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_test_turb_waves_SOURCES = main_test_turb_waves.cc $(SPBL_COMMON_SOURCES)
main_test_turb_waves_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_test_parker_spiral_SOURCES = main_test_parker_spiral.cc $(SPBL_COMMON_SOURCES)
main_test_parker_spiral_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_test_init_cond_records_SOURCES = main_test_init_cond_records.cc $(SPBL_COMMON_SOURCES)
main_test_init_cond_records_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_test_dipole_periods_SOURCES = main_test_dipole_periods.cc $(SPBL_COMMON_SOURCES)
main_test_dipole_periods_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_test_pa_distro_isotrop_SOURCES = main_test_pa_distro_isotrop.cc $(SPBL_COMMON_SOURCES)
main_test_pa_distro_isotrop_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_test_pa_scatt_SOURCES = main_test_pa_scatt.cc $(SPBL_COMMON_SOURCES)
main_test_pa_scatt_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_test_perp_diff_SOURCES = main_test_perp_diff.cc $(SPBL_COMMON_SOURCES)
main_test_perp_diff_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_test_full_diff_SOURCES = main_test_full_diff.cc $(SPBL_COMMON_SOURCES)
main_test_full_diff_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_test_modulation_cartesian_parker_SOURCES = main_test_modulation_cartesian_parker.cc $(SPBL_COMMON_SOURCES)
main_test_modulation_cartesian_parker_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_postprocess_modulation_cartesian_parker_SOURCES = main_postprocess_modulation_cartesian_parker.cc $(SPBL_COMMON_SOURCES)
main_postprocess_modulation_cartesian_parker_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_generate_cartesian_solarwind_background_SOURCES = main_generate_cartesian_solarwind_background.cc $(SPBL_COMMON_SOURCES)
main_generate_cartesian_solarwind_background_LDADD = $(MPI_LIBS) $(GSL_LIBS)
