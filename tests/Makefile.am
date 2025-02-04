## Process this file with automake to produce Makefile.in

AM_CXXFLAGS = $(MPI_CFLAGS) -DGEO_DEBUG
#AM_CXXFLAGS = $(MPI_CFLAGS)
AM_LDFLAGS = $(MPI_LIBS)
BATL_HOME = /home/vladimir/Science/Global3D/Michigan/BATL-master

AM_FCFLAGS = -J$(BATL_HOME)/share/include

bin_PROGRAMS = main_test_trajectory1 \
               main_test_trajectory_mswave \
               main_test_simulation1 \
               main_test_simulation2 \
               main_test_mpi_config \
               main_test_units \
               main_plot_backgrounds \
               main_test_background_bochum \
               main_test_background_waves \
               main_test_background_mswave \
               main_test_background_batl \
               main_plot_solarwind_batl \
               main_tesstest \
               main_simulation_lee \
               main_testread \
               main_test_trajectory_batl \
               main_test_alfven_turb \
               main_tesstest \
               main_drawsector \
               main_test_exchange_site \
               main_test_shared_site \
               main_test_grid_block \
               main_test_buffered_block

SPBL_COMMON_DIR = ../common
SPBL_SOURCE_DIR = ../src
SPBL_GEODESIC_DIR = ../geodesic
SPBL_FLUID_DIR = ../fluid
SPBL_GEOMETRY_DIR = ../geometry

main_test_trajectory1_SOURCES = main_test_trajectory1.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_SOURCE_DIR)/trajectory_guiding_scatt.cc \
   $(SPBL_SOURCE_DIR)/trajectory_guiding_scatt.hh \
   $(SPBL_SOURCE_DIR)/trajectory_guiding.cc \
   $(SPBL_SOURCE_DIR)/trajectory_guiding.hh \
   $(SPBL_SOURCE_DIR)/trajectory_base.cc \
   $(SPBL_SOURCE_DIR)/trajectory_base.hh \
   $(SPBL_SOURCE_DIR)/background_uniform.cc \
   $(SPBL_SOURCE_DIR)/background_uniform.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.hh \
   $(SPBL_SOURCE_DIR)/boundary_base.cc \
   $(SPBL_SOURCE_DIR)/boundary_base.hh \
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
   $(SPBL_SOURCE_DIR)/traj_config.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

main_test_trajectory1_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_test_trajectory_mswave_SOURCES = main_test_trajectory_mswave.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_SOURCE_DIR)/trajectory_lorentz.cc \
   $(SPBL_SOURCE_DIR)/trajectory_lorentz.hh \
   $(SPBL_SOURCE_DIR)/trajectory_base.cc \
   $(SPBL_SOURCE_DIR)/trajectory_base.hh \
   $(SPBL_SOURCE_DIR)/background_mswave.cc \
   $(SPBL_SOURCE_DIR)/background_mswave.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_SOURCE_DIR)/boundary_time.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.hh \
   $(SPBL_SOURCE_DIR)/boundary_base.cc \
   $(SPBL_SOURCE_DIR)/boundary_base.hh \
   $(SPBL_SOURCE_DIR)/initial_space.cc \
   $(SPBL_SOURCE_DIR)/initial_space.hh \
   $(SPBL_SOURCE_DIR)/initial_momentum.cc \
   $(SPBL_SOURCE_DIR)/initial_momentum.hh \
   $(SPBL_SOURCE_DIR)/initial_base.cc \
   $(SPBL_SOURCE_DIR)/initial_base.hh \
   $(SPBL_SOURCE_DIR)/traj_config.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

main_test_trajectory_mswave_LDADD = $(GSL_LIBS)

main_test_simulation1_SOURCES = main_test_simulation1.cc \
   $(SPBL_SOURCE_DIR)/simulation.cc \
   $(SPBL_SOURCE_DIR)/simulation.hh \
   $(SPBL_SOURCE_DIR)/distribution_other.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_SOURCE_DIR)/trajectory_guiding_scatt.cc \
   $(SPBL_SOURCE_DIR)/trajectory_guiding_scatt.hh \
   $(SPBL_SOURCE_DIR)/trajectory_guiding.cc \
   $(SPBL_SOURCE_DIR)/trajectory_guiding.hh \
   $(SPBL_SOURCE_DIR)/trajectory_base.cc \
   $(SPBL_SOURCE_DIR)/trajectory_base.hh \
   $(SPBL_SOURCE_DIR)/background_uniform.cc \
   $(SPBL_SOURCE_DIR)/background_uniform.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.hh \
   $(SPBL_SOURCE_DIR)/boundary_space.cc \
   $(SPBL_SOURCE_DIR)/boundary_space.hh \
   $(SPBL_SOURCE_DIR)/boundary_base.cc \
   $(SPBL_SOURCE_DIR)/boundary_base.hh \
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
   $(SPBL_SOURCE_DIR)/traj_config.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

main_test_simulation1_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_test_simulation2_SOURCES = main_test_simulation2.cc \
   $(SPBL_SOURCE_DIR)/simulation.cc \
   $(SPBL_SOURCE_DIR)/simulation.hh \
   $(SPBL_SOURCE_DIR)/distribution_other.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.hh \
   $(SPBL_SOURCE_DIR)/distribution_templated.cc \
   $(SPBL_SOURCE_DIR)/distribution_templated.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_SOURCE_DIR)/trajectory_guiding_diff_scatt.cc \
   $(SPBL_SOURCE_DIR)/trajectory_guiding_diff_scatt.hh \
   $(SPBL_SOURCE_DIR)/trajectory_guiding_diff.cc \
   $(SPBL_SOURCE_DIR)/trajectory_guiding_diff.hh \
   $(SPBL_SOURCE_DIR)/trajectory_guiding.cc \
   $(SPBL_SOURCE_DIR)/trajectory_guiding.hh \
   $(SPBL_SOURCE_DIR)/trajectory_base.cc \
   $(SPBL_SOURCE_DIR)/trajectory_base.hh \
   $(SPBL_SOURCE_DIR)/background_uniform.cc \
   $(SPBL_SOURCE_DIR)/background_uniform.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.hh \
   $(SPBL_SOURCE_DIR)/boundary_space.cc \
   $(SPBL_SOURCE_DIR)/boundary_space.hh \
   $(SPBL_SOURCE_DIR)/boundary_base.cc \
   $(SPBL_SOURCE_DIR)/boundary_base.hh \
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
   $(SPBL_SOURCE_DIR)/traj_config.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

main_test_simulation2_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_test_mpi_config_SOURCES = main_test_mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh

main_test_mpi_config_LDADD = $(MPI_LIBS)

main_test_units_SOURCES = main_test_units.cc \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

main_plot_backgrounds_SOURCES = main_plot_backgrounds.cc \
   $(SPBL_SOURCE_DIR)/background_vlism_bochum.cc \
   $(SPBL_SOURCE_DIR)/background_vlism_bochum.hh \
   $(SPBL_SOURCE_DIR)/background_dipole.cc \
   $(SPBL_SOURCE_DIR)/background_dipole.hh \
   $(SPBL_SOURCE_DIR)/background_cylindrical_obstacle.cc \
   $(SPBL_SOURCE_DIR)/background_cylindrical_obstacle.hh \
   $(SPBL_SOURCE_DIR)/background_magnetized_cylinder.cc \
   $(SPBL_SOURCE_DIR)/background_magnetized_cylinder.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/multi_index.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/definitions.cc \
   $(SPBL_COMMON_DIR)/definitions.hh

main_plot_backgrounds_LDADD = $(GSL_LIBS)

main_test_background_bochum_SOURCES = main_test_background_bochum.cc \
   $(SPBL_SOURCE_DIR)/background_vlism_bochum.cc \
   $(SPBL_SOURCE_DIR)/background_vlism_bochum.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

main_test_background_bochum_LDADD = $(GSL_LIBS)

main_test_background_batl_SOURCES = main_test_background_batl.cc \
   $(SPBL_SOURCE_DIR)/background_batl.cc \
   $(SPBL_SOURCE_DIR)/background_batl.hh \
   $(SPBL_SOURCE_DIR)/background_cartesian.cc \
   $(SPBL_SOURCE_DIR)/background_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/server_batl.cc \
   $(SPBL_SOURCE_DIR)/server_batl.hh \
   $(SPBL_SOURCE_DIR)/server_cartesian.cc \
   $(SPBL_SOURCE_DIR)/server_cartesian.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_SOURCE_DIR)/block_batl.cc \
   $(SPBL_SOURCE_DIR)/block_batl.hh \
   $(SPBL_SOURCE_DIR)/block_cartesian.cc \
   $(SPBL_SOURCE_DIR)/block_cartesian.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/reader_cartesian.cc \
   $(SPBL_SOURCE_DIR)/reader_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/spectrum_interface.f90 \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

main_test_background_batl_LDADD = -L$(BATL_HOME)/lib -lREADAMR -lTIMING -lSHARE -lgfortran -lmpi_mpifh

main_plot_solarwind_batl_SOURCES = main_plot_solarwind_batl.cc \
   $(SPBL_SOURCE_DIR)/background_batl.cc \
   $(SPBL_SOURCE_DIR)/background_batl.hh \
   $(SPBL_SOURCE_DIR)/background_cartesian.cc \
   $(SPBL_SOURCE_DIR)/background_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_server.cc \
   $(SPBL_SOURCE_DIR)/background_server.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/server_batl.cc \
   $(SPBL_SOURCE_DIR)/server_batl.hh \
   $(SPBL_SOURCE_DIR)/server_cartesian.cc \
   $(SPBL_SOURCE_DIR)/server_cartesian.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_SOURCE_DIR)/block_batl.cc \
   $(SPBL_SOURCE_DIR)/block_batl.hh \
   $(SPBL_SOURCE_DIR)/block_cartesian.cc \
   $(SPBL_SOURCE_DIR)/block_cartesian.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/reader_cartesian.cc \
   $(SPBL_SOURCE_DIR)/reader_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/spectrum_interface.f90 \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

main_plot_solarwind_batl_LDADD = -L$(BATL_HOME)/lib -lREADAMR -lTIMING -lSHARE -lgfortran -lmpi_mpifh

main_test_trajectory_batl_SOURCES = main_test_trajectory_batl.cc \
   $(SPBL_SOURCE_DIR)/distribution_templated.cc \
   $(SPBL_SOURCE_DIR)/distribution_templated.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_SOURCE_DIR)/trajectory_lorentz.cc \
   $(SPBL_SOURCE_DIR)/trajectory_lorentz.hh \
   $(SPBL_SOURCE_DIR)/trajectory_guiding.cc \
   $(SPBL_SOURCE_DIR)/trajectory_guiding.hh \
   $(SPBL_SOURCE_DIR)/trajectory_fieldline.cc \
   $(SPBL_SOURCE_DIR)/trajectory_fieldline.hh \
   $(SPBL_SOURCE_DIR)/trajectory_base.cc \
   $(SPBL_SOURCE_DIR)/trajectory_base.hh \
   $(SPBL_SOURCE_DIR)/background_uniform.cc \
   $(SPBL_SOURCE_DIR)/background_uniform.hh \
   $(SPBL_SOURCE_DIR)/background_batl.cc \
   $(SPBL_SOURCE_DIR)/background_batl.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/server_batl.cc \
   $(SPBL_SOURCE_DIR)/server_batl.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/block_batl.cc \
   $(SPBL_SOURCE_DIR)/block_batl.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.hh \
   $(SPBL_SOURCE_DIR)/boundary_space.cc \
   $(SPBL_SOURCE_DIR)/boundary_space.hh \
   $(SPBL_SOURCE_DIR)/boundary_base.cc \
   $(SPBL_SOURCE_DIR)/boundary_base.hh \
   $(SPBL_SOURCE_DIR)/initial_space.cc \
   $(SPBL_SOURCE_DIR)/initial_space.hh \
   $(SPBL_SOURCE_DIR)/initial_momentum.cc \
   $(SPBL_SOURCE_DIR)/initial_momentum.hh \
   $(SPBL_SOURCE_DIR)/initial_base.cc \
   $(SPBL_SOURCE_DIR)/initial_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_SOURCE_DIR)/spectrum_interface.f90 \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

main_test_trajectory_batl_LDADD = -L$(BATL_HOME)/lib -lREADAMR -lTIMING -lSHARE -lgfortran -lmpi_mpifh $(GSL_LIBS)

main_testread_SOURCES = main_testread.cc \
   $(SPBL_SOURCE_DIR)/spectrum_interface.f90

main_testread_LDADD = -L$(BATL_HOME)/lib -lREADAMR -lTIMING -lSHARE -lgfortran -lmpi_mpifh

main_test_background_waves_SOURCES = main_test_background_waves.cc \
   $(SPBL_SOURCE_DIR)/background_waves.cc \
   $(SPBL_SOURCE_DIR)/background_waves.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/multi_index.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/turb_prop.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/definitions.cc

main_test_background_waves_LDADD = $(GSL_LIBS)

main_test_background_mswave_SOURCES = main_test_background_mswave.cc \
   $(SPBL_SOURCE_DIR)/background_mswave.cc \
   $(SPBL_SOURCE_DIR)/background_mswave.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

main_test_background_mswave_LDADD = $(GSL_LIBS)

main_simulation_lee_SOURCES = main_simulation_lee.cc \
   $(SPBL_SOURCE_DIR)/simulation.cc \
   $(SPBL_SOURCE_DIR)/simulation.hh \
   $(SPBL_SOURCE_DIR)/distribution_other.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_SOURCE_DIR)/trajectory_lorentz.cc \
   $(SPBL_SOURCE_DIR)/trajectory_lorentz.hh \
   $(SPBL_SOURCE_DIR)/trajectory_base.cc \
   $(SPBL_SOURCE_DIR)/trajectory_base.hh \
   $(SPBL_SOURCE_DIR)/background_waves.cc \
   $(SPBL_SOURCE_DIR)/background_waves.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_SOURCE_DIR)/boundary_time.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.hh \
   $(SPBL_SOURCE_DIR)/boundary_base.cc \
   $(SPBL_SOURCE_DIR)/boundary_base.hh \
   $(SPBL_SOURCE_DIR)/initial_space.cc \
   $(SPBL_SOURCE_DIR)/initial_space.hh \
   $(SPBL_SOURCE_DIR)/initial_momentum.cc \
   $(SPBL_SOURCE_DIR)/initial_momentum.hh \
   $(SPBL_SOURCE_DIR)/initial_base.cc \
   $(SPBL_SOURCE_DIR)/initial_base.hh \
   $(SPBL_SOURCE_DIR)/traj_config.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

main_simulation_lee_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_test_alfven_turb_SOURCES = main_test_alfven_turb.cc \
   $(SPBL_SOURCE_DIR)/simulation.cc \
   $(SPBL_SOURCE_DIR)/simulation.hh \
   $(SPBL_SOURCE_DIR)/trajectory_lorentz.cc \
   $(SPBL_SOURCE_DIR)/trajectory_lorentz.hh \
   $(SPBL_SOURCE_DIR)/trajectory_base.cc \
   $(SPBL_SOURCE_DIR)/trajectory_base.hh \
   $(SPBL_SOURCE_DIR)/background_cartesian.cc \
   $(SPBL_SOURCE_DIR)/background_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_waves.cc \
   $(SPBL_SOURCE_DIR)/background_waves.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/server_cartesian.cc \
   $(SPBL_SOURCE_DIR)/server_cartesian.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_SOURCE_DIR)/block_cartesian.cc \
   $(SPBL_SOURCE_DIR)/block_cartesian.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/reader_cartesian.cc \
   $(SPBL_SOURCE_DIR)/reader_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/boundary_momentum.cc \
   $(SPBL_SOURCE_DIR)/boundary_momentum.hh \
   $(SPBL_SOURCE_DIR)/boundary_space.cc \
   $(SPBL_SOURCE_DIR)/boundary_space.hh \
   $(SPBL_SOURCE_DIR)/boundary_time.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.hh \
   $(SPBL_SOURCE_DIR)/boundary_base.cc \
   $(SPBL_SOURCE_DIR)/boundary_base.hh \
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
   $(SPBL_SOURCE_DIR)/server_exceptions.hh \
   $(SPBL_COMMON_DIR)/workload_manager.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/turb_prop.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/multi_index.cc \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/definitions.cc

main_test_alfven_turb_LDADD = $(MPI_LIBS) $(GSL_LIBS)


main_drawsector_SOURCES = main_drawsector.cc \
   $(SPBL_GEODESIC_DIR)/drawable_sector.cc \
   $(SPBL_GEODESIC_DIR)/drawable_sector.hh \
   $(SPBL_GEODESIC_DIR)/geodesic_sector.cc \
   $(SPBL_GEODESIC_DIR)/geodesic_sector.hh \
   $(SPBL_GEODESIC_DIR)/polygonal_addressing.hh \
   $(SPBL_GEODESIC_DIR)/polyhedron.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh


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

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_exchange_site_SOURCES = main_test_exchange_site.cc \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/exchange_site.hh \
   $(SPBL_COMMON_DIR)/communication_site.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

main_test_exchange_site_LDADD =

#-----------------------------------------------------------------------------------------------------------------------------------------------------

main_test_shared_site_SOURCES = main_test_shared_site.cc \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/shared_site.hh \
   $(SPBL_COMMON_DIR)/communication_site.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

main_test_shared_site_LDADD =

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
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

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
   $(SPBL_COMMON_DIR)/exchange_site.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/gpu_config.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

main_test_buffered_block_LDADD = $(GSL_LIBS)

