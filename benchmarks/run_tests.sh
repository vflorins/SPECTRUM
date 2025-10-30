#!/bin/bash

# File to hold benchmark results
results_file="benchmark_results.txt"

# Parallelization parameters
n_cpus=$(grep "cpu cores" /proc/cpuinfo | uniq |  awk '{print $4}')
long_sim=100000
short_sim=10000
long_batch_size=1000
short_batch_size=100

# Flags to select which tests to run
# VISUAL TESTS
dipole_visual_test=false
turb_waves_test=false
parker_spiral_test=false
init_cond_record_test=false
# QUANTITATIVE TESTS
dipole_drifts_test=false
pa_distro_iso_test=false
pa_scatt_test=false
perp_diff_test=false
full_diff_test=false
diff_shock_acc_test=true
modulation_cartesian_parker_spiral=true

# Function to go up one directory and configure code
function configure {
	cd ..
	if test "$5" = "SELF"
	then
		./configure CXXFLAGS="-Ofast" --with-mpi=openmpi --with-execution=$1 \
			--with-trajectory=$2 --with-time_flow=$3 --with-rkmethod=$4 \
			--with-server=$5
	else
		./configure CXXFLAGS="-Ofast" --with-mpi=openmpi --with-execution=$1 \
			--with-trajectory=$2 --with-time_flow=$3 --with-rkmethod=$4 \
			--with-server=$5 --with-server_interp_order=$6 --with-server_num_gcs=$7
	fi
	cd -
}

# Function to make and run test
function make_and_run {
	make $2
	mkdir -p $1
	cd $1
	if [ $3 -eq 1 ]
	then
		../$2 >> ../$results_file
	else
		mpirun -np $3 ../$2 $4 $5 >> ../$results_file
	fi
	cd -
}

# Function to return exit code
function report_if_failed {
	if [ $1 -ne 0 ]
	then
	echo "" >> $results_file
	echo $2 >> $results_file
	echo "=========================================================" >> $results_file
	echo "***failed*** with error code ${1}" >> $results_file
	echo "=========================================================" >> $results_file
	echo "" >> $results_file
	fi
}

# Before running this script, execute the following in the root directory:
# 1) module load mpi
# 2) autoreconf
# 3) automake --add-missing
echo -n "Benchmark results " > $results_file
date >> $results_file
echo "------------------------------------------------------------" >> $results_file

# DIPOLE FIELD VISUALIZATION
if $dipole_visual_test
then	
	configure SERIAL LORENTZ FORWARD 0 SELF
	make_and_run output_data main_test_dipole_visualization 1
fi
report_if_failed $? "DIPOLE FIELD VISUALIZATION"

# TURBULENCE VIA SUPERPOSITION OF WAVES
if $turb_waves_test
then
	configure SERIAL FIELDLINE FORWARD 29 SELF
	make_and_run output_data main_test_turb_waves 1
fi
report_if_failed $? "TURBULENCE VIA SUPERPOSITION OF WAVES"

# PARKER SPIRAL SOLAR WIND
if $parker_spiral_test
then
	configure SERIAL FIELDLINE FORWARD 29 SELF
	make_and_run output_data main_test_parker_spiral 1
fi
report_if_failed $? "PARKER SPIRAL SOLAR WIND"

# INITIAL CONDITION RECORDS TEST
if $init_cond_record_test
then	
	configure PARALLEL FOCUSED FORWARD 29 SELF
	make_and_run output_data main_test_init_cond_records $n_cpus $short_sim $short_batch_size
fi
report_if_failed $? "INITIAL CONDITION RECORDS TEST"

# DIPOLE FIELD DRIFT PERIODS
if $dipole_drifts_test
then	
	configure SERIAL GUIDING FORWARD 29 SELF
	make_and_run output_data main_test_dipole_periods 1
fi
report_if_failed $? "DIPOLE FIELD DRIFT PERIODS"

# PITCH ANGLE DISTRIBUTION ISOTROPIZATION
if $pa_distro_iso_test
then
	configure PARALLEL GUIDING_SCATT FORWARD 29 SELF
	make_and_run output_data main_test_pa_distro_isotrop $n_cpus $long_sim $long_batch_size
fi
report_if_failed $? "PITCH ANGLE DISTRIBUTION ISOTROPIZATION"

# PITCH ANGLE SCATTERING
if $pa_scatt_test
then
	configure PARALLEL GUIDING_SCATT FORWARD 29 SELF
	make_and_run output_data main_test_pa_scatt $n_cpus $short_sim $short_batch_size
fi
report_if_failed $? "PITCH ANGLE SCATTERING"

# PERPENDICULAR DIFFUSION
if $perp_diff_test
then
	configure PARALLEL GUIDING_DIFF FORWARD 29 SELF
	make_and_run output_data main_test_perp_diff $n_cpus $long_sim $long_batch_size
fi
report_if_failed $? "PERPENDICULAR DIFFUSION"

# FULL (PERP+PARA) DIFFUSION
if $full_diff_test
then	
	configure PARALLEL PARKER FORWARD 29 SELF
	make_and_run output_data main_test_full_diff $n_cpus $long_sim $long_batch_size
fi
report_if_failed $? "FULL (PERP+PARA) DIFFUSION"

# DIFFUSIVE SHOCK ACCELERATION
if $diff_shock_acc_test
then	
	configure PARALLEL PARKER_SOURCE BACKWARD 0 SELF
	make_and_run output_data main_test_diffusive_shock_acceleration $n_cpus $long_sim $long_batch_size
	report_if_failed $? "DIFFUSIVE SHOCK ACCELERATION"
	make_and_run output_data main_postprocess_diffusive_shock_acceleration 1
	report_if_failed $? "POST-PROCESSING DIFFUSIVE SHOCK ACCELERATION"
fi

# MODULATION WITH CARTESIAN PARKER SPIRAL
if $modulation_cartesian_parker_spiral
then
	configure PARALLEL PARKER BACKWARD 25 CARTESIAN 1 0
	make_and_run cartesian_backgrounds main_generate_cartesian_solarwind_background 1
	report_if_failed $? "CARTESIAN PARKER SPIRAL BACKGROUND GENERATION"
	make_and_run output_data main_test_modulation_cartesian_parker $n_cpus $long_sim $long_batch_size
	report_if_failed $? "MODULATION CARTESIAN PARKER SPIRAL"
	make_and_run output_data main_postprocess_modulation_cartesian_parker 1
	report_if_failed $? "POST-PROCESSING MODULATION CARTESIAN PARKER SPIRAL"
fi
