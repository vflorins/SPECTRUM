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
dipole_drifts_test=true
pa_distro_iso_test=true
pa_scatt_test=true
perp_diff_test=true
full_diff_test=true
turb_waves_test=true
parker_sprial_test=true
init_cond_record_test=true
modulation_cartesian_parker_spiral=false

# Function to go up one directory and configure code
function configure {
	cd ..
	./configure CXXFLAGS="-Ofast" --with-mpi=openmpi --with-trajectory=$1 --with-rkmethod=$2 --with-server=$3
	cd -
}

# Function to make and run test
function make_and_run {
	make $1
	if [ $2 -eq 1 ]
	then
		./$1 >> $results_file
	else
		mpirun -np $2 $1 $3 $4 >> $results_file
	fi
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
mkdir -p "output_data"
echo -n "Benchmark results " > $results_file
date >> $results_file
echo "------------------------------------------------------------" >> $results_file

# DIPOLE FIELD DRIFT PERIODS
if $dipole_drifts_test
then	
	configure GUIDING 29 SELF
	# configure FOCUSED 29 SELF
	# configure LORENTZ 29 SELF
	make_and_run main_test_dipole_periods 1
fi
report_if_failed $? "DIPOLE FIELD DRIFT PERIODS"

# PITCH ANGLE DISTRIBUTION ISOTROPIZATION
if $pa_distro_iso_test
then
	configure GUIDING_SCATT 29 SELF
	# configure GUIDING_DIFF_SCATT 29 SELF
	# configure FOCUSED_SCATT 29 SELF
	# configure FOCUSED_DIFF_SCATT 29 SELF
	make_and_run main_test_pa_distro_isotrop $n_cpus $long_sim $long_batch_size
fi
report_if_failed $? "PITCH ANGLE DISTRIBUTION ISOTROPIZATION"

# PITCH ANGLE SCATTERING
if $pa_scatt_test
then
	configure GUIDING_SCATT 29 SELF
	# configure GUIDING_DIFF_SCATT 29 SELF
	# configure FOCUSED_SCATT 29 SELF
	# configure FOCUSED_DIFF_SCATT 29 SELF
	make_and_run main_test_pa_scatt $n_cpus $short_sim $short_batch_size
fi
report_if_failed $? "PITCH ANGLE SCATTERING"

# PERPENDICULAR DIFFUSION
if $perp_diff_test
then
	configure GUIDING_DIFF 29 SELF
	# configure GUIDING_DIFF_SCATT 29 SELF
	# configure FOCUSED_DIFF 29 SELF
	# configure FOCUSED_DIFF_SCATT 29 SELF
	make_and_run main_test_perp_diff $n_cpus $long_sim $long_batch_size
fi
report_if_failed $? "PERPENDICULAR DIFFUSION"

# FULL (PERP+PARA) DIFFUSION
if $full_diff_test
then	
	configure PARKER 29 SELF
	make_and_run main_test_full_diff $n_cpus $long_sim $long_batch_size
fi
report_if_failed $? "FULL (PERP+PARA) DIFFUSION"

# TURBULENCE VIA SUPERPOSITION OF WAVES
if $turb_waves_test
then
	configure FIELDLINE 29 SELF
	make_and_run main_test_turb_waves 1
fi
report_if_failed $? "TURBULENCE VIA SUPERPOSITION OF WAVES"

# PARKER SPIRAL SOLAR WIND
if $parker_sprial_test
then
	configure FIELDLINE 29 SELF
	# configure LORENTZ 29 SELF
	# configure FOCUSED 29 SELF
	# configure GUIDING 29 SELF
	make_and_run main_test_parker_spiral 1
fi
report_if_failed $? "PARKER SPIRAL SOLAR WIND"

# INITIAL CONDITION RECORDS TEST
if $init_cond_record_test
then	
	# configure GUIDING 29 SELF
	configure FOCUSED 29 SELF
	# configure LORENTZ 29 SELF
	make_and_run main_test_init_cond_records $n_cpus $short_sim $short_batch_size
fi
report_if_failed $? "INITIAL CONDITION RECORDS TEST"

# MODULATION WITH CARTESIAN PARKER SPIRAL
if $modulation_cartesian_parker_spiral
then
	configure PARKER 29 CARTESIAN
	# configure PARKER 29 SELF
	make_and_run main_test_modulation_cartesian_parker $n_cpus $long_sim $long_batch_size
	report_if_failed $? "MODULATION CARTESIAN PARKER SPIRAL"
	make_and_run main_postprocess_modulation_cartesian_parker 1
	report_if_failed $? "POST-PROCESSING MODULATION CARTESIAN PARKER SPIRAL"
fi