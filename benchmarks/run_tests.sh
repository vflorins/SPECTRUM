#!/bin/bash

# File to hold benchmark results
results_file="benchmark_results.txt"

# Parallelization parameters
long_sim=100000
short_sim=10000
long_batch_size=100
short_batch_size=10

# Number of CPUs
n_cpus_per_socket=$(grep 'core id' /proc/cpuinfo | sort -u | wc -l)
n_sockets=$(grep 'physical id' /proc/cpuinfo | sort -u | wc -l)
n_cpus_max=$(($n_cpus_per_socket * $n_sockets))
n_cpus=${2:-$n_cpus_max}
if [[ $n_cpus -gt $n_cpus_max ]]
then
	n_cpus=$n_cpus_max
fi

# If no input is provided, explain usage
if [ "$#" -eq 0 ] || [ "${1}" = "-h" ]
then
   echo "Script to run SPECTRUM benchmark tests."
   echo "Usage: ./run_tests.sh [t] [n]"
   echo ""
   echo "   Value of t indicates which test(s) to run:"
   echo "--------------------------------------------------------------"
   echo "      1: DIPOLE FIELD VISUALIZATION (1 CPU)"
   echo "      2: TURBULENCE VIA SUPERPOSITION OF WAVES (1 CPU)"
   echo "      3: DIPOLE FIELD DRIFT PERIODS (1 CPU)"
   echo "      4: INITIAL CONDITION RECORDS TEST (2+ CPU)"
   echo "      5: PITCH ANGLE DISTRIBUTION ISOTROPIZATION (2+ CPUs)"
   echo "      6: PITCH ANGLE SCATTERING (2+ CPUs)"
   echo "      7: PERPENDICULAR DIFFUSION (2+ CPUs)"
   echo "      8: FULL (PERP+PARA) DIFFUSION (2+ CPUs)"
   echo "      9: DIFFUSIVE SHOCK ACCELERATION (2+ CPUs)"
   echo "     10: PARKER SPIRAL MAGNETIC FIELD (3 CPUs)"
   echo "     11: PARKER MODULATION WITH CARTESIAN BACKGROUND (3+ CPUs)"
   echo "    all: ALL BENCHMARK TESTS"
   echo "--------------------------------------------------------------"
   echo ""
   echo "   Value of n indicates how many processors to use for parallel tests."
   echo ""
fi

# If input `all` is provided, run all tests
if [ "${1}" = "all" ]
then
	run_all_tests=true
else
	run_all_tests=false
fi

# Function to print header
# $1: header title
# $2: destination file
function print_header {
	echo -n "${1} --- " > "${2}"
	date >> "${2}"
	echo "------------------------------------------------------------" >> "${2}"
}

# Function to return exit code
# $1: exit code of previous command
# $2: stage to report
# $3: name of log file
# $4: flag to output to stdout as well
# $5: name of test
# $6: duration of code execution in seconds
function report_if_failed {
	echo "" >> "${3}"
	echo $2 >> "${3}"
	echo "========================================" >> "${3}"
	if [ $1 -ne 0 ]
	then
		echo "***failed*** with error code ${1}" >> "${3}"
	else
		echo "completed" >> "${3}"
	fi
	echo "========================================" >> "${3}"
	echo "" >> "${3}"
	if $4
	then
		echo ""
		echo $5
		echo "========================================"
		if [ $1 -ne 0 ]
		then
			echo "did not complete... check log files"
		else
			echo "completed... check output plot"
			echo "Run-time = ${duration} s."
		fi
		echo "========================================"
		echo ""
	fi
}

# Function to go up one directory and configure test
# $1: log file
# $2-8: configure parameters
function configure_test {
	cd ..
   configure_args="--with-execution=$2 --with-trajectory=$3 --with-time_flow=$4 --with-rkmethod=$5 --with-server=$6 --enable-gsl"
   if [ "$2" = "PARALLEL" ]; then
      configure_args="--with-mpi=openmpi $configure_args"
   fi
	if [ "$6" != "SELF" ]; then
      configure_args="$configure_args --with-server_interp_order=$7 --with-server_num_gcs=$8"
   fi
	./configure $configure_args 1>> "benchmarks/logs/${1}" 2>> "benchmarks/logs/${1}"
	report_if_failed $? "CONFIGURATION" "benchmarks/logs/${1}" false
	cd - 1>> "benchmarks/logs/${1}"
}

# Function to make test
# $1: test to make
# $2: log file
function make_test {
	make $1 1>> "logs/${2}" 2>> "logs/${2}"
	report_if_failed $? "COMPILATION" "logs/${2}" false
}

# Function to run test
# $1: folder in which to run test
# $2: test to run
# $3: log file
# $4: number of processors to use
# $5: number of trajectories
# $6: trajectory batch size
function run_test {
	mkdir -p $1
	cd $1
	start_time=$SECONDS
	if [ $4 -eq 1 ]
	then
		../$2 1>> "../${results_file}" 2>> "../logs/${3}"
	else
		mpirun -np $4 ../$2 $5 $6 1>> "../${results_file}" 2>> "../logs/${3}"
	fi
	duration=$((SECONDS - start_time))
	report_if_failed $? "EXECUTION" "../logs/${3}" true $2 $duration
	cd - 1>> "../logs/${3}"
}

# Function to plot test results
# $1: plotting script
# $2: path to data
# $3: path to plots
# $4: log file
function plot_test {
	mkdir -p $3
	python "plotting_scripts/${1}" $2 $3 1>> "logs/${4}" 2>> "logs/${4}"
	report_if_failed $? "PLOTTING" "logs/${4}" false
}

# Before running this script, execute the following in the root directory:
# 1) module load mpi
# 2) autoreconf
# 3) automake --add-missing

# Print header to the results file
print_header "BENCHMARK RESULTS" $results_file

# Make log directory
mkdir -p logs

# DIPOLE FIELD VISUALIZATION
if [ "${1}" = "1" ] || $run_all_tests
then
	test_title="DIPOLE FIELD VISUALIZATION"
	log_file="log_test_dipole_visual.txt"
	echo "Currently running ${test_title} benchmark on 1 CPU"
	print_header "${test_title}" "logs/${log_file}"
	configure_test $log_file SERIAL FIELDLINE FORWARD 0 SELF
	make_test main_test_dipole_visual $log_file
	run_test output_data main_test_dipole_visual $log_file 1
	plot_test plot_test_dipole_visual.py output_data/ output_plots/ $log_file
fi

# TURBULENCE VIA SUPERPOSITION OF WAVES
if [ "${1}" = "2" ] || $run_all_tests
then
	test_title="TURBULENCE VIA SUPERPOSITION OF WAVES"
	log_file="log_test_turb_waves.txt"
	echo "Currently running ${test_title} benchmark on 1 CPU"
	print_header "${test_title}" "logs/${log_file}"
	configure_test $log_file SERIAL FIELDLINE FORWARD 29 SELF
	make_test main_test_turb_waves $log_file
	run_test output_data main_test_turb_waves $log_file 1
	plot_test plot_test_turb_waves.py output_data/ output_plots/ $log_file
fi

# DIPOLE FIELD DRIFT PERIODS
if [ "${1}" = "3" ] || $run_all_tests
then
	test_title="DIPOLE FIELD DRIFT PERIODS"
	log_file="log_test_dipole_periods.txt"
	echo "Currently running ${test_title} benchmark on 1 CPU"
	print_header "${test_title}" "logs/${log_file}"
	configure_test $log_file SERIAL GUIDING FORWARD 29 SELF
	make_test main_test_dipole_periods $log_file
	run_test output_data main_test_dipole_periods $log_file 1
	plot_test plot_test_dipole_periods.py output_data/ output_plots/ $log_file
fi

# INITIAL CONDITION RECORDS
if [ "${1}" = "4" ] || $run_all_tests
then
	test_title="INITIAL CONDITION RECORDS"
	n_cpus_req=2
	if [ "${n_cpus}" -ge "${n_cpus_req}" ]
	then
		log_file="log_test_init_cond_records.txt"
		echo "Currently running ${test_title} benchmark on ${n_cpus} CPUs"
		print_header "${test_title}" "logs/${log_file}"
		configure_test $log_file PARALLEL FOCUSED FORWARD 29 SELF
		make_test main_test_init_cond_records $log_file
		run_test output_data main_test_init_cond_records $log_file $n_cpus $short_sim $short_batch_size
		plot_test plot_test_init_cond_records.py output_data/ output_plots/ $log_file
	else
		echo "${test_title} benchmark cannot run with less than ${n_cpus_req} CPUs."
		echo ""
	fi
fi

# PITCH ANGLE DISTRIBUTION ISOTROPIZATION
if [ "${1}" = "5" ] || $run_all_tests
then
	test_title="PITCH ANGLE DISTRIBUTION ISOTROPIZATION"
	n_cpus_req=2
	if [ "${n_cpus}" -ge "${n_cpus_req}" ]
	then
		log_file="log_test_pa_distro_isotrop.txt"
		echo "Currently running ${test_title} benchmark on ${n_cpus} CPUs"
		print_header "${test_title}" "logs/${log_file}"
		configure_test $log_file PARALLEL GUIDING_SCATT FORWARD 29 SELF
		make_test main_test_pa_distro_isotrop $log_file
		run_test output_data main_test_pa_distro_isotrop $log_file $n_cpus $long_sim $long_batch_size
		plot_test plot_test_pa_distro_isotrop.py output_data/ output_plots/ $log_file
	else
		echo "${test_title} benchmark cannot run with less than ${n_cpus_req} CPUs."
		echo ""
	fi
fi

# PITCH ANGLE SCATTERING
if [ "${1}" = "6" ] || $run_all_tests
then
	test_title="PITCH ANGLE SCATTERING"
	n_cpus_req=2
	if [ "${n_cpus}" -ge "${n_cpus_req}" ]
	then
		log_file="log_test_pa_scatt.txt"
		echo "Currently running ${test_title} benchmark on ${n_cpus} CPUs"
		print_header "${test_title}" "logs/${log_file}"
		configure_test $log_file PARALLEL GUIDING_SCATT FORWARD 29 SELF
		make_test main_test_pa_scatt $log_file
		run_test output_data main_test_pa_scatt $log_file $n_cpus $short_sim $short_batch_size
		plot_test plot_test_pa_scatt.py output_data/ output_plots/ $log_file
	else
		echo "${test_title} benchmark cannot run with less than ${n_cpus_req} CPUs."
		echo ""
	fi
fi

# PERPENDICULAR DIFFUSION
if [ "${1}" = "7" ] || $run_all_tests
then
	test_title="PERPENDICULAR DIFFUSION"
	n_cpus_req=2
	if [ "${n_cpus}" -ge "${n_cpus_req}" ]
	then
		log_file="log_test_perp_diff.txt"
		echo "Currently running ${test_title} benchmark on ${n_cpus} CPUs"
		print_header "${test_title}" "logs/${log_file}"
		configure_test $log_file PARALLEL GUIDING_DIFF FORWARD 29 SELF
		make_test main_test_perp_diff $log_file
		run_test output_data main_test_perp_diff $log_file $n_cpus $long_sim $long_batch_size
		plot_test plot_test_perp_diff.py output_data/ output_plots/ $log_file
	else
		echo "${test_title} benchmark cannot run with less than ${n_cpus_req} CPUs."
		echo ""
	fi
fi

# FULL (PERP+PARA) DIFFUSION
if [ "${1}" = "8" ] || $run_all_tests
then
	test_title="FULL DIFFUSION"
	n_cpus_req=2
	if [ "${n_cpus}" -ge "${n_cpus_req}" ]
	then
		log_file="log_test_full_diff.txt"
		echo "Currently running ${test_title} benchmark on ${n_cpus} CPUs"
		print_header "${test_title}" "logs/${log_file}"
		configure_test $log_file PARALLEL PARKER FORWARD 29 SELF
		make_test main_test_full_diff $log_file
		run_test output_data main_test_full_diff $log_file $n_cpus $long_sim $long_batch_size
		plot_test plot_test_full_diff.py output_data/ output_plots/ $log_file
	else
		echo "${test_title} benchmark cannot run with less than ${n_cpus_req} CPUs."
		echo ""
	fi
fi

# DIFFUSIVE SHOCK ACCELERATION
if [ "${1}" = "9" ] || $run_all_tests
then
	test_title="DIFFUSIVE SHOCK ACCELERATION"
	n_cpus_req=2
	if [ "${n_cpus}" -ge "${n_cpus_req}" ]
	then
		log_file="log_test_diff_shock_acc.txt"
		echo "Currently running ${test_title} benchmark on ${n_cpus} CPUs"
		print_header "${test_title}" "logs/${log_file}"
		configure_test $log_file PARALLEL PARKER_SOURCE BACKWARD 0 SELF
		make_test main_test_diff_shock_acc $log_file
		run_test output_data main_test_diff_shock_acc $log_file $n_cpus $long_sim $long_batch_size
		
		test_title="DIFFUSIVE SHOCK ACCELERATION POST-PROCESS"
		log_file="log_postprocess_diff_shock_acc.txt"
		echo "Currently running ${test_title} benchmark on 1 CPU"
		print_header "${test_title}" "logs/${log_file}"
		make_test main_postprocess_diff_shock_acc $log_file
		run_test output_data main_postprocess_diff_shock_acc $log_file 1
		plot_test plot_test_diff_shock_acc.py output_data/ output_plots/ $log_file
	else
		echo "${test_title} benchmark cannot run with less than ${n_cpus_req} CPUs."
		echo ""
	fi
fi

# PARKER SPIRAL MAGNETIC FIELD
if [ "${1}" = "10" ] || $run_all_tests
then
	test_title="PARKER SPIRAL MAGNETIC FIELD"
	n_cpus_req=3
	if [ "${n_cpus}" -ge "${n_cpus_req}" ]
	then
		test_title="CARTESIAN PARKER SPIRAL GENERATION"
		log_file="log_generate_cartesian_parker_spiral.txt"
		echo "Currently running ${test_title} benchmark on 1 CPU"
		print_header "${test_title}" "logs/${log_file}"
		configure_test $log_file PARALLEL FIELDLINE FORWARD 25 CARTESIAN 1 0
		make_test main_generate_cartesian_parker_spiral $log_file
		run_test cartesian_backgrounds main_generate_cartesian_parker_spiral $log_file 1
		
		test_title="CARTESIAN PARKER SPIRAL"
		log_file="log_test_cartesian_parker_spiral.txt"
		echo "Currently running ${test_title} benchmark on 3 CPUs"
		print_header "${test_title}" "logs/${log_file}"
		make_test main_test_cartesian_parker_spiral $log_file
		run_test output_data main_test_cartesian_parker_spiral $log_file 3 1 1

	   test_title="SOLARWIND PARKER SPIRAL"
		log_file="log_test_solarwind_parker_spiral.txt"
		echo "Currently running ${test_title} benchmark on 1 CPU"
		print_header "${test_title}" "logs/${log_file}"
		configure_test $log_file SERIAL FIELDLINE FORWARD 25 SELF
		make_test main_test_solarwind_parker_spiral $log_file
		run_test output_data main_test_solarwind_parker_spiral $log_file 1
		plot_test plot_test_parker_spiral.py output_data/ output_plots/ $log_file
	else
		echo "${test_title} benchmark cannot run with less than ${n_cpus_req} CPUs."
		echo ""
	fi
fi

# PARKER MODULATION WITH CARTESIAN BACKGROUND
if [ "${1}" = "11" ] || $run_all_tests
then
	test_title="PARKER MODULATION WITH CARTESIAN BACKGROUND"
	n_cpus_req=3
	if [ "${n_cpus}" -ge "${n_cpus_req}" ]
	then
		test_title="CARTESIAN SOLAR WIND BACKGROUND GENERATION"
		log_file="log_generate_cartesian_solarwind_background.txt"
		echo "Currently running ${test_title} benchmark on 1 CPU"
		print_header "${test_title}" "logs/${log_file}"
		configure_test $log_file PARALLEL PARKER BACKWARD 0 CARTESIAN 1 0
		make_test main_generate_cartesian_solarwind_background $log_file
		run_test cartesian_backgrounds main_generate_cartesian_solarwind_background $log_file 1
		
		test_title="MODULATION CARTESIAN SOLAR WIND"
		log_file="log_test_modulation_cartesian_parker.txt"
		echo "Currently running ${test_title} benchmark on ${n_cpus} CPUs"
		print_header "${test_title}" "logs/${log_file}"
		make_test main_test_modulation_cartesian_parker $log_file
		run_test output_data main_test_modulation_cartesian_parker $log_file $n_cpus $long_sim $long_batch_size
		
		test_title="MODULATION CARTESIAN SOLAR WIND POST-PROCESS"
		log_file="log_postprocess_modulation_cartesian_parker.txt"
		echo "Currently running ${test_title} benchmark on 1 CPU"
		print_header "${test_title}" "logs/${log_file}"
		make_test main_postprocess_modulation_cartesian_parker $log_file
		run_test output_data main_postprocess_modulation_cartesian_parker $log_file 1
		plot_test plot_test_modulation_cartesian_parker.py output_data/ output_plots/ $log_file
	else
		echo "${test_title} benchmark cannot run with less than ${n_cpus_req} CPUs."
		echo ""
	fi
fi
