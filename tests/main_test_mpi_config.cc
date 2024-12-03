// mpicxx -I.. -Wall -DGEO_DEBUG -O0 -g -std=c++20 -o main_test_mpi_config main_test_mpi_config.cc ../common/mpi_config.cc

#include <iostream>
#include <iomanip>

#include "common/mpi_config.hh"

using namespace Spectrum;

int main(int argc, char** argv)
{
// Create an instance of the class to initialize all MPI objects
   MPI_Config mpi_config(argc, argv);

// From this point on the code that includes "mpi_config.hh" can access all data members of the "MPI_Config" type without referring to any instance
   TestMPIConfig();
};

