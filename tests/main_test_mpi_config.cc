// mpicxx -I.. -Wall -DGEO_DEBUG -O0 -g -std=c++20 -o main_test_mpi_config main_test_mpi_config.cc ../common/mpi_config.cc

#include <iostream>
#include <iomanip>

#include "common/mpi_config.hh"

using namespace Spectrum;

// Create an instance of the class to initialize all MPI objects
void InitializeMPI(int argc, char** argv)
{
   static MPI_Config mpi_config(argc, argv);
};

int main(int argc, char** argv)
{
   InitializeMPI(argc, argv);
   TestMPIConfig();
};

