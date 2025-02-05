// mpicxx -I.. -Wall -std=c++20 -DGEO_DEBUG -DUSE_SILO -g -O0 -o main_test_shared_site main_test_shared_site.cc ../common/mpi_config.cc
// mpirun -n 10 --oversubscribe --host localhost:10 ./main_test_shared_site

#include "common/shared_site.hh"

using namespace Spectrum;

const int test_rank = 2;

int main(int argc, char *argv[])
{
   MPI_Config mpi_config(argc, argv);
   TestShared(test_rank);
};

