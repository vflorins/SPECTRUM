// mpicxx -I.. -Wall -std=c++20 -DGEO_DEBUG -DUSE_SILO -g -O0 -o main_test_exchange_site main_test_exchange_site.cc ../common/mpi_config.cc
// mpirun -n 10 --oversubscribe --host localhost:10 ./main_test_exchange_site

#include "common/exchange_site.hh"

using namespace Spectrum;

const int test_rank = 1;

int main(int argc, char *argv[])
{
   MPI_Config mpi_config(argc, argv);
   TestExchange(test_rank);
};

