// mpicxx -I.. -Wall -Wno-comment -DGEO_DEBUG -DUSE_MPI -g -O0 -o main_test_exchange_site main_test_exchange_site.cc ../common/mpi_config.cc
// mpirun -n 10 --host localhost:10 ./main_test_exchange_site

#include <iostream>
#include <iomanip>
#include "common/exchange_site.hh"

using namespace Spectrum;

const int test_rank = 2;

int main(int argc, char *argv[])
{
   std::shared_ptr<MPI_Config> mpi_config = std::make_shared<MPI_Config>(argc, argv);
   TestExchange(mpi_config, test_rank);
   return 0;
};

