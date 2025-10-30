#include "geodesic/buffered_block.hh"
#include "geodesic/traversable_tesselation.hh"
#include "geodesic/domain_partition.hh"

using namespace Spectrum;

const double rmin = 1.0;
const double rmax = 5.0;

const int face_div = 3;
const int ghost_width = 2;
const int ghost_height = 2;
const int n_shells = 10;

int main(int argc, char *argv[])
{
   MPI_Config mpi_config(argc, argv);

// Set up a radial map
   DataContainer container;
   container.Clear();
   container.Insert(rmin);
   container.Insert(rmax);
   std::shared_ptr<DistanceBase> dist_map = std::make_shared<DistanceExponential>();
   dist_map->SetupObject(container);

// Create the simulation object
//   DomainPartition<GridBlock<VERTS_PER_FACE>> domain_partition(n_shells, dist_map, face_div);
//   DomainPartition<StenciledBlock<VERTS_PER_FACE>> domain_partition(n_shells, dist_map, face_div);
   DomainPartition<BufferedBlock<VERTS_PER_FACE, int>> domain_partition(n_shells, dist_map, face_div);

   domain_partition.PrintTopology();
   domain_partition.Save(0);
};
