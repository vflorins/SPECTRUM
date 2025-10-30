#include "geodesic/stenciled_block.hh"
#include "geodesic/traversable_tesselation.hh"
#include "geodesic/domain_partition.hh"

using namespace Spectrum;

const double rmin = 0.3;
const double rmax = 5.0;

const int face_div = 5;
const int n_shells = 256;

int main(int argc, char *argv[])
{
   std::string fname;
   MPI_Config mpi_config(argc, argv);

// Set up a radial map
   DataContainer container;
   container.Clear();
   container.Insert(rmin);
   container.Insert(rmax);
   container.Insert(2.0);
   std::shared_ptr<DistanceBase> dist_map = std::make_shared<DistancePowerLaw>();
   dist_map->SetupObject(container);

   DomainPartition<StenciledBlock<VERTS_PER_FACE>> domain_partition(n_shells, dist_map, face_div);
   domain_partition.PrintTopology();

   for (auto block = 0; block < domain_partition.GetNBlocks(); block++) {
      fname = "zone_centers_" + std::to_string(domain_partition.GetBlock(block).GetIndex()) + ".dat";
      domain_partition.GetBlock(block).PrintZoneCentroids(fname);
   };
};
