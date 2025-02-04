// mpicxx -I.. -Wall -Wno-comment -std=c++20 -DGEO_DEBUG -g -O0 -o main_test_buffered_block main_test_buffered_block.cc ../geometry/distance_map.cc ../geodesic/spherical_tesselation.cc ../geodesic/requestable_tesselation.cc ../geodesic/traversable_tesselation.cc ../geodesic/geodesic_sector.cc ../geodesic/grid_block.cc ../geodesic/stenciled_block.cc ../geodesic/grid_block.cc ../common/params.cc ../common/data_container.cc ../common/vectors.cc -lgsl -lgslcblas

#include "geodesic/buffered_block.hh"
#include "geodesic/traversable_tesselation.hh"

using namespace Spectrum;

//#define POLY_TYPE POLY_HEXAHEDRON
#define POLY_TYPE POLY_ICOSAHEDRON
#define MAX_FACE_DIVISION 5

#if (POLY_TYPE == POLY_HEXAHEDRON)
#define verts_per_face 4
#else
#define verts_per_face 3
#endif

const double rmin = 1.0;
const double rmax = 2.0;

const int sector_division = 0;
const int face_division = 3;
const int ghost_width = 2;
const int sector_width = Pow2(face_division - sector_division);

const int shells_per_slab = 8;
const int n_slabs = 3;
const int ghost_height = 2;

int main(int argc, char *argv[])
{
   int edges_per_vert = EdgesAtVert<verts_per_face>();
   int n_sect_parts_actual, sectors[edges_per_vert];
   bool corners[verts_per_face];
   bool borders[2];
   int n_sectors, n_blocks, existing_sectors[2];

   MPI_Config mpi_config(argc, argv);
   int my_rank = MPI_Config::glob_comm_rank;
   bool force_both_blocks_rank0 = false;

   PrintMessage(__FILE__, __LINE__, "Size of glob commm is " + std::to_string(MPI_Config::glob_comm_size), MPI_Config::is_master);

// Can only test with 1 or 2 processes
   if(MPI_Config::work_comm_size > 2) return 1;

// Create a Tesselation object
   TraversableTesselation<POLY_TYPE, MAX_FACE_DIVISION> tesselation;
   n_sectors = tesselation.NFaces(sector_division);
   n_blocks = n_slabs * n_sectors;

// Set up a radial map
   DataContainer container;
   container.Clear();
   container.Insert(rmin);
   container.Insert(rmax);
   std::shared_ptr<DistanceBase> dist_map = std::make_shared<DistanceExponential>();
   dist_map->SetupObject(container);

// Currently both blocks are in the same slab - change this code to test TFACE, TEDGE, and VERTX exchanges
   int existing_slab = 1;
   borders[0] = (existing_slab == 0 ? true : false);
   borders[1] = (existing_slab == n_slabs - 1 ? true : false);

// Pick one vertex and two of its face neighbors - this allows testing RFACE or REDGE exchanges
   int base_vert = 5;
   n_sect_parts_actual = tesselation.FaceNeighborsOfVert(sector_division, base_vert, sectors);
   existing_sectors[0] = sectors[2];
   existing_sectors[1] = sectors[3];

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Create an array of blocks, as in DomainPartition::DomainPartition()
//----------------------------------------------------------------------------------------------------------------------------------------------------

   std::vector<BufferedBlock<verts_per_face, int>> blocks_local;
   int* face_index;
   int* vert_index;
   GeoVector* vcart;

// Always allocate on rank 0, conditionally on rank 1
   if ((my_rank == 0) || !force_both_blocks_rank0) {
      blocks_local.emplace_back();
      blocks_local.back().SetDimensions(sector_width, ghost_width, shells_per_slab, ghost_height, false);
      blocks_local.back().SetIndex(existing_slab * n_sectors + existing_sectors[my_rank]);

      face_index = new int[blocks_local.back().TotalFaces()];
      vert_index = new int[blocks_local.back().TotalVerts()];
      vcart = new GeoVector[blocks_local.back().TotalVerts()];

      tesselation.GetAllInsideFaceNative(sector_division, existing_sectors[my_rank], face_division, ghost_width, face_index, vert_index, corners);
      tesselation.FillVertCoordArrays(blocks_local.back().TotalVerts(), vert_index, vcart);
      blocks_local.back().AssociateMesh(existing_slab / (double)n_slabs, (existing_slab + 1.0) / (double)n_slabs, corners, borders, vcart, dist_map, false);

      delete[] face_index;
      delete[] vert_index;
      delete[] vcart;
      std::cerr << "Block number " << blocks_local.back().GetIndex() << " created\n";
   }

// Serial run or both blocks forced on one process
   if ((MPI_Config::work_comm_size == 1) || force_both_blocks_rank0) {
      blocks_local.emplace_back();
      blocks_local.back().SetDimensions(sector_width, ghost_width, shells_per_slab, ghost_height, false);
      blocks_local.back().SetIndex(existing_slab * n_sectors + existing_sectors[1]);

      face_index = new int[blocks_local.back().TotalFaces()];
      vert_index = new int[blocks_local.back().TotalVerts()];
      vcart = new GeoVector[blocks_local.back().TotalVerts()];

      tesselation.GetAllInsideFaceNative(sector_division, existing_sectors[1], face_division, ghost_width, face_index, vert_index, corners);
      tesselation.FillVertCoordArrays(blocks_local.back().TotalVerts(), vert_index, vcart);
      blocks_local.back().AssociateMesh(existing_slab / (double)n_slabs, (existing_slab + 1.0) / (double)n_slabs, corners, borders, vcart, dist_map, false);

      delete[] face_index;
      delete[] vert_index;
      delete[] vcart;
      std::cerr << "Block number " << blocks_local.back().GetIndex() << " created\n";
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Create _all_ exchange sites of type REDGE as in DomainPartition::SetUpExchangeSites()
//----------------------------------------------------------------------------------------------------------------------------------------------------

   int n_sites_sect, buf_size, sidx, part;
   int ranks[2 * edges_per_vert], labels[2 * edges_per_vert], verts[verts_per_face];
   std::vector<std::shared_ptr<ExchangeSite<int>>> exch_sites;
   std::vector<int> exch_sites_local;
   
// A process might have no blocks, but the site still needs to be registered with correct buffer size
   BufferedBlock<verts_per_face, int> block_temp;
   block_temp.SetDimensions(sector_width, ghost_width, shells_per_slab, ghost_height, false);
   buf_size = block_temp.GetBufferSize(GEONBR_REDGE);
   PrintMessage(__FILE__, __LINE__, "Buffer size is " + std::to_string(buf_size), MPI_Config::is_master);

   sidx = 0;
   for (auto slab = 0; slab < n_slabs; slab++) {
      for (auto vert = 0; vert < tesselation.NVerts(sector_division); vert++) {
         exch_sites.emplace_back(std::make_shared<ExchangeSite<int>>());
         n_sect_parts_actual = tesselation.FaceNeighborsOfVert(sector_division, vert, sectors);

// Loop over participating sectors - only two exist, the rest are assigned MPI_PROC_NULL.
         part = 0;
         for (auto it = 0; it < n_sect_parts_actual; it++) {
            labels[part] = slab * tesselation.NFaces(sector_division) + sectors[it];

// The first block is always on process 0, while the second could be on either
            if (sectors[it] == existing_sectors[0]) ranks[part] = 0;
            else if (sectors[it] == existing_sectors[1]) {
               if ((MPI_Config::work_comm_size == 1) || force_both_blocks_rank0) ranks[part] = 0;
               else ranks[part] = 1;
            }

// All the other parts are missing, but a valid rank is required, so we use rank 0
            else ranks[part] = 0;
            part++;
         };
            
// Initialize the exchange site and add it to the list if it is local to this process
         exch_sites.back()->SetUpProperties(sidx, n_sect_parts_actual, buf_size, ranks, labels);
         if (exch_sites.back()->GetCommSize() > 0) exch_sites_local.push_back(sidx);
         exch_sites.back()->FillBuffers(-2);
         sidx++;
      };
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Inform the blocks about their exchange sites as in DomainPartition::ExportExchangeSites()
//----------------------------------------------------------------------------------------------------------------------------------------------------

   std::vector<std::shared_ptr<ExchangeSite<int>>> exch_sites_block;

   if(MPI_Config::is_master) std::cerr << std::endl;

   for (auto block = 0; block < blocks_local.size(); block++) {
      exch_sites_block.clear();
      n_sites_sect = tesselation.VertNeighborsOfFace(sector_division, blocks_local[block].GetIndex() % n_sectors, verts);

      for (auto iv = 0; iv < n_sites_sect; iv++) {
         sidx = existing_slab * tesselation.NVerts(sector_division) + verts[iv];
         exch_sites_block.push_back(exch_sites[sidx]);

         if(MPI_Config::is_master) {
            std::cerr << "Adding exchange site " << sidx << " to block " << blocks_local[block].GetIndex() << std::endl;
         };

      };
      blocks_local[block].ImportExchangeSites(GEONBR_REDGE, exch_sites_block);
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Perform the exchange as in DomainPartition::ExchangeAll(); most sites are empty and should not do anything at all
//----------------------------------------------------------------------------------------------------------------------------------------------------

   std::cerr << std::endl;
   for (auto site = 0; site < blocks_local[0].GetExchangeSiteCount(GEONBR_REDGE); site++) {
      for (part = 0; part < blocks_local[0].GetPartCount(GEONBR_REDGE); part++) {
         std::cerr << "Site " << site << " Part " << part << " face map:\n";
         blocks_local[0].PrintBufferMap(GEONBR_REDGE, site, part);
      };
   };

   blocks_local[0].FillWithIndexData();
   blocks_local[1].FillWithIndexData();
   blocks_local[0].PrintContents();

   for (auto block = 0; block < blocks_local.size(); block++) blocks_local[block].PackBuffers(GEONBR_REDGE);
   for (auto sidx : exch_sites_local) {
      exch_sites[sidx]->Exchange();
   };
   for (auto block = 0; block < blocks_local.size(); block++) blocks_local[block].UnPackBuffers(GEONBR_REDGE);

   blocks_local[0].PrintContents();

};

