/*!
\file domain_partition.cc
\brief Implements the class responsible for partitioning of the simuilation into blocks
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include <pthread.h>
#ifdef USE_SILO
#include <silo.h>
//#include <pmpio.h>
#endif

#include "domain_partition.hh"

namespace Spectrum {

/*!
\author Vladimir Florinski
\date 06/20/2024
\param[in[ n_shells Total number of shells
\param[in[ dist_map Distance scaling function
\param[in[ face_div Face division
*/
DomainPartition::DomainPartition(int n_shells, const std::shared_ptr<DistanceBase> dist_map, double face_div, const std::shared_ptr<MPI_Config> mpi_config_in)
{
   int shells_per_slab, sector_width, bidx;
   double slab_height_fp;

   n_shells_global = n_shells;
   face_division = face_div;
   distance_map = dist_map;
   mpi_config = mpi_config_in;

   n_faces_global = tesselation.NFaces(face_division);

// Perform domain decomposition based on sector size
   slab_height_fp = n_shells_global;
   sector_width = Pow2(face_division - sector_division);
   n_sectors = tesselation.NFaces(sector_division);
   sector_division = 0;
   n_slabs = 1;
   
   while((n_sectors * n_slabs < mpi_config->n_workers) && (slab_height_fp >= min_block_height) && (sector_width >= min_block_width)) {
      if(sector_width < slab_height_fp) {
         n_slabs++;
         slab_height_fp = n_shells_global / n_slabs;
      }
      else {
         sector_division++;
         n_sectors = tesselation.NFaces(sector_division);
         sector_width = Pow2(face_division - sector_division);
      };
   };
   thinner_slab_height = n_shells_global / n_slabs;
   n_thinner_slabs = n_slabs * (thinner_slab_height + 1) - n_shells_global;
   n_faces_global = tesselation.NFaces(face_division);
   n_blocks = n_slabs * n_sectors;

// Allocate space for slab interfaces and set up the radial limits
   slab_interfaces = new double[n_slabs + 1];
   slab_interfaces[0] = distance_map->GetPhysical(0.0);
   slab_interfaces[n_slabs] = distance_map->GetPhysical(1.0);

// Assign blocks in a round-robin fashion. The vector creates copies of an empty template block; the actual memory allocation happens within the vector.
   bidx = 0;
   BufferedBlock<VERTS_PER_FACE> block_templ;
   for(auto slab = 0; slab < n_slabs; slab++) {
      shells_per_slab = (slab < n_thinner_slabs ? thinner_slab_height : thinner_slab_height + 1);
      if(slab) slab_interfaces[slab] = slab_interfaces[slab - 1] + shells_per_slab / n_shells_global;

      for(auto sector = 0; sector < n_sectors; sector++) {
         if(MakePeriodic(bidx, mpi_config->work_comm_size) == mpi_config->work_comm_rank) {
            blocks_local.push_back(block_templ);
            blocks_local.back().SetDimensions(sector_width, GHOST_WIDTH, shells_per_slab, GHOST_HEIGHT, false);
            blocks_local.back().SetIndex(bidx);
         };
         bidx++;
      };
   };

// Print the mesh summary
   if(mpi_config->is_master) {
      std::cerr << std::endl;
      std::cerr << "-------------------------------------------------------------------------------\n";
      std::cerr << "Mesh summary\n";
      std::cerr << "-------------------------------------------------------------------------------\n";
      std::cerr << "T-face division               |" << std::setw(15) << face_division << std::endl;
      std::cerr << "Sector division               |" << std::setw(15) << sector_division << std::endl;
      std::cerr << "Number of slabs               |" << std::setw(15) << n_slabs << std::endl;
      std::cerr << "Number of sectors             |" << std::setw(15) << n_sectors << std::endl;
      std::cerr << "Blocks in the mesh            |" << std::setw(15) << n_slabs * n_sectors << std::endl;
      std::cerr << "Total shells                  |" << std::setw(15) << n_shells_global << std::endl;
      std::cerr << "Total faces                   |" << std::setw(15) << n_faces_global << std::endl;
      std::cerr << "Total zones                   |" << std::setw(15) << n_shells_global * n_faces_global << std::endl;
      std::cerr << "-------------------------------------------------------------------------------\n";
   };
};

/*!
\author Vladimir Florinski
\date 06/21/2024
*/
DomainPartition::~DomainPartition()
{
   delete[] slab_interfaces;
};

/*!
\author Vladimir Florinski
\date 06/21/2024
\param[in] pos Vector position
\return The block global index, -1 if vector is outside the boundaries
*/
int DomainPartition::LocateBlock(const GeoVector& pos) const
{
   int slab, sect;
   
// Locate the slab
   slab = LocateInArray(0, n_slabs, slab_interfaces, pos.Norm(), false);
   if(slab < 0) return -1;

// Locate the sector using the tesselation
   sect = tesselation.Locate(sector_division, UnitVec(pos));

   return GetBlockIndex(slab, sect);
};

/*!
\author Vladimir Florinski
\date 06/28/2024
*/
void DomainPartition::SetUpExchangeSites(void)
{
   int site;
   NeighborType ntype;

// Figure out the number of sites of each type
   exch_site_count[GEONBR_TFACE] = (n_slabs + 1) * tesselation.NFaces(sector_division);
   exch_site_count[GEONBR_RFACE] =  n_slabs      * tesselation.NEdges(sector_division);
   exch_site_count[GEONBR_TEDGE] = (n_slabs + 1) * tesselation.NEdges(sector_division);
   exch_site_count[GEONBR_REDGE] =  n_slabs      * tesselation.NVerts(sector_division);
   exch_site_count[GEONBR_VERTX] = (n_slabs + 1) * tesselation.NVerts(sector_division);

// TODO Create the exchange sites
   ExchangeSite<ConservedVariables> exch_site_templ;
   for(ntype = GEONBR_TFACE; ntype <= GEONBR_VERTX; GEO_INCR(ntype, NeighborType)) {
      for(site = 0; site < exch_site_count[ntype]; site++) {
         exch_sites[ntype].push_back(exch_site_templ);








         
//         exch_sites[ntype].back().SetUpProperties(.....);

// This is a local site, add it to the list
         if(exch_sites[ntype].back().site_comm != MPI_COMM_NULL) exch_site_idx[ntype].push_back(site);
      };
   };

// TODO Generate list of exchange sites for each block


// TODO Inform the blocks about their exchange sites

};


/*!
\author Vladimir Florinski
\date 06/27/2024
*/
void DomainPartition::ExchangeAll(void)
{
   NeighborType ntype;

// All work is done by the exchange sites
   for(ntype = GEONBR_TFACE; ntype <= GEONBR_VERTX; GEO_INCR(ntype, NeighborType)) {
      for(auto idx = 0; idx < exch_site_idx[ntype].size(); idx++) {
         exch_sites[ntype][exch_site_idx[ntype][idx]].Exchange();
      };
   };
};

};
