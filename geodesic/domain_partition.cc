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
template <typename datatype, typename blocktype>
DomainPartition<datatype, blocktype>::DomainPartition(int n_shells, const std::shared_ptr<DistanceBase> dist_map, double face_div,
                                                      const std::shared_ptr<MPI_Config> mpi_config_in)
{
   int slab, slab_i, sector, shells_per_slab, sector_width, bidx;
   double slab_height_fp;
   bool borders[2], corners[verts_per_face];

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
   for(slab_i = 1; slab_i < n_slabs; slab_i++) {
      shells_per_slab = (slab_i <= n_thinner_slabs ? thinner_slab_height : thinner_slab_height + 1);
      slab_interfaces[slab_i] = slab_interfaces[slab_i - 1] + (double)shells_per_slab / (double)n_shells_global;
   };

   block_ranks = new int[n_blocks];
   int* face_index;
   int* vert_index;
   GeoVector* vcart;

// Assign blocks in a round-robin fashion. The vector creates copies of an empty template block; the actual memory allocation happens within the vector.
   bidx = 0;
   blocktype block_templ;
   for(slab = 0; slab < n_slabs; slab++) {
      shells_per_slab = (slab < n_thinner_slabs ? thinner_slab_height : thinner_slab_height + 1);
      borders[0] = (slab == n_slabs - 1 ? 1 : 0);
      borders[1] = (slab == 0 ? 1 : 0);

      for(sector = 0; sector < n_sectors; sector++) {
         block_ranks[bidx] = MakePeriodic(bidx, mpi_config->work_comm_size);
         if(block_ranks[bidx] == mpi_config->work_comm_rank) {
            blocks_local.push_back(block_templ);
            blocks_local.back().SetDimensions(sector_width, GHOST_WIDTH, shells_per_slab, GHOST_HEIGHT, false);
            blocks_local.back().SetIndex(bidx);

            face_index = new int[blocks_local.back().TotalFaces()];
            vert_index = new int[blocks_local.back().TotalVerts()];
            vcart = new GeoVector[blocks_local.back().TotalVerts()];
            
            tesselation.GetAllInsideFaceNative(sector_division, sector, face_division, GHOST_WIDTH, face_index, vert_index, corners);
            tesselation.FillVertCoordArrays(blocks_local.back().TotalVerts(), vert_index, vcart);
            blocks_local.back().AssociateMesh(slab_interfaces[slab], slab_interfaces[slab + 1], corners, borders, vcart, distance_map);

            delete[] face_index;
            delete[] vert_index;
            delete[] vcart;
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
template <typename datatype, typename blocktype>
DomainPartition<datatype, blocktype>::~DomainPartition()
{
   delete[] slab_interfaces;
   delete[] block_ranks;
};

/*!
\author Vladimir Florinski
\date 06/21/2024
\param[in] pos Vector position
\return The block global index, -1 if vector is outside the boundaries
*/
template <typename datatype, typename blocktype>
int DomainPartition<datatype, blocktype>::LocateBlock(const GeoVector& pos) const
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
template <typename datatype, typename blocktype>
void DomainPartition<datatype, blocktype>::SetUpExchangeSites(void)
{
   int sector_width, n_slab_parts_actual, n_sect_parts_actual, buf_size, slab, slab_i, face, edge, vert, is, it, sidx, part;
   ExchangeSite<datatype> exch_site_templ;

// Initialize two buffered blocks, one for thinner slabs and one for thicker slabs.
   sector_width = Pow2(face_division - sector_division);
   blocktype block_temp_thin, block_temp_thick;
   block_temp_thin .SetDimensions(sector_width, GHOST_WIDTH, thinner_slab_height    , GHOST_HEIGHT, false);
   block_temp_thick.SetDimensions(sector_width, GHOST_WIDTH, thinner_slab_height + 1, GHOST_HEIGHT, false);

// Make the arrays large enough for any site
   int faces[edges_per_vert], ranks[2 * edges_per_vert], labels[2 * edges_per_vert];

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TFACE
//----------------------------------------------------------------------------------------------------------------------------------------------------

// The TFACE buffer size does not depend on the slab thickness; can use either thin or thick
   buf_size = block_temp_thin.GetBufferSize(GEONBR_TFACE);
   n_sect_parts_actual = n_sect_parts[GEONBR_TFACE];

// Loop on slab interfaces
   sidx = 0;
   for(slab_i = 0; slab_i <= n_slabs; slab_i++) {

// Boundaries
      if((slab_i == 0) || (slab_i == n_slabs)) n_slab_parts_actual = n_slab_parts[GEONBR_TFACE] - 1;
      else n_slab_parts_actual = n_slab_parts[GEONBR_TFACE];

// Loop on faces
      for(face = 0; face < tesselation.NFaces(sector_division); face++) {
         exch_sites[GEONBR_TFACE].push_back(exch_site_templ);

// Loop over participating slabs
         part = 0;
         for(is = 0; is < n_slab_parts_actual; is++) {

// Loop over participating sectors
            for(it = 0; it < n_sect_parts_actual; it++) {
               if(slab_i == 0) labels[part] = face;
               else labels[part] = (slab_i + is - 1) * tesselation.NFaces(sector_division) + face;
               ranks[part] = block_ranks[labels[part]];
               part++;
            };
         };
            
// Initialize the exchange site and add it to the list if it is local to this process
         exch_sites[GEONBR_TFACE].back().SetUpProperties(sidx, n_slab_parts_actual * n_sect_parts_actual, buf_size, ranks, labels, mpi_config);
         if(exch_sites[GEONBR_TFACE].back().site_comm != MPI_COMM_NULL) exch_site_idx[GEONBR_TFACE].push_back(sidx);
         sidx++;
      };
   };         

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RFACE
//----------------------------------------------------------------------------------------------------------------------------------------------------

// This site has only one participating slab
   n_slab_parts_actual = n_slab_parts[GEONBR_RFACE];

// Loop on slabs
   sidx = 0;
   for(slab = 0; slab < n_slabs; slab++) {

// The buffer size depends on the slab thickness
      buf_size = (slab < n_thinner_slabs ? block_temp_thin.GetBufferSize(GEONBR_RFACE) : block_temp_thick.GetBufferSize(GEONBR_RFACE));

// Loop on edges
      for(edge = 0; edge < tesselation.NEdges(sector_division); edge++) {
         exch_sites[GEONBR_RFACE].push_back(exch_site_templ);

// Request a list of face neighbors of the edge
         n_sect_parts_actual = tesselation.FaceNeighborsOfEdge(sector_division, edge, faces);

// Loop over participating slabs
         part = 0;
         for(is = 0; is < n_slab_parts_actual; is++) {

// Loop over participating sectors
            for(it = 0; it < n_sect_parts_actual; it++) {
               labels[part] = (slab + is) * tesselation.NFaces(sector_division) + faces[it];
               ranks[part] = block_ranks[labels[part]];
               part++;
            };
         };

// Initialize the exchange site and add it to the list if it is local to this process
         exch_sites[GEONBR_RFACE].back().SetUpProperties(sidx, n_slab_parts_actual * n_sect_parts_actual, buf_size, ranks, labels, mpi_config);
         if(exch_sites[GEONBR_RFACE].back().site_comm != MPI_COMM_NULL) exch_site_idx[GEONBR_RFACE].push_back(sidx);
         sidx++;
     };
  };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TEDGE
//----------------------------------------------------------------------------------------------------------------------------------------------------

// The TEDGE buffer size does not depend on the slab thickness; can use either thin or thick
   buf_size = block_temp_thin.GetBufferSize(GEONBR_TEDGE);
   n_sect_parts_actual = n_sect_parts[GEONBR_TEDGE];

// Loop on slab interfaces
   sidx = 0;
   for(slab_i = 0; slab_i <= n_slabs; slab_i++) {

// Boundaries
      if((slab_i == 0) || (slab_i == n_slabs)) n_slab_parts_actual = n_slab_parts[GEONBR_TFACE] - 1;
      else n_slab_parts_actual = n_slab_parts[GEONBR_TFACE];

// Loop on edges
      for(edge = 0; edge < tesselation.NEdges(sector_division); edge++) {
         exch_sites[GEONBR_TEDGE].push_back(exch_site_templ);

// Request a list of face neighbors of the edge
         n_sect_parts_actual = tesselation.FaceNeighborsOfEdge(sector_division, edge, faces);

// Loop over participating slabs
         part = 0;
         for(is = 0; is < n_slab_parts_actual; is++) {

// Loop over participating sectors
            for(it = 0; it < n_sect_parts_actual; it++) {
               if(slab_i == 0) labels[part] = faces[it];
               else labels[part] = (slab_i + is - 1) * tesselation.NFaces(sector_division) + faces[it];
               ranks[part] = block_ranks[labels[part]];
               part++;
            };
         };
            
// Initialize the exchange site and add it to the list if it is local to this process
         exch_sites[GEONBR_TEDGE].back().SetUpProperties(sidx, n_slab_parts_actual * n_sect_parts_actual, buf_size, ranks, labels, mpi_config);
         if(exch_sites[GEONBR_TEDGE].back().site_comm != MPI_COMM_NULL) exch_site_idx[GEONBR_TEDGE].push_back(sidx);
         sidx++;
      };
   };         

//----------------------------------------------------------------------------------------------------------------------------------------------------
// REDGE
//----------------------------------------------------------------------------------------------------------------------------------------------------

// This site has only one participating slab
   n_slab_parts_actual = n_slab_parts[GEONBR_REDGE];

// Loop on slabs
   sidx = 0;
   for(slab = 0; slab < n_slabs; slab++) {

// The buffer size depends on the slab thickness
      buf_size = (slab < n_thinner_slabs ? block_temp_thin.GetBufferSize(GEONBR_REDGE) : block_temp_thick.GetBufferSize(GEONBR_REDGE));

// Loop on vertices
      for(vert = 0; vert < tesselation.NVerts(sector_division); vert++) {
         exch_sites[GEONBR_REDGE].push_back(exch_site_templ);

// Request a list of face neighbors of the vertex
         n_sect_parts_actual = tesselation.FaceNeighborsOfVert(sector_division, vert, faces);

// Loop over participating slabs
         part = 0;
         for(is = 0; is < n_slab_parts_actual; is++) {

// Loop over participating sectors
            for(it = 0; it < n_sect_parts_actual; it++) {
               labels[part] = (slab + is) * tesselation.NFaces(sector_division) + faces[it];
               ranks[part] = block_ranks[labels[part]];
               part++;
            };
         };
            
// Initialize the exchange site and add it to the list if it is local to this process
         exch_sites[GEONBR_REDGE].back().SetUpProperties(sidx, n_slab_parts_actual * n_sect_parts_actual, buf_size, ranks, labels, mpi_config);
         if(exch_sites[GEONBR_REDGE].back().site_comm != MPI_COMM_NULL) exch_site_idx[GEONBR_REDGE].push_back(sidx);
         sidx++;
      };
   };         

//----------------------------------------------------------------------------------------------------------------------------------------------------
// VERTX
//----------------------------------------------------------------------------------------------------------------------------------------------------

// The VERTX buffer size does not depend on the slab thickness; can use either thin or thick
   buf_size = block_temp_thin.GetBufferSize(GEONBR_VERTX);
   n_sect_parts_actual = n_sect_parts[GEONBR_VERTX];

// Loop on slab interfaces
   sidx = 0;
   for(slab_i = 0; slab_i <= n_slabs; slab_i++) {

// Boundaries
      if((slab_i == 0) || (slab_i == n_slabs)) n_slab_parts_actual = n_slab_parts[GEONBR_VERTX] - 1;
      else n_slab_parts_actual = n_slab_parts[GEONBR_VERTX];

// Loop on vertices
      for(vert = 0; vert < tesselation.NVerts(sector_division); vert++) {
         exch_sites[GEONBR_VERTX].push_back(exch_site_templ);

// Request a list of face neighbors of the vertex
         n_sect_parts_actual = tesselation.FaceNeighborsOfVert(sector_division, vert, faces);

// Loop over participating slabs
         part = 0;
         for(is = 0; is < n_slab_parts_actual; is++) {

// Loop over participating sectors
            for(it = 0; it < n_sect_parts_actual; it++) {
               if(slab_i == 0) labels[part] = faces[it];
               else labels[part] = (slab_i + is - 1) * tesselation.NFaces(sector_division) + faces[it];
               ranks[part] = block_ranks[labels[part]];
               part++;
            };
         };
            
// Initialize the exchange site and add it to the list if it is local to this process
         exch_sites[GEONBR_VERTX].back().SetUpProperties(sidx, n_slab_parts_actual * n_sect_parts_actual, buf_size, ranks, labels, mpi_config);
         if(exch_sites[GEONBR_VERTX].back().site_comm != MPI_COMM_NULL) exch_site_idx[GEONBR_VERTX].push_back(sidx);
         sidx++;
      };
   };         
};

/*!
\author Vladimir Florinski
\date 07/02/2024
*/
template <typename datatype, typename blocktype>
void DomainPartition<datatype, blocktype>::ExportExchangeSites(void)
{
   int slab, sector, bidx, sidx, n_sites_sect, is, ie, iv;
   unsigned int block;
   int edges[verts_per_face], verts[verts_per_face];
   std::vector<ExchangeSite<ConservedVariables>*> exch_sites_in;

   for(block = 0; block < blocks_local.size(); block++) {
      bidx = blocks_local[block].GetIndex();
      slab = GetSlab(bidx);
      sector = GetSector(bidx);

// TFACE
      exch_sites_in.clear();
      for(is = 0; is < n_slab_parts[GEONBR_TFACE]; is++) {
         sidx = (slab + is) * tesselation.NFaces(sector_division) + sector;
         exch_sites_in.push_back(&exch_sites[GEONBR_TFACE][sidx]);
      };
      blocks_local[bidx].ImportExchangeSites(GEONBR_TFACE, exch_sites_in);

// RFACE
      exch_sites_in.clear();
      n_sites_sect = tesselation.EdgeNeighborsOfFace(sector_division, sector, edges);
      for(is = 0; is < n_slab_parts[GEONBR_RFACE]; is++) {
         for(ie = 0; ie < n_sites_sect; ie++) {
            sidx = (slab + is) * tesselation.NEdges(sector_division) + edges[ie];
            exch_sites_in.push_back(&exch_sites[GEONBR_RFACE][sidx]);
         };
      };
      blocks_local[bidx].ImportExchangeSites(GEONBR_RFACE, exch_sites_in);

// TEDGE
      exch_sites_in.clear();
      n_sites_sect = tesselation.EdgeNeighborsOfFace(sector_division, sector, edges);
      for(is = 0; is < n_slab_parts[GEONBR_TEDGE]; is++) {
         for(ie = 0; ie < n_sites_sect; ie++) {
            sidx = (slab + is) * tesselation.NEdges(sector_division) + edges[ie];
            exch_sites_in.push_back(&exch_sites[GEONBR_TEDGE][sidx]);
         };
      };
      blocks_local[bidx].ImportExchangeSites(GEONBR_TEDGE, exch_sites_in);

// REDGE
      exch_sites_in.clear();
      n_sites_sect = tesselation.VertNeighborsOfFace(sector_division, sector, verts);
      for(is = 0; is < n_slab_parts[GEONBR_REDGE]; is++) {
         for(iv = 0; iv < n_sites_sect; iv++) {
            sidx = (slab + is) * tesselation.NVerts(sector_division) + verts[iv];
            exch_sites_in.push_back(&exch_sites[GEONBR_REDGE][sidx]);
         };
      };
      blocks_local[bidx].ImportExchangeSites(GEONBR_REDGE, exch_sites_in);
      
// VERTX
      exch_sites_in.clear();
      n_sites_sect = tesselation.VertNeighborsOfFace(sector_division, sector, verts);
      for(is = 0; is < n_slab_parts[GEONBR_VERTX]; is++) {
         for(iv = 0; iv < n_sites_sect; iv++) {
            sidx = (slab + is) * tesselation.NVerts(sector_division) + verts[iv];
            exch_sites_in.push_back(&exch_sites[GEONBR_VERTX][sidx]);
         };
      };
      blocks_local[bidx].ImportExchangeSites(GEONBR_VERTX, exch_sites_in);
   };
};

/*!
\author Vladimir Florinski
\date 06/27/2024
*/
template <typename datatype, typename blocktype>
void DomainPartition<datatype, blocktype>::ExchangeAll(void)
{
   unsigned int idx;
   NeighborType ntype;

// All work is done by the exchange sites
   for(ntype = GEONBR_TFACE; ntype <= GEONBR_VERTX; GEO_INCR(ntype, NeighborType)) {
      for(idx = 0; idx < exch_site_idx[ntype].size(); idx++) {
         exch_sites[ntype][exch_site_idx[ntype][idx]].Exchange();
      };
   };
};

template class DomainPartition<ConservedVariables, BufferedBlock<VERTS_PER_FACE>>;

};
