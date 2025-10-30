/*!
\file domain_partition.cc
\brief Implements the class responsible for partitioning of the simuilation into blocks
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "geodesic/domain_partition.hh"

#ifdef USE_SILO
#include <pmpio.h>
#endif

namespace Spectrum {

/*!
\author Vladimir Florinski
\date 03/21/2025
\param[in[ n_shells   Total number of shells
\param[in[ dist_map   Distance scaling function
\param[in[ face_div   Face division
*/
template <typename blocktype>
DomainPartition<blocktype>::DomainPartition(int n_shells, const std::shared_ptr<DistanceBase> dist_map, double face_div)
{
   int shells_per_slab, sector_width, bidx;
   double slab_height_fp;
   bool borders[2], corners[verts_per_face];

   n_shells_global = n_shells;
   face_division = face_div;
   distance_map = dist_map;

   n_faces_global = tesselation.NFaces(face_division);

// Perform domain decomposition based on sector size
   slab_height_fp = n_shells_global;
   sector_division = 0;
   sector_width = Pow2(face_division - sector_division);
   n_sectors = tesselation.NFaces(sector_division);
   n_slabs = 1;

// This adjusts the sector and slab sizes to find the optimal values (roughly equal number of zones in each dimension)
   while ((n_sectors * n_slabs < MPI_Config::n_workers) && (slab_height_fp >= min_height_to_ghost * GHOST_HEIGHT) && (sector_width >= min_width_to_ghost * GHOST_WIDTH)) {
      if (sector_width < slab_height_fp) {
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
   slab_interfaces[0] = 0.0;
   slab_interfaces[n_slabs] = 1.0;
   for (auto slab_i = 1; slab_i < n_slabs; slab_i++) {
      shells_per_slab = (slab_i <= n_thinner_slabs ? thinner_slab_height : thinner_slab_height + 1);
      slab_interfaces[slab_i] = slab_interfaces[slab_i - 1] + (double)shells_per_slab / (double)n_shells_global;
   };

   block_ranks = new int[n_blocks];
   blocks_in_rank = new int[MPI_Config::work_comm_size];
   memset(blocks_in_rank, 0x0, MPI_Config::work_comm_size * sizeof(int));
   int* face_index;
   int* vert_index;
   GeoVector* vcart;

// This is a loop on _all_ blocks in the simulation. Each process grabs blocks out of the pool in a round-robin fashion. The vector "blocks_local" emplaces empty blocks; the actual memory allocation happens in "SetDimensions()" and "AssociateMesh()".
   bidx = 0;
   for (auto slab = 0; slab < n_slabs; slab++) {
      shells_per_slab = (slab < n_thinner_slabs ? thinner_slab_height : thinner_slab_height + 1);
      borders[0] = (slab == n_slabs - 1 ? 1 : 0);
      borders[1] = (slab == 0 ? 1 : 0);

      for (auto sector = 0; sector < n_sectors; sector++) {
         block_ranks[bidx] = MakePeriodic(bidx, MPI_Config::work_comm_size);
         if (block_ranks[bidx] == MPI_Config::work_comm_rank) {
            blocks_local.emplace_back();
            blocks_local.back().SetIndex(bidx);
            blocks_local.back().SetDimensions(sector_width, GHOST_WIDTH, shells_per_slab, GHOST_HEIGHT, false);

// Inefficient to allocate/delete memory for each block, but we must obtain the array sizes from the block. In principle this can accommodate blocks of different sizes.
            face_index = new int[blocks_local.back().TotalFaces()];
            vert_index = new int[blocks_local.back().TotalVerts()];
            vcart = new GeoVector[blocks_local.back().TotalVerts()];

            tesselation.GetAllInsideFaceNative(sector_division, sector, face_division, GHOST_WIDTH, face_index, vert_index, corners);
            tesselation.FillVertCoordArrays(blocks_local.back().TotalVerts(), vert_index, vcart);
            blocks_local.back().AssociateMesh(slab_interfaces[slab], slab_interfaces[slab + 1], corners, borders, vcart, distance_map, false);

            delete[] face_index;
            delete[] vert_index;
            delete[] vcart;
         };
         blocks_in_rank[block_ranks[bidx++]]++;
      };
   };

// Print the mesh summary
   if(MPI_Config::is_master) {
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
\date 03/14/2025
*/
template <typename blocktype>
DomainPartition<blocktype>::~DomainPartition()
{
   delete[] slab_interfaces;
   delete[] block_ranks;
   delete[] blocks_in_rank;
};

/*!
\author Vladimir Florinski
\date 03/14/2025
\param[in] pos Vector position
\return The block global index, -1 if vector is outside the boundaries
*/
template <typename blocktype>
int DomainPartition<blocktype>::LocateBlock(const GeoVector& pos) const
{
   int slab, sect;
   
// Locate the slab
   slab = LocateInArray(0, n_slabs, slab_interfaces, pos.Norm(), false);
   if (slab < 0) return -1;

// Locate the sector using the tesselation
   sect = tesselation.Locate(sector_division, UnitVec(pos));
   return GetBlockIndex(slab, sect);
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 03/14/2025
*/
template <typename blocktype>
template <typename datatype, std::enable_if_t<!std::is_void<datatype>::value, bool>>
void DomainPartition<blocktype>::CreateExchangeSites(void)
{
   int sector_width, n_slab_parts_actual, n_sect_parts_actual, buf_size, sidx, part;

// Initialize two buffered blocks, one for thinner slabs and one for thicker slabs - we only need them to call GetBufferSize().
   sector_width = Pow2(face_division - sector_division);
   blocktype block_temp_thin, block_temp_thick;
   block_temp_thin .SetDimensions(sector_width, GHOST_WIDTH, thinner_slab_height    , GHOST_HEIGHT, false);
   block_temp_thick.SetDimensions(sector_width, GHOST_WIDTH, thinner_slab_height + 1, GHOST_HEIGHT, false);

// Make the arrays large enough for any site. Note that "labels" are basically global block indices.
   int sectors[edges_per_vert], ranks[2 * edges_per_vert], labels[2 * edges_per_vert];

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TFACE
//----------------------------------------------------------------------------------------------------------------------------------------------------

// The TFACE buffer size does not depend on the slab thickness; can use either thin or thick
   buf_size = block_temp_thin.GetBufferSize(GEONBR_TFACE);
   n_sect_parts_actual = n_sect_parts[GEONBR_TFACE];

// Loop on slab interfaces
   sidx = 0;
   for (auto slab_i = 0; slab_i <= n_slabs; slab_i++) {

// Boundaries
      n_slab_parts_actual = n_slab_parts[GEONBR_TFACE];
      if ((slab_i == 0) || (slab_i == n_slabs)) n_slab_parts_actual--;

// Loop on faces
      for (auto face = 0; face < tesselation.NFaces(sector_division); face++) {
         exch_sites[GEONBR_TFACE].emplace_back(std::make_shared<ExchangeSite<datatype>>());

// Loop over participating slabs
         part = 0;
         for (auto is = 0; is < n_slab_parts_actual; is++) {

// Loop over participating sectors
            for (auto it = 0; it < n_sect_parts_actual; it++) {
               if (slab_i == 0) labels[part] = face;
               else labels[part] = (slab_i + is - 1) * tesselation.NFaces(sector_division) + face;
               ranks[part] = block_ranks[labels[part]];
               part++;
            };
         };
            
// Initialize the exchange site and add it to the list if it is local to this process
         exch_sites[GEONBR_TFACE].back()->SetUpProperties(sidx, n_slab_parts_actual * n_sect_parts_actual, buf_size, ranks, labels);
         if (exch_sites[GEONBR_TFACE].back()->GetCommSize() > 0) exch_sites_local[GEONBR_TFACE].push_back(sidx);
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
   for (auto slab = 0; slab < n_slabs; slab++) {

// The buffer size depends on the slab thickness
      buf_size = (slab < n_thinner_slabs ? block_temp_thin.GetBufferSize(GEONBR_RFACE) : block_temp_thick.GetBufferSize(GEONBR_RFACE));

// Loop on edges
      for (auto edge = 0; edge < tesselation.NEdges(sector_division); edge++) {
         exch_sites[GEONBR_RFACE].emplace_back(std::make_shared<ExchangeSite<datatype>>());

// Request a list of face neighbors of the edge
         n_sect_parts_actual = tesselation.FaceNeighborsOfEdge(sector_division, edge, sectors);

// Loop over participating slabs
         part = 0;
         for (auto is = 0; is < n_slab_parts_actual; is++) {

// Loop over participating sectors
            for (auto it = 0; it < n_sect_parts_actual; it++) {
               labels[part] = (slab + is) * tesselation.NFaces(sector_division) + sectors[it];
               ranks[part] = block_ranks[labels[part]];
               part++;
            };
         };

// Initialize the exchange site and add it to the list if it is local to this process
         exch_sites[GEONBR_RFACE].back()->SetUpProperties(sidx, n_slab_parts_actual * n_sect_parts_actual, buf_size, ranks, labels);
         if (exch_sites[GEONBR_RFACE].back()->GetCommSize() > 0) exch_sites_local[GEONBR_RFACE].push_back(sidx);
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
   for (auto slab_i = 0; slab_i <= n_slabs; slab_i++) {

// Boundaries
      if ((slab_i == 0) || (slab_i == n_slabs)) n_slab_parts_actual = n_slab_parts[GEONBR_TFACE] - 1;
      else n_slab_parts_actual = n_slab_parts[GEONBR_TFACE];

// Loop on edges
      for (auto edge = 0; edge < tesselation.NEdges(sector_division); edge++) {
         exch_sites[GEONBR_TEDGE].emplace_back(std::make_shared<ExchangeSite<datatype>>());

// Request a list of face neighbors of the edge
         n_sect_parts_actual = tesselation.FaceNeighborsOfEdge(sector_division, edge, sectors);

// Loop over participating slabs
         part = 0;
         for (auto is = 0; is < n_slab_parts_actual; is++) {

// Loop over participating sectors
            for (auto it = 0; it < n_sect_parts_actual; it++) {
               if (slab_i == 0) labels[part] = sectors[it];
               else labels[part] = (slab_i + is - 1) * tesselation.NFaces(sector_division) + sectors[it];
               ranks[part] = block_ranks[labels[part]];
               part++;
            };
         };
            
// Initialize the exchange site and add it to the list if it is local to this process
         exch_sites[GEONBR_TEDGE].back()->SetUpProperties(sidx, n_slab_parts_actual * n_sect_parts_actual, buf_size, ranks, labels);
         if (exch_sites[GEONBR_TEDGE].back()->GetCommSize() > 0) exch_sites_local[GEONBR_TEDGE].push_back(sidx);
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
   for (auto slab = 0; slab < n_slabs; slab++) {

// The buffer size depends on the slab thickness
      buf_size = (slab < n_thinner_slabs ? block_temp_thin.GetBufferSize(GEONBR_REDGE) : block_temp_thick.GetBufferSize(GEONBR_REDGE));

// Loop on vertices
      for (auto vert = 0; vert < tesselation.NVerts(sector_division); vert++) {
         exch_sites[GEONBR_REDGE].emplace_back(std::make_shared<ExchangeSite<datatype>>());

// Request a list of face neighbors of the vertex
         n_sect_parts_actual = tesselation.FaceNeighborsOfVert(sector_division, vert, sectors);

// Loop over participating slabs
         part = 0;
         for (auto is = 0; is < n_slab_parts_actual; is++) {

// Loop over participating sectors
            for (auto it = 0; it < n_sect_parts_actual; it++) {
               labels[part] = (slab + is) * tesselation.NFaces(sector_division) + sectors[it];
               ranks[part] = block_ranks[labels[part]];
               part++;
            };
         };
            
// Initialize the exchange site and add it to the list if it is local to this process
         exch_sites[GEONBR_REDGE].back()->SetUpProperties(sidx, n_slab_parts_actual * n_sect_parts_actual, buf_size, ranks, labels);
         if (exch_sites[GEONBR_REDGE].back()->GetCommSize() > 0) exch_sites_local[GEONBR_REDGE].push_back(sidx);
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
   for (auto slab_i = 0; slab_i <= n_slabs; slab_i++) {

// Boundaries
      if ((slab_i == 0) || (slab_i == n_slabs)) n_slab_parts_actual = n_slab_parts[GEONBR_VERTX] - 1;
      else n_slab_parts_actual = n_slab_parts[GEONBR_VERTX];

// Loop on vertices
      for (auto vert = 0; vert < tesselation.NVerts(sector_division); vert++) {
         exch_sites[GEONBR_VERTX].emplace_back(std::make_shared<ExchangeSite<datatype>>());

// Request a list of face neighbors of the vertex
         n_sect_parts_actual = tesselation.FaceNeighborsOfVert(sector_division, vert, sectors);

// Loop over participating slabs
         part = 0;
         for (auto is = 0; is < n_slab_parts_actual; is++) {

// Loop over participating sectors
            for (auto it = 0; it < n_sect_parts_actual; it++) {
               if (slab_i == 0) labels[part] = sectors[it];
               else labels[part] = (slab_i + is - 1) * tesselation.NFaces(sector_division) + sectors[it];
               ranks[part] = block_ranks[labels[part]];
               part++;
            };
         };
            
// Initialize the exchange site and add it to the list if it is local to this process
         exch_sites[GEONBR_VERTX].back()->SetUpProperties(sidx, n_slab_parts_actual * n_sect_parts_actual, buf_size, ranks, labels);
         if (exch_sites[GEONBR_VERTX].back()->GetCommSize() > 0) exch_sites_local[GEONBR_VERTX].push_back(sidx);
         sidx++;
      };
   };         
};

/*!
\author Vladimir Florinski
\date 07/02/2024
*/
template <typename blocktype>
template <typename datatype, std::enable_if_t<!std::is_void<datatype>::value, bool>>
void DomainPartition<blocktype>::ExportExchangeSites(void)
{
   int slab, sector, sidx, n_sites_sect;
   int edges[verts_per_face], verts[verts_per_face];
   std::vector<std::shared_ptr<ExchangeSite<datatype>>> exch_sites_block;

// This loop generates an array "exch_sites_block" for each block in this process from the global "exch_sites" array.
   for (auto& block : blocks_local) {
      slab = GetSlab(block.GetIndex());
      sector = GetSector(block.GetIndex());

// TFACE
      exch_sites_block.clear();
      for (auto is = 0; is < n_slab_parts[GEONBR_TFACE]; is++) {
         sidx = (slab + is) * tesselation.NFaces(sector_division) + sector;
         exch_sites_block.push_back(exch_sites[GEONBR_TFACE][sidx]);
      };
      block.ImportExchangeSites(GEONBR_TFACE, exch_sites_block);

// RFACE and TEDGE
      n_sites_sect = tesselation.EdgeNeighborsOfFace(sector_division, sector, edges);
      for (auto ntype : {GEONBR_RFACE, GEONBR_TEDGE}) {
         exch_sites_block.clear();
         for (auto is = 0; is < n_slab_parts[ntype]; is++) {
            for (auto ie = 0; ie < n_sites_sect; ie++) {
               sidx = (slab + is) * tesselation.NEdges(sector_division) + edges[ie];
               exch_sites_block.push_back(exch_sites[ntype][sidx]);
            };
         };
         block.ImportExchangeSites(ntype, exch_sites_block);
      };

// REDGE and VERTX
      n_sites_sect = tesselation.VertNeighborsOfFace(sector_division, sector, verts);
      for (auto ntype : {GEONBR_REDGE, GEONBR_VERTX}) {
         exch_sites_block.clear();
         for (auto is = 0; is < n_slab_parts[ntype]; is++) {
            for (auto iv = 0; iv < n_sites_sect; iv++) {
               sidx = (slab + is) * tesselation.NVerts(sector_division) + verts[iv];
               exch_sites_block.push_back(exch_sites[ntype][sidx]);
            };
         };
         block.ImportExchangeSites(ntype, exch_sites_block);
      };
   };
};

/*!
\author Vladimir Florinski
\date 06/27/2024
*/
template <typename blocktype>
template <typename datatype, std::enable_if_t<!std::is_void<datatype>::value, bool>>
void DomainPartition<blocktype>::ExchangeAll(int test_block)
{
   for (auto ntype = GEONBR_TFACE; ntype <= GEONBR_VERTX; GEO_INCR(ntype, NeighborType)) {

// Pack the buffers - loop location should allow for some parallelism
      for (auto& block : blocks_local) block.PackBuffers(ntype, test_block);
      
// Perform the exchange
      for (auto sidx : exch_sites_local[ntype]) exch_sites[ntype][sidx]->Exchange();

// Unpack the buffers - test with a different loop location
      for (auto& block : blocks_local) block.UnPackBuffers(ntype, test_block);
   };
};

#ifdef USE_SILO

/*!
\author Vladimir Florinski
\date 03/14/2025
*/
template <typename blocktype>
void DomainPartition<blocktype>::Save(int stamp)
{
// Character strings holding the numerical values of the time slice, file, directory, and block indices
   char time_num[time_length + 1];
   char file_num[file_length + 1];
   char dirc_num[dirc_length + 1];
   std::string file_name, dirc_name;

// Initialize the SILO parallel subsystem
   DBfile* data_file;
   PMPIO_baton_t* bat = PMPIO_Init(n_silofiles, PMPIO_WRITE, MPI_Config::work_comm, pmpio_wtag, PMPIO_DefaultCreate, PMPIO_DefaultOpen, PMPIO_DefaultClose, NULL);

// Compose the file and directory names based on our PMPIO group and time stamp
   sprintf(time_num, time_format.c_str(), stamp);
   sprintf(file_num, file_format.c_str(), PMPIO_GroupRank(bat, MPI_Config::work_comm_rank));
   sprintf(dirc_num, dirc_format.c_str(), PMPIO_RankInGroup(bat, MPI_Config::work_comm_rank));
   file_name = datafile_base + file_num + ft_separator + time_num + geofile_ext;
   dirc_name = dirc_base + dirc_num;

   data_file = (DBfile*)PMPIO_WaitForBaton(bat, file_name.c_str(), dirc_name.c_str());

// Each process writes its blocks sequentially
   for (auto& block : blocks_local) block.WriteSilo(data_file, false);
   PMPIO_HandOffBaton(bat, data_file);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Generate the SILO assembly file
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Name arrays used for both multimeshes and multivars
   DBfile* assembly_file;
   char** mmvnames;
   int* mmvtypes;
   char blck_num[blck_length + 1];
   std::string assembled_name, long_name;
   
   if (MPI_Config::is_master) {
      mmvnames = new char*[n_blocks];
      for (auto bidx = 0; bidx < n_blocks; bidx++) mmvnames[bidx] = new char[path_length];
      mmvtypes = new int[n_blocks];

// Form block, directory, and file names
      for (auto bidx = 0; bidx < n_blocks; bidx++) {
         sprintf(file_num, file_format.c_str(), PMPIO_GroupRank(bat, block_ranks[bidx]));
         sprintf(dirc_num, dirc_format.c_str(), PMPIO_RankInGroup(bat, block_ranks[bidx]));
         sprintf(blck_num, blck_format.c_str(), bidx);
         file_name = datafile_base + file_num + ft_separator + time_num + geofile_ext;
         dirc_name = dirc_base + dirc_num;
         long_name = file_name + ":/" + dirc_name + "/" + um_base + "_" + blck_num;
         strcpy(mmvnames[bidx], long_name.c_str());
         mmvtypes[bidx] = DB_UCDMESH;
      };

      assembled_name = assembled_base + time_num + geofile_ext;
      assembly_file = DBCreate(assembled_name.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
      DBPutMultimesh(assembly_file, asmb_mesh_name.c_str(), n_blocks, mmvnames, mmvtypes, NULL);

// Form block, directory, and file names
      for (auto bidx = 0; bidx < n_blocks; bidx++) {
         sprintf(file_num, file_format.c_str(), PMPIO_GroupRank(bat, block_ranks[bidx]));
         sprintf(dirc_num, dirc_format.c_str(), PMPIO_RankInGroup(bat, block_ranks[bidx]));
         sprintf(blck_num, blck_format.c_str(), bidx);
         file_name = datafile_base + file_num + ft_separator + time_num + geofile_ext;
         dirc_name = dirc_base + dirc_num;
         long_name = file_name + ":/" + dirc_name + "/" + um_base + "_" + blck_num;
         strcpy(mmvnames[bidx], long_name.c_str());
         mmvtypes[bidx] = DB_UCDVAR;
      };

      DBPutMultivar(assembly_file, asmb_var_name.c_str(), n_blocks, mmvnames, mmvtypes, NULL);
      DBClose(assembly_file);

      for (auto bidx = 0; bidx < n_blocks; bidx++) delete[] mmvnames[bidx];
      delete[] mmvnames;
      delete[] mmvtypes;
   };

// Close the PMPIO subsystem. One could in principle do this right after the blocks have been written, but then we would not have the PMPIO ranks and groups that are used to name the meshes and variables that are entered into the assembly file.
   PMPIO_Finish(bat);
};

#endif

#ifdef GEO_DEBUG

/*!
\author Vladimir Florinski
\date 03/16/2025
*/
template <typename blocktype>
void DomainPartition<blocktype>::PrintTopology(void) const
{
   for (auto proc = 0; proc < MPI_Config::work_comm_size; proc++) {
      if (proc == MPI_Config::work_comm_rank) {
         std::cerr << "Printing process " << proc << " topology\n";
         for (auto& block : blocks_local) {
            std::cerr << "   Block index " << block.GetIndex() << std::endl;
         };
      };
      MPI_Barrier(MPI_Config::work_comm);
   };
};

/*!
\author Vladimir Florinski
\date 03/20/2025
*/
template <typename blocktype>
template <typename datatype, std::enable_if_t<!std::is_void<datatype>::value, bool>>
void DomainPartition<blocktype>::PrintExchSites(void) const
{
   if (!MPI_Config::is_master) return;
   for (auto ntype = GEONBR_TFACE; ntype <= GEONBR_VERTX; GEO_INCR(ntype, NeighborType)) { 
      std::cerr << "Printing exchange sites of type " << neighbor_text[ntype] << ":\n";
      for (const auto& site : exch_sites[ntype]) {
         std::cerr << "   Site " << std::setw(5) << site->GetIndex() << ": ";
         std::cerr << std::setw(5) << site->GetPartCount() << " parts ";
         std::cerr << std::setw(5) << site->GetCommSize() << " processes\n";
      };

// Unpack the buffers - test with a different loop location
      for (auto block = 0; block < blocks_local.size(); block++) blocks_local[block].UnPackBuffers(ntype);
   };
};

//template void DomainPartition<BufferedBlock<VERTS_PER_FACE, int>>::PrintExchSites(void) const;

#endif

//template class DomainPartition<GridBlock<VERTS_PER_FACE>>;
template class DomainPartition<StenciledBlock<VERTS_PER_FACE>>;
//template class DomainPartition<BufferedBlock<VERTS_PER_FACE, int>>;
//template void DomainPartition<BufferedBlock<VERTS_PER_FACE, int>>::CreateExchangeSites(void);
//template void DomainPartition<BufferedBlock<VERTS_PER_FACE, int>>::ExportExchangeSites(void);
//template void DomainPartition<BufferedBlock<VERTS_PER_FACE, int>>::ExchangeAll(int test_block);

};
