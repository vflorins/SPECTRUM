/*!
\file domain_partition.hh
\brief Declares the class responsible for partitioning of the simuilation into blocks
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DOMAIN_PARTITION_HH
#define SPECTRUM_DOMAIN_PARTITION_HH

#include <memory>

#include <config.h>

#include <common/exchange_site.hh>
#include <geodesic/stenciled_block.hh>
#include <geodesic/traversable_tesselation.hh>
#include <geodesic/neighbors.hh>

namespace Spectrum {

// TODO these should be set by configure
#define POLY_TYPE POLY_ICOSAHEDRON
#define MAX_FACE_DIVISION 6
#define GHOST_WIDTH 2
#define GHOST_HEIGHT 2

constexpr int VERTS_PER_FACE = (POLY_TYPE == POLY_HEXAHEDRON ? 4 : 3);

#ifdef USE_SILO

//! The number of files for PMPIO
constexpr int n_silofiles = 4;

//! Write tag for PMPIO
constexpr int pmpio_wtag = 1001;

//! Read tag for PMPIO
constexpr int pmpio_rtag = 1002;

//! The length, in characters,  of the time stamp
constexpr int time_length = 3;

//! The length, in characters, of the file index
constexpr int file_length = 3;

//! The length, in characters, of the directory index
constexpr int dirc_length = 5;

//! Full path name length
constexpr int path_length = 256;

//! The separator of file and time indices
constexpr char ft_separator = '_';

//! File name prefix
const std::string datafile_base = "pltf";

//! File name extension
const std::string geofile_ext = ".silo";

//! Pattern for directory naming within each file for PMPIO
const std::string dirc_base = "piornk";

//! The assembled file name
const std::string assembled_base = "imgf";

//! The assembled mesh name
const std::string asmb_mesh_name = "geomesh";

//! The assembled variable name
const std::string asmb_var_name = "geovar";

//! The format for the time stamp
const std::string time_format = "%0" + std::to_string(time_length) + 'i';

//! The format for the file index
const std::string file_format = "%0" + std::to_string(file_length) + 'i';

//! The format for the directory index
const std::string dirc_format = "%0" + std::to_string(dirc_length) + 'i';

#endif

/*!
\brief Divide elements among almost groups with a difference in size of 0 or 1.
\param[in] n Number of elements
\param[in] m Number of groups
\return Elements in a smaller group (first), number of smaller groups (second)
*/
SPECTRUM_DEVICE_FUNC inline std::pair<int, int> DistributeEvenly(int n, int m)
{
   std::pair<int, int> res;
   res.first = n / m;
   res.second = m * (res.first + 1) - n;
   return res;
};

/*!
\brief A class controlling the partitioning of a simulation into computational units
\author Vladimir Florinski

This is a base class of a simulation control program. It assembles the mesh out of individual grid blocks and tells the blocks which data operations to perform. Only a single instance per process is permitted.
*/
template <typename blocktype>
class DomainPartition
{
protected:

//! Whether this is a parallel or a serial run
   bool is_parallel = false;

//! Number of vertices per face
   static constexpr int verts_per_face = VERTS_PER_FACE;

//! Number of edges meeting at a vertex
   static constexpr int edges_per_vert = EdgesAtVert<verts_per_face>();

//! Number of participating sectors per site
   static constexpr int n_sect_parts[N_NBRTYPES] = {1, 2, 2, edges_per_vert, edges_per_vert};

//! Sector division
   int sector_division;

//! Total number of slabs (global)
   int n_slabs;

//! Number of _thinner_ slabs (global)
   int n_thinner_slabs;

//! Total number of sectors on the sphere (global)
   int n_sectors;

//! Number of blocks in the mesh (global)
   int n_blocks;

//! Face division
   int face_division;

//! Number of r-shells in the mesh (global)
   int n_shells_global;

//! Number of t-faces in the mesh (global)
   int n_faces_global;

//! Size of the _thinner_ slab
   int thinner_slab_height;

//! Radii of all slab interfaces
   double* slab_interfaces = nullptr;

//! Ranks of processes owning each block
   int* block_ranks = nullptr;

//! Number of blocks in each rank
   int* blocks_in_rank = nullptr;

//! Global array of exchange sites. This works even if "datatype" is void because "ExchangeSite" is specialized to do nothing in that case.
   std::vector<std::shared_ptr<ExchangeSite<typename blocktype::datatype>>> exch_sites[N_NBRTYPES];

//! Local blocks
   std::vector<blocktype> blocks_local;

#if RUNTYPE == PARTICLES

   std::vector<blocktype>* blocks_ptrs;

   void GetBlockPtrs(void);

#endif

//! Local exchange sites indices
   std::vector<int> exch_sites_local[N_NBRTYPES];

//! Tesselation object
   TraversableTesselation<POLY_TYPE, MAX_FACE_DIVISION> tesselation;

//! Distance map object
   std::shared_ptr<DistanceBase> distance_map;

//! Default constructor
   DomainPartition(void) = default;

//! Compute the global index of the block from slab and sector
   int GetBlockIndex(int slab, int sector) const;

//! Compute the slab index from block index
   int GetSlab(int bidx) const;

//! Compute the sector index from block index
   int GetSector(int bidx) const;

//! Find the global block index containing the given vector
   int LocateBlock(const GeoVector& pos) const;

public:

//! Constructor with arguments
   DomainPartition(int n_shells, const std::shared_ptr<DistanceBase> dist_map, double face_div);

//! Copy constructor - deleted because only a single partition can exist
   DomainPartition(const DomainPartition& other) = delete;

//! Destructor
   ~DomainPartition();

//! Create the exchange sites (only enabled if "datatype" is available)
   template <typename datatype = blocktype::datatype, std::enable_if_t<!std::is_void<datatype>::value, bool> = true>
   void CreateExchangeSites(void);

//! Assign exchange site objects to blocks (only enabled if "datatype" is available)
   template <typename datatype = blocktype::datatype, std::enable_if_t<!std::is_void<datatype>::value, bool> = true>
   void ExportExchangeSites(void);

//! Perform the exchange (only enabled if "datatype" is available)
   template <typename datatype = blocktype::datatype, std::enable_if_t<!std::is_void<datatype>::value, bool> = true>
   void ExchangeAll(int test_block = -1);

//! Return the number of blocks
   int GetNBlocks(void) const {return blocks_local.size();};

//! Save the data
   void Save(int stamp);

#ifdef GEO_DEBUG

//! Print the block assignment to ranks
   void PrintTopology(void) const;

//! Print information about every exchange site
   template <typename datatype = blocktype::datatype, std::enable_if_t<!std::is_void<datatype>::value, bool> = true>
   void PrintExchSites(void) const;

//! Return reference to the block
   blocktype& GetBlock(int block) {return blocks_local[block];};

#endif
};

/*!
\param[in] slab   Slab index
\param[in] sector Sector index
\return The global block index
*/
template <typename blocktype>
inline int DomainPartition<blocktype>::GetBlockIndex(int slab, int sector) const
{
   return slab * n_sectors + sector;
};

/*!
\param[in] bidx Block index
\return Slab index
*/
template <typename blocktype>
inline int DomainPartition<blocktype>::GetSlab(int bidx) const
{
   return bidx / n_sectors;
};

/*!
\param[in] bidx Block index
\return Sector index
*/
template <typename blocktype>
inline int DomainPartition<blocktype>::GetSector(int bidx) const
{
   return bidx % n_sectors;
};

};

#endif
