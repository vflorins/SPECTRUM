/*!
\file domain_partition.hh
\brief Declares the class responsible for partitioning of the simuilation into blocks
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DOMAIN_PARTITION_HH
#define SPECTRUM_DOMAIN_PARTITION_HH

#include "geodesic/traversable_tesselation.hh"
#include "geodesic/buffered_block.hh"
#include "common/exchange_site.hh"

namespace Spectrum {

// TODO these should be set by configure
#define POLY_TYPE POLY_ICOSAHEDRON
#define MAX_FACE_DIVISION 5

#if (POLY_TYPE == POLY_HEXAHEDRON)
#define VERTS_PER_FACE 4
#else
#define VERTS_PER_FACE 3
#endif

#define GHOST_WIDTH 2
#define GHOST_HEIGHT 2

// TODO Put all SILO stuff in a derived class

#ifdef USE_SILO

//! The number of files for PMPIO
const int n_silofiles = 4;

//! The length, in characters,  of the time stamp
const int time_length = 3;

//! The length, in characters, of the file index
const int file_length = 3;

//! The length, in characters, of the directory index
const int dirc_length = 5;

//! File name prefix
const std::string datafile_base = "pltf";

//! The separator of file and time indices
const char ft_separator = '_';

//! File name extension
const std::string geofile_ext = ".silo";

//! Pattern for directory naming within each file for PMPIO
const std::string dirc_base = "piornk";

//! The assembled file name
const std::string assembled_base = "imgf";

//! The assembled mesh name
const std::string asmb_mesh_name = "geomesh";

//! Write tag for PMPIO
const int pmpio_wtag = 1001;

//! Read tag for PMPIO
const int pmpio_rtag = 1002;

//! Full path name length
const int path_length = 256;

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
template <typename datatype, typename blocktype>
class DomainPartition
{
protected:

//! Whether this is a parallel or a serial run
   bool is_parallel = false;

//! Number of vertices per face
   static constexpr int verts_per_face = VERTS_PER_FACE;

//! Number of edges meeting at a vertex
   static constexpr int edges_per_vert = EdgesAtVert<verts_per_face>();

//! Number of participating slabs per site (= number of exchange sites per slab)
   const int n_slab_parts[N_NBRTYPES] = {2, 1, 2, 1, 2};

//! Number of participating sectors per site
   const int n_sect_parts[N_NBRTYPES] = {1, 2, 2, edges_per_vert, edges_per_vert};

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

//! Global array of exchange sites
   std::vector<std::shared_ptr<ExchangeSite<datatype>>> exch_sites[N_NBRTYPES];

//! Local simulation blocks 
   std::vector<blocktype> blocks_local;

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

//! Set up block geometric properties
   void SetUpBlockGeometry(void);

//! Set up the exchange sites
   void SetUpExchangeSites(void);

//! Assign exchange site objects to blocks
   void ExportExchangeSites(void);

//! Perform the exchange
   void ExchangeAll(void);
};

/*!
\param[in] slab   Slab index
\param[in] sector Sector index
\return The global block index
*/
template <typename datatype, typename blocktype>
inline int DomainPartition<datatype, blocktype>::GetBlockIndex(int slab, int sector) const
{
   return slab * n_sectors + sector;
};

/*!
\param[in] bidx Block index
\return Slab index
*/
template <typename datatype, typename blocktype>
inline int DomainPartition<datatype, blocktype>::GetSlab(int bidx) const
{
   return bidx / n_slabs;
};

/*!
\param[in] bidx Block index
\return Sector index
*/
template <typename datatype, typename blocktype>
inline int DomainPartition<datatype, blocktype>::GetSector(int bidx) const
{
   return bidx % n_sectors;
};

};

#endif

