/*!
\file server_types.hh
\brief Defines data structures for server-worker communication in a 3D grid computational environment
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SERVER_TYPES_HH
#define SPECTRUM_SERVER_TYPES_HH

#include "common/mpi_config.hh"
#include "common/vectors.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RequestInfo class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\brief Class to handle non-blocking receives
\author Juan G Alonso Guzman

This class houses all the variables necessary to implement non-blocking communications
*/
struct RequestInfo
{

//! MPI Request array
   MPI_Request* mpi_req = nullptr;

//! CPU rank array of processes that have completed requests
   int* cpu_rank = nullptr;

//! Count of CPUs that have completed requests
   int count = 0;

//! Constructor with a parameter
   RequestInfo(int size);

//! Destructor
   ~RequestInfo();
};


//----------------------------------------------------------------------------------------------------------------------------------------------------
// Inquiry class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\brief Data inquiry type
\author Vladimir Florinski
*/
struct Inquiry {

//! Type of inquiry: "0" is by ID, "1" is by position
   int type;

//! Node index if requesting by ID
   int node;

//! Position if requesting by position
   GeoVector pos;
};


//----------------------------------------------------------------------------------------------------------------------------------------------------
// Stencil class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\brief Interpolation stencil for 3D grid topology
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
*/
template <int stencil_n_elements_>
struct Stencil {

//! Number of elements in the stencil
   static constexpr int n_elements = stencil_n_elements_;

//! List of blocks (can also be used to store global nodes)
   int blocks[n_elements];

//! List of zones
   MultiIndex zones[n_elements];

//! Weights
   double weights[n_elements];

//! Partial derivatives of weights
   double derivatives[3*n_elements];

//! Print the list of blocks, zones, and weights
   void Print(void) const;

};



//----------------------------------------------------------------------------------------------------------------------------------------------------
// Block class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A three-dimensional cuboid region partitioned into zones
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

A topologically Cartesian three-dimensional hexahedral block identified by the "node" member. It uses linear storage model for compatibility with external solvers. A block stores its neighbor node indices and a set of physical variables at each interior zone.
*/
template <typename HConfig_>
class Block {
public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;

   // todo DataFields for DataBackground types
   using DataFields = Config::DataFields;

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/30/2025
*/
   static constexpr MultiIndex compute_block_size() {
      MultiIndex out;
      if (num_ghost_cells > 0)
         out = Config::block_size + 2 * num_ghost_cells;
      return out;
   };

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/30/2025
*/
   static constexpr GeoVector compute_ghost_to_phys_ratio() {
      GeoVector out;
      out[0] = (double)num_ghost_cells / (double)block_size.i;
      out[1] = (double)num_ghost_cells / (double)block_size.j;
      out[2] = (double)num_ghost_cells / (double)block_size.k;
      return out;
   }


/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/21/2025
\param[in] zone Zone's multi-index
\return Location in the array
*/
   static constexpr int fields_index(const MultiIndex& zone) {
      return zone.i + block_size.i * (zone.j + block_size.j * zone.k);
   }


/*!
\author Vladimir Florinski
\date 06/03/2020
\param[in] zone Zone's multi-index
\return True if the zone is in the block
*/
   static constexpr bool ZoneInside(const MultiIndex& zone) {
      return !((zone.i < 0) || (zone.i >= block_size.i)
          || (zone.j < 0) || (zone.j >= block_size.j)
          || (zone.k < 0) || (zone.k >= block_size.k));
   };

public:

   //! Number of ghost cells per side (persistent)
   static constexpr int num_ghost_cells = Config::num_ghost_cells;

//! Number of zones on each side of a block (persistent)
// todo multiindex config type
   static constexpr MultiIndex block_size = compute_block_size();

//! Ratio of ghost cells per side to cells per dimension
   static constexpr GeoVector ghost_to_phys_ratio = compute_ghost_to_phys_ratio();

//! Maximum number of neighbors per dimension (persistent)
   static constexpr int max_neighbors_per_dim = Config::max_neighbors_per_dim;

//! Maximum number of neighbors per dimension minus one (persistent)
   static constexpr int max_neighbors_per_dim_mns = Config::max_neighbors_per_dim - 1;

//! Largest possible number of neighbors (persistent)
   static constexpr int max_neighbors = Config::max_neighbors;

//! Largest possible number of neighbor levels (persistent)
   static constexpr int max_neighbor_levels = Config::max_neighbor_levels;

public:
   /*
    * TODO final review for constexpr
    *
    */

   //! Node index
   int node = -1;

//! Vector to the corner of the zone with minimum coordinates
   GeoVector face_min;

//! Vector to the corner of the zone with maximum coordinates
   GeoVector face_max;

//! Vector to the corner of the physical zone with minimum coordinates
   GeoVector face_min_phys;

//! Vector to the corner of the physical zone with maximum coordinates
   GeoVector face_max_phys;

//! Vector to the zone center closest to face_min
   GeoVector cent_min;

//! Vector to the zone center closest to face_max
   GeoVector cent_max;

//! Vector to the center of the block
   GeoVector center;

//! Zone size in code units
   GeoVector zone_length;

//! Neighbor indices in a linear array
   int* neighbor_nodes = nullptr;

//! Neighbor refinement levels in a linear array
   int* neighbor_levels = nullptr;

//! Variables in a linear array
   DataFields* fields = nullptr;

public:

//! Default constructor
   Block(void);

//! Destructor
   ~Block();

   //! Copy constructor
   Block(const Block& other);

//! Move constructor
   Block(Block&& other) noexcept;

   //! Assignment operator
   Block& operator =(const Block& other);

//! Move assignment operator
   Block& operator =(Block&& other);

protected:

//! Compute class data members
// todo probably not constexpr bc of possibility of block location/dims coming from data reader
   void ConfigureProperties(void);

public:

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/19/2025
\param[in] zone Zone
\return The fields datatype found at the zone
*/
   template <typename HConfig>
   DataFields operator[](const MultiIndex& zone) const
   {
      return fields[fields_index(zone)];
   };

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/19/2025
\param[in] pos  Position
\return The fields datatype found at the position
*/
   template <typename HConfig>
   DataFields operator[](const GeoVector& pos) const
   {
      MultiIndex zone = GetZone(pos);
      return fields[fields_index(zone)];
   };

public: // inline:

//! Set the block's node
   void SetNode(int node_in);

//! Return the block's index
   int GetNode(void) const;

//! Return the address of "node" (for MPI)
   int* GetNodeAddress(void);

//! Set up the corners of the block
   void SetDimensions(const GeoVector& corner1, const GeoVector& corner2);

//! Return the vector to one of the corner
   GeoVector GetFaceMin(void) const;

//! Return the vector to the opposite corner
   GeoVector GetFaceMax(void) const;

//! Return the address of "face_min" (for MPI)
   GeoVector* GetFaceMinAddress(void);

//! Return the address of "face_max" (for MPI)
   GeoVector* GetFaceMaxAddress(void);

//! Return the vector to one of the corner
   GeoVector GetFaceMinPhys(void) const;

//! Return the vector to the opposite corner
   GeoVector GetFaceMaxPhys(void) const;

//! Return the address of "face_min" (for MPI)
   GeoVector* GetFaceMinPhysAddress(void);

//! Return the address of "face_max" (for MPI)
   GeoVector* GetFaceMaxPhysAddress(void);

//! Return the zone size
   GeoVector GetZoneLength(void) const;

//! Return the number of zones in the block
   constexpr int GetZoneCount(void) const;

//! Return the block size
   MultiIndex constexpr GetBlockSize(void);

//! Return the max number of neighbors
   int GetNeighborCount(void) const;

//! Return the address of "fields"
//   double* GetVariablesAddress(void);

//! Return the number of fields
//   int GetVariableCount(void) const;

//! Return the address of "neighbor_nodes"
   int* GetNeighborNodesAddress(void);

//! Return the max number of neighbor levels
   int GetNeighborLevelCount(void) const;

//! Return the address of "neighbor_levels"
   int* GetNeighborLevelsAddress(void);

public: // API:

//! Return the index of the zone owning the position
   MultiIndex GetZone(const GeoVector& pos) const;

   //! Return the LN zone center and offset for a given point
   MultiIndex GetQuadrant(const GeoVector& pos) const;

//! Return the refinement level of a neighbor - multi-index version
   int GetNeighborLevel(const MultiIndex& level_idx) const;

//! Return the node of a neighbor - multi-index version
   int GetNeighborNode(const MultiIndex& node_idx) const;

//! Return the LN zone center and offset for a given point
   void GetZoneOffset(const GeoVector& pos, MultiIndex& zone, GeoVector& offset) const;

//! Check if a position is inside the block
   bool PositionInside(const GeoVector& pos) const;

   //! Load the dimensions of the block based on face_max, face_min
   void LoadDimensions(double unit_length_block);

   //! Print fields in block
   void PrintFields(void) const requires (HConfig::build_type == Config::BuildMode::debug);

//! Print the indices of all neighbor blocks
   void PrintNeighbors(void) const requires (HConfig::build_type == Config::BuildMode::debug);

//! Load coordinates into the first three fields (for testing only)
   void LoadCoords(void) requires (HConfig::build_type == Config::BuildMode::debug);

};



//----------------------------------------------------------------------------------------------------------------------------------------------------
// Block inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/08/2020
\param[in] node_in Node of this block
*/
template <typename HConfig>
inline void Block<HConfig>::SetNode(int node_in)
{
   node = node_in;
};

/*!
\author Vladimir Florinski
\date 06/08/2020
\return The index (node) of the block
*/
template <typename HConfig>
inline int Block<HConfig>::GetNode(void) const
{
   return node;
};

/*!
\author Vladimir Florinski
\date 01/25/2023
\return The address of "node"
*/
template <typename HConfig>
inline int* Block<HConfig>::GetNodeAddress(void)
{
   return &node;
};

/*!
\author Vladimir Florinski
\date 01/25/2023
\param[in] corner1 Coordinates of the "inner" corner of the block
\param[in] corner2 Coordinates of the "outer" corner of the block
*/
template <typename HConfig>
inline void Block<HConfig>::SetDimensions(const GeoVector& corner1, const GeoVector& corner2)
{
   face_min = corner1;
   face_max = corner2;
};

/*!
\author Vladimir Florinski
\date 01/25/2023
\return Coordinates of the "inner" corner of the block
*/
template <typename HConfig>
inline GeoVector Block<HConfig>::GetFaceMin(void) const
{
   return face_min;
};

/*!
\author Vladimir Florinski
\date 01/25/2023
\return Coordinates of the "outer" corner of the block
*/
template <typename HConfig>
inline GeoVector Block<HConfig>::GetFaceMax(void) const
{
   return face_max;
};

/*!
\author Vladimir Florinski
\date 01/25/2023
\return The address of "face_min"
*/
template <typename HConfig>
inline GeoVector* Block<HConfig>::GetFaceMinAddress(void)
{
   return &face_min;
};

/*!
\author Vladimir Florinski
\date 01/25/2023
\return The address of "face_max"
*/
template <typename HConfig>
inline GeoVector* Block<HConfig>::GetFaceMaxAddress(void)
{
   return &face_max;
};

/*!
\author Juan G Alonso Guzman
\date 08/04/2023
\return Coordinates of the physical "inner" corner of the block
*/
template <typename HConfig>
inline GeoVector Block<HConfig>::GetFaceMinPhys(void) const
{
   return face_min_phys;
};

/*!
\author Juan G Alonso Guzman
\date 08/04/2023
\return Coordinates of the physical "outer" corner of the block
*/
template <typename HConfig>
inline GeoVector Block<HConfig>::GetFaceMaxPhys(void) const
{
   return face_max_phys;
};

/*!
\author Juan G Alonso Guzman
\date 08/04/2023
\return The address of "face_min_phys"
*/
template <typename HConfig>
inline GeoVector* Block<HConfig>::GetFaceMinPhysAddress(void)
{
   return &face_min_phys;
};

/*!
\author Juan G Alonso Guzman
\date 08/04/2023
\return The address of "face_max_phys"
*/
template <typename HConfig>
inline GeoVector* Block<HConfig>::GetFaceMaxPhysAddress(void)
{
   return &face_max_phys;
};

/*!
\author Vladimir Florinski
\date 06/08/2020
\return The size of a zone
*/
template <typename HConfig>
inline GeoVector Block<HConfig>::GetZoneLength(void) const
{
   return zone_length;
};

/*!
\author Vladimir Florinski
\date 11/30/2022
\return The number of zones in the block
*/
template <typename HConfig>
inline constexpr int Block<HConfig>::GetZoneCount(void) const
{
   return block_size.Prod();
};

/*!
\author Vladimir Florinski
\date 10/28/2022
\return The size of the block
*/
template <typename HConfig>
inline constexpr MultiIndex Block<HConfig>::GetBlockSize(void)
{
   return block_size;
};

/*!
\author Vladimir Florinski
\date 10/27/2022
return Number of elements in "neighbor_nodes"
*/
template <typename HConfig>
inline int Block<HConfig>::GetNeighborCount(void) const
{
   return max_neighbors;
};

/*!
\author Vladimir Florinski
\date 10/27/2022
return Pointer to the neighbor block storage
*/
template <typename HConfig>
inline int* Block<HConfig>::GetNeighborNodesAddress(void)
{
   return neighbor_nodes;
};

/*!
\author Vladimir Florinski
\date 10/28/2022
return Number of elements in "neighbor_levels"
*/
template <typename HConfig>
inline int Block<HConfig>::GetNeighborLevelCount(void) const
{
   return max_neighbor_levels;
};

/*!
\author Vladimir Florinski
\date 10/27/2022
return Pointer to the neighbor block storage
*/
template <typename HConfig>
inline int* Block<HConfig>::GetNeighborLevelsAddress(void)
{
   return neighbor_levels;
};


};

#include "server_types.cc"

#endif
