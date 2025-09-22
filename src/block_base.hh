/*!
\file block_base.hh
\brief Defines a base class to operate on a 3D grid block
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BLOCK_BASE_HH
#define SPECTRUM_BLOCK_BASE_HH

//#include "config.h" // needed so that SERVER_TYPE is defined in order to typedef BlockType
#include "common/vectors.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BlockBase class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A three-dimensional region partitioned into zones
\author Vladimir Florinski

A topologically Cartesian three-dimensional hexahedral block identified by the "node" member. It uses linear storage model for compatibility with external solvers. A block stores its neighbor node indices and a set of physical variables at each interior zone.
*/
template <typename HConfig_>
class BlockBase {
public:

   using HConfig = HConfig_;

protected:

//! Number of variables (persistent)
   int n_variables;

//! Number of zones on each side of a block (persistent)
   MultiIndex block_size;

//! Number of ghost cells per side (persistent)
   int num_ghost_cells = 0;

//! Ratio of ghost cells per side to cells per dimension
   GeoVector ghost_to_phys_ratio;

//! Maximum number of neighbors per dimension (persistent)
   int max_neighbors_per_dim;

//! Maximum number of neighbors per dimension minus one (persistent)
   int max_neighbors_per_dim_mns;

//! Largest possible number of neighbors (persistent)
   int max_neighbors;

//! Largest possible number of neighbor levels (persistent)
   int max_neighbor_levels;

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
   double* variables = nullptr;

//! Default constructor (protected, class not designed to be instantiated)
   BlockBase(void) = default;

//! Copy constructor
   BlockBase(const BlockBase& other);

//! Move constructor
   BlockBase(BlockBase&& other);

//! Location of a variable in the "variables" linear array - multi-index version
   int DataLoc(int vidx, const MultiIndex& zone) const;

//! Location of a variable in the "variables" linear array - three index version
   int DataLoc(int vidx, int i, int j, int k) const;

//! Test if a zone is inside the block
   bool ZoneInsideBlock(const MultiIndex& zone) const;

public:

//! Destructor
   virtual ~BlockBase(void);

//! Assignment operator
   virtual BlockBase& operator =(const BlockBase& other);

//! Move assignment operator
   virtual BlockBase& operator =(BlockBase&& other);

//! Set number of ghost cells
   void SetGhostCells(int num_ghost_cells_in);

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

//! Obtain the extents of the block from external library
   virtual void LoadDimensions(double unit_length_block);

//! Request all neighbors from external library
   virtual void LoadNeighbors(void) = 0;

//! Load all variables into the block from external library
   virtual void LoadVariables(void) = 0;

//! Compute derived class data members
   void ConfigureProperties(void);

//! Return the zone size
   GeoVector GetZoneLength(void) const;

//! Return the number of zones in the block
   int GetZoneCount(void) const;

//! Return the block size
   MultiIndex GetBlockSize(void) const;

//! Return the index of the zone owning the position
   MultiIndex GetZone(const GeoVector& pos) const;

//! Return the LN zone center and offset for a given point
   void GetZoneOffset(const GeoVector& pos, MultiIndex& zone, GeoVector& offset) const;

//! Check if a position is inside the block
   bool PositionInside(const GeoVector& pos) const;

//! Return the LN zone center and offset for a given point
   virtual MultiIndex GetQuadrant(const GeoVector& pos) const = 0;

//! Retrieve a variable in a given zone
   double GetValue(const MultiIndex& zone, int vidx) const;

//! Retrieve a variable at the nearest cell center
   double GetValue(const GeoVector& pos, int vidx) const;

//! Return the max number of neighbors
   int GetNeighborCount(void) const;

//! Return the number of variables
   int GetVariableCount(void) const;

//! Return the address of "neighbor_nodes"
   int* GetNeighborNodesAddress(void);

//! Return the max number of neighbor levels
   int GetNeighborLevelCount(void) const;

//! Return the address of "neighbor_levels"
   int* GetNeighborLevelsAddress(void);

//! Return the refinement level of a neighbor - multi-index version
   virtual int GetNeighborLevel(const MultiIndex& level_idx) const = 0;

//! Return the refinement level of a neighbor - three index version
   virtual int GetNeighborLevel(int i, int j, int k) const = 0;

//! Return the node of a neighbor - multi-index version
   virtual int GetNeighborNode(const MultiIndex& node_idx) const = 0;

//! Return the node of a neighbor - three index version
   virtual int GetNeighborNode(int i, int j, int k) const = 0;

//! Return the address of "variables"
   double* GetVariablesAddress(void);

//! Print variables in block
   virtual void PrintVariables(void) const = 0;

//! Print the indices of all neighbor blocks
   void PrintNeighbors(void) requires (HConfig::build_type == BuildMode::debug) {};

};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BlockBase inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/08/2020
\param[in] node_in Node of this block
*/
template <typename HConfig>
inline void BlockBase<HConfig>::SetNode(int node_in)
{
   node = node_in;
};

/*!
\author Vladimir Florinski
\date 06/08/2020
\return The index (node) of the block
*/
template <typename HConfig>
inline int BlockBase<HConfig>::GetNode(void) const
{
   return node;
};

/*!
\author Vladimir Florinski
\date 01/25/2023
\return The address of "node"
*/
template <typename HConfig>
inline int* BlockBase<HConfig>::GetNodeAddress(void)
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
inline void BlockBase<HConfig>::SetDimensions(const GeoVector& corner1, const GeoVector& corner2)
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
inline GeoVector BlockBase<HConfig>::GetFaceMin(void) const
{
   return face_min;
};

/*!
\author Vladimir Florinski
\date 01/25/2023
\return Coordinates of the "outer" corner of the block
*/
template <typename HConfig>
inline GeoVector BlockBase<HConfig>::GetFaceMax(void) const
{
   return face_max;
};

/*!
\author Vladimir Florinski
\date 01/25/2023
\return The address of "face_min"
*/
template <typename HConfig>
inline GeoVector* BlockBase<HConfig>::GetFaceMinAddress(void)
{
   return &face_min;
};

/*!
\author Vladimir Florinski
\date 01/25/2023
\return The address of "face_max"
*/
template <typename HConfig>
inline GeoVector* BlockBase<HConfig>::GetFaceMaxAddress(void)
{
   return &face_max;
};

/*!
\author Juan G Alonso Guzman
\date 08/04/2023
\return Coordinates of the physical "inner" corner of the block
*/
template <typename HConfig>
inline GeoVector BlockBase<HConfig>::GetFaceMinPhys(void) const
{
   return face_min_phys;
};

/*!
\author Juan G Alonso Guzman
\date 08/04/2023
\return Coordinates of the physical "outer" corner of the block
*/
template <typename HConfig>
inline GeoVector BlockBase<HConfig>::GetFaceMaxPhys(void) const
{
   return face_max_phys;
};

/*!
\author Juan G Alonso Guzman
\date 08/04/2023
\return The address of "face_min_phys"
*/
template <typename HConfig>
inline GeoVector* BlockBase<HConfig>::GetFaceMinPhysAddress(void)
{
   return &face_min_phys;
};

/*!
\author Juan G Alonso Guzman
\date 08/04/2023
\return The address of "face_max_phys"
*/
template <typename HConfig>
inline GeoVector* BlockBase<HConfig>::GetFaceMaxPhysAddress(void)
{
   return &face_max_phys;
};

/*!
\author Vladimir Florinski
\date 06/08/2020
\return The size of a zone
*/
template <typename HConfig>
inline GeoVector BlockBase<HConfig>::GetZoneLength(void) const
{
   return zone_length;
};

/*!
\author Vladimir Florinski
\date 11/30/2022
\return The number of zones in the block
*/
template <typename HConfig>
inline int BlockBase<HConfig>::GetZoneCount(void) const
{
   return block_size.Prod();
};

/*!
\author Vladimir Florinski
\date 10/28/2022
\return The size of the block
*/
template <typename HConfig>
inline MultiIndex BlockBase<HConfig>::GetBlockSize(void) const
{
   return block_size;
};

/*!
\author Vladimir Florinski
\date 10/27/2022
return Number of elements in "neighbor_nodes"
*/
template <typename HConfig>
inline int BlockBase<HConfig>::GetNeighborCount(void) const
{
   return max_neighbors;
};

/*!
\author Vladimir Florinski
\date 10/28/2022
return Number of variables
*/
template <typename HConfig>
inline int BlockBase<HConfig>::GetVariableCount(void) const
{
   return n_variables;
};

/*!
\author Vladimir Florinski
\date 10/27/2022
return Pointer to the neighbor block storage
*/
template <typename HConfig>
inline int* BlockBase<HConfig>::GetNeighborNodesAddress(void)
{
   return neighbor_nodes;
};

/*!
\author Vladimir Florinski
\date 10/28/2022
return Number of elements in "neighbor_levels"
*/
template <typename HConfig>
inline int BlockBase<HConfig>::GetNeighborLevelCount(void) const
{
   return max_neighbor_levels;
};

/*!
\author Vladimir Florinski
\date 10/27/2022
return Pointer to the neighbor block storage
*/
template <typename HConfig>
inline int* BlockBase<HConfig>::GetNeighborLevelsAddress(void)
{
   return neighbor_levels;
};

/*!
\author Vladimir Florinski
\date 10/27/2022
return Pointer to the variables storage
*/
template <typename HConfig>
inline double* BlockBase<HConfig>::GetVariablesAddress(void)
{
   return variables;
};

};

#endif
