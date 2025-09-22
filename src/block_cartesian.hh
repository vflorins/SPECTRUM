/*!
\file block_cartesian.hh
\brief Defines a class to operate on uniform cartesian blocks
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BLOCK_CARTESIAN_HH
#define SPECTRUM_BLOCK_CARTESIAN_HH

#include "block_base.hh"

namespace Spectrum {

//! Number of variables per zone
const int n_variables_cartesian = 9;

//! Largest possible number of neighbors
const int max_neighbors_cartesian = 27;

//! Largest possible number of neighbor levels
const int max_neighbor_levels_cartesian = 27;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BlockCartesian class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A three-dimensional cuboid region partitioned into zones
\author Juan G Alonso Guzman

This class is designed to facilitate the interaction with the Cartesian library. It corresponds to a single block in a simulation. Information stored include block dimensions, neighbor list, and data. 
*/
template <typename HConfig_>
class BlockCartesian : public BlockBase<HConfig_> {
public:

   using HConfig = HConfig_;

protected:

//! Allocate memory for variables and neighbors
   void AllocateMemory(void);

public:

//! Default constructor
   BlockCartesian(void);

//! Destructor
   ~BlockCartesian() override = default;

//! Obtain the extents of the block from Cartesian
   void LoadDimensions(double unit_length_block) override;

//! Request all neighbors from Cartesian
   void LoadNeighbors(void) override;

//! Load all variables into the block from Cartesian
   void LoadVariables(void) override;

//! Return the LN zone center and offset for a given point
   MultiIndex GetQuadrant(const GeoVector& pos) const override;

//! Return the refinement level of a neighbor - multi-index version
   int GetNeighborLevel(const MultiIndex& level_idx) const override;

//! Return the refinement level of a neighbor - three index version
   int GetNeighborLevel(int i, int j, int k) const override;

//! Return the node of a neighbor - multi-index version
   int GetNeighborNode(const MultiIndex& node_idx) const override;

//! Return the node of a neighbor - three index version
   int GetNeighborNode(int i, int j, int k) const override;

//! Print variables in block
   void PrintVariables(void) const override;

#ifdef GEO_DEBUG

//! Print the indices of all neighbor blocks
   void PrintNeighbors(void) override;

//! Load coordinates into the first three variables (for testing only)
   void LoadCoords(void);

#endif

};

////! Block type
//#if SERVER_TYPE == SERVER_CARTESIAN
//typedef BlockCartesian BlockType;
//#endif

};

#endif
