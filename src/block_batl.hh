/*!
\file block_batl.hh
\brief Defines a class to operate on BATL blocks
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BLOCK_BATL_HH
#define SPECTRUM_BLOCK_BATL_HH

#include "block_cartesian.hh"

namespace Spectrum {

//! Number of variables per zone
const int n_variables_batl = 10;

//! Size of the block in each dimension
const MultiIndex block_size_batl(8, 8, 8);

//! Number of neighbors per dimension (depends only on the refinement ratio)
const int max_neighbors_per_dim_batl = 4;

//! Largest possible number of neighbors (depends only on the refinement ratio)
const int max_neighbors_batl = 64;

//! Largest possible number of neighbor levels (depends only on the refinement ratio)
const int max_neighbor_levels_batl = 27;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Interface for the Fortran routines in "spectrum_interface.f90"
//----------------------------------------------------------------------------------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif

/*!
\brief Obtain the pointers to the neighbor arrays for a node
\param[in]  node            Node
\param[out] neighbor_nodes  Entry into the BATL array of neighbor nodes
\param[out] neighbor_levels Entry into the BATL array of neighbor levels
*/
//void spectrum_get_all_neighbor_nodes(int node, int** neighbor_nodes, int** neighbor_levels);

/*!
\brief Obtain copies of the neighbor arrays for a node
\param[in]  node            Node
\param[out] neighbor_nodes  Array of neighbor nodes
\param[out] neighbor_levels Array of neighbor levels
*/
void spectrum_get_all_neighbor_copies(int node, int* neighbor_nodes, int* neighbor_levels);

/*!
\brief Return the coordinates of two opposite corners of the block
\param[in]  node     Node
\param[out] face_min Coordinates of the lower corner
\param[out] face_max Coordinates of the upper corner
*/
void spectrum_get_block_corners(int node, double* face_min, double* face_max);

/*!
\brief Obtain copies of the variables associated with a node
\param[in]  node      Node
\param[out] variables Array of variables
*/
void spectrum_get_block_data(int node, double* variables);

#ifdef __cplusplus
}
#endif

/*!
\brief Convert a node multi-index to a level milti-index
\param[in] node_idx Node index (0-3)
\return Level index: 0->0, 1->1, 2->1, 3->2
*/
inline MultiIndex NodeToLevel(const MultiIndex& node_idx)
{
   return (node_idx + 1) / 2;
};

/*!
\brief Convert a level multi-index to a node milti-index
\param[in] level_idx Level index (0-2)
\return Node index: 0->0, 1->2, 2->3
*/
inline MultiIndex LevelToNode(const MultiIndex& level_idx)
{
   return 2 * level_idx - level_idx / 2;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BlockBATL class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A three-dimensional cuboid region partitioned into zones
\author Vladimir Florinski
\author Juan G Alonso Guzman

This class is designed to facilitate the interaction with the BATL library. It corresponds to a single block in a simulation. Information stored include block dimensions, neighbor list, and data. 
*/
class BlockBATL : public BlockCartesian {

public:

//! Default constructor
   BlockBATL(void);

//! Destructor
   ~BlockBATL() override = default;

//! Obtain the extents of the block from BATL
   void LoadDimensions(double unit_length_block) override;

//! Request all neighbors from BATL
   void LoadNeighbors(void) override;

//! Load all variables into the block from BATL
   void LoadVariables(void) override;

//! Return the refinement level of a neighbor - multi-index version
   virtual int GetNeighborLevel(const MultiIndex& level_idx) const override;

//! Return the refinement level of a neighbor - three index version
   virtual int GetNeighborLevel(int i, int j, int k) const override;
};

////! Block type
//#if SERVER_TYPE == SERVER_BATL
//typedef BlockBATL BlockType;
//#endif

};

#endif
