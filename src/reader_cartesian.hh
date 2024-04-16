/*!
\file reader_cartesian.hh
\brief Defines a global structure of data reader for a uniform Cartesian grid
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_READER_CARTESIAN_HH
#define SPECTRUM_READER_CARTESIAN_HH

#include "common/vectors.hh"

namespace Spectrum {

//! Size of the block in each dimension
const MultiIndex block_size_cartesian(4, 4, 4);

//! Number of neighbors per dimension
const int max_neighbors_per_dim_cartesian = 3;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ReaderCartesian structure declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

struct ReaderCartesian {

//! Number of blocks per dimension
   MultiIndex Nblocks;

//! Number of zones per dimension
   MultiIndex domain_size;

//! Minimum domain coordinates
   GeoVector domain_min;

//! Maximum domain coordinates
   GeoVector domain_max;

//! Number of variables per zone
   int Nvar;

//! Number of variables per block
   int data_block_size;

//! Length of block in Cartesian units
   GeoVector block_length;

//! Length of zone in Cartesian units
   GeoVector zone_length;

//! Array with ALL Cartesian variables organized by block
   double* variables_by_block = nullptr;

//! Array with ALL Cartesian variables organized by dimension
   double* variables_by_dimen = nullptr;

};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ReaderCartesian structure global functions declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Read header containing Cartesian file metadata
void ReadCartesianHeader(const char* filename, int fname_len, int verbose);

//! Read file containing Cartesian data
void ReadCartesianData(const char* filename, int fname_len, int read_header, int verbose);

//! De-allocate global Cartesian structure arrays
void ReadCartesianClean(void);

//! Get domain coordinate limits
void ReadCartesianGetDomain(double* domain_min_out, double* domain_max_out);

//! Get variables from position
void ReadCartesianGetBlockData(const double* pos, double* vars, int* found);

//! Get node ID from position
void ReadCartesianGetNode(const double* pos, int* node_id);

//! Get block corners from node ID
void ReadCartesianGetBlockCorners(int node, double* face_min, double* face_max);

//! Get neighbor ID's from node ID
void ReadCartesianGetNodeNeighbors(int node, int* neighbor_nodes, int* neighbor_levels);

//! Get block variables from node ID
void ReadCartesianGetBlockData(int node, double* block_vars);

};

#endif
