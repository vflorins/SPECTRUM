/*!
\file neighbors.hh
\brief Defines the neighbor types
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include <string>

#ifndef SPECTRUM_NEIGHBORS_HH
#define SPECTRUM_NEIGHBORS_HH

namespace Spectrum {

//! Number of neighbor types
#define N_NBRTYPES 5

//! Neighbor types
enum NeighborType {
   GEONBR_TFACE,
   GEONBR_RFACE,
   GEONBR_TEDGE,
   GEONBR_REDGE,
   GEONBR_VERTX
};

//! Text description of neighbor types
constexpr std::string neighbor_text[] = {"t-face", "r-face", "t-edge", "r-edge", "vertex"};

//! Number of participating slabs per site (= number of exchange sites per slab)
constexpr int n_slab_parts[N_NBRTYPES] = {2, 1, 2, 1, 2};

};

#endif