/*!
\file grid_block.cc
\brief Implements the grid block class
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include <cstring>
#include <utility>

#include "common/physics.hh"
#include "geodesic/grid_block.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// GridBlock public methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/17/2024
*/
template <int verts_per_face>
GridBlock<verts_per_face>::GridBlock(void)
                         : SphericalSlab(),
                           GeodesicSector<verts_per_face>()
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Default constructing a GridBlock\n";
#endif
#endif
};

/*!
\author Vladimir Florinski
\date01/08/2025
\param[in] other Object to initialize from
\note The copy constructor for "SphericalSlab" calls its "SetDimensions()" method
\note The copy constructor for "GeodesicSector" calls its "SetDimensions()" method
*/
template <int verts_per_face>
GridBlock<verts_per_face>::GridBlock(const GridBlock& other)
                         : SphericalSlab(static_cast<const SphericalSlab&>(other)),
                           GeodesicSector<verts_per_face>(static_cast<const GeodesicSector<verts_per_face>&>(other))
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Copy constructing a GridBlock\n";
#endif
#endif

   if (other.side_length == -1) return;
   SetDimensions(other.side_length, other.ghost_width, other.n_shells, other.ghost_height, true);
   SetIndex(other.block_index);

   if (!other.mesh_associated) return;
   AssociateMesh(other.xi_in[0], other.xi_in[n_shells_withghost], other.corner_type, other.border_type, other.block_vert_cart, other.dist_map, true);
};

/*!
\author Vladimir Florinski
\date 01/08/2025
\param[in] other Object to move into this
\note The move constructor for "SphericalSlab" calls its "SetDimensions()" method
*/
template <int verts_per_face>
GridBlock<verts_per_face>::GridBlock(GridBlock&& other) noexcept
                         : SphericalSlab(std::move(static_cast<SphericalSlab&&>(other))),
                           GeodesicSector<verts_per_face>(std::move(static_cast<GeodesicSector<verts_per_face>&&>(other)))
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Move constructing a GridBlock (moving the content)\n";
#endif
#endif

   if (other.side_length == -1) {
      PrintMessage(__FILE__, __LINE__, "Move constructor called, but the dimension of the moved object was not set", true);
      return;
   };

// Copy the radial grid properties
   Rmin = other.Rmin;
   Rmax = other.Rmax;
   dxi = other.dxi;
   
// Move the radial grid
   xi_in = other.xi_in;
   r_in = other.r_in;
   other.xi_in = nullptr;
   other.r_in = nullptr;

   dist_map = std::move(other.dist_map);

// Copy the geodesic mesh properties
   SetIndex(other.block_index);
   mesh_associated = other.mesh_associated;
   n_singular = other.n_singular;
   std::memcpy(border_type, other.border_type, 2 * sizeof(bool));
   std::memcpy(corner_type, other.corner_type, verts_per_face * sizeof(bool));

// Move the duplicate element arrays
   std::memcpy(dup_vert, other.dup_vert, verts_per_face * sizeof(int**));
   other.dup_vert[0] = nullptr;

// Move the vertex coordinates
   block_vert_cart = other.block_vert_cart;
   other.block_vert_cart = nullptr;
   
#ifdef USE_SILO

   n_verts_silo = other.n_verts_silo;
   n_faces_silo = other.n_faces_silo;

// Move the SILO index arrays
   vert_to_silo = other.vert_to_silo;
   face_to_silo = other.face_to_silo;
   silo_to_vert = other.silo_to_vert;
   silo_to_face = other.silo_to_face;
   other.vert_to_silo = nullptr;
   other.face_to_silo = nullptr;
   other.silo_to_vert = nullptr;
   other.silo_to_face = nullptr;

#endif
};

/*!
\author Vladimir Florinski
\date 01/08/2025
\param[in] other Object to move into this
\note The move constructor for "SphericalSlab" calls its "SetDimensions()" method
*/
template <int verts_per_face>
GridBlock<verts_per_face>::GridBlock(GridBlock&& other)
                         : SphericalSlab(std::move(static_cast<SphericalSlab&&>(other))),
                           GeodesicSector<verts_per_face>(std::move(static_cast<GeodesicSector<verts_per_face>&&>(other)))
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Move constructing a GridBlock\n";
#endif
#endif

   Setup();
   if (other.side_length == -1) return;

// Copy the radial grid properties
   Rmin = other.Rmin;
   Rmax = other.Rmax;
   LogRmax_Rmin = other.LogRmax_Rmin;
   drp_ratio = other.drp_ratio;
   dxi = other.dxi;
   
// Move the radial grid
   xi_in = other.xi_in;
   r_in = other.r_in;
   r2_in = other.r2_in;
   r3_in = other.r3_in;
   rp_in = other.rp_in;
   dr = other.dr;
   r_mp = other.r_mp;
   drp = other.drp;
   other.xi_in = nullptr;
   other.r_in = nullptr;
   other.r2_in = nullptr;
   other.r3_in = nullptr;
   other.rp_in = nullptr;
   other.dr = nullptr;
   other.r_mp = nullptr;
   other.drp = nullptr;

   dist_map = std::move(other.dist_map);

// Copy the geodesic mesh properties
   block_index = other.block_index;
   mesh_associated = other.mesh_associated;
   n_singular = other.n_singular;
   std::memcpy(border_type, other.border_type, 2 * sizeof(bool));
   std::memcpy(corner_type, other.corner_type, verts_per_face * sizeof(bool));
   other.block_index = -1;

// Move the duplicate element arrays
   std::memcpy(dup_vert, other.dup_vert, verts_per_face * sizeof(int**));
   std::memcpy(dup_edge, other.dup_edge, verts_per_face * sizeof(int**));
   missing_faces = other.missing_faces;
   other.dup_vert[0] = nullptr;
   other.dup_edge[0] = nullptr;
   other.missing_faces = nullptr;

// Move the vertex coordinates
   block_vert_cart = other.block_vert_cart;
   other.block_vert_cart = nullptr;
   
#ifdef USE_SILO

// Move the SILO index arrays
   vert_to_silo = other.vert_to_silo;
   face_to_silo = other.face_to_silo;
   silo_to_vert = other.silo_to_vert;
   silo_to_face = other.silo_to_face;
   other.vert_to_silo = nullptr;
   other.face_to_silo = nullptr;
   other.silo_to_vert = nullptr;
   other.silo_to_face = nullptr;

#endif
};

/*!
\author Vladimir Florinski
\date 05/08/2024
\param[in] width  Length of the side, without ghost cells
\param[in] wgohst Width of the ghost cell layer outside the sector
\param[in] height Hight of the block, without ghost cells
\param[in] hghost Number of ghost shells outside the slab
*/
template <int verts_per_face>
GridBlock<verts_per_face>::GridBlock(int width, int wghost, int height, int hghost)
                                              : SphericalSlab(height, hghost),
                                                GeodesicSector<verts_per_face>(width, wghost)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Argument constructing a GridBlock\n";
#endif
#endif

   SetDimensions(width, wghost, height, hghost, true);
};

/*!
\author Vladimir Florinski
\date 05/16/2018
*/
template <int verts_per_face>
GridBlock<verts_per_face>::~GridBlock(void)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Destructing a GridBlock\n";
#endif
#endif

   FreeStorage();
};

/*!
\author Vladimir Florinski
\date 05/08/2024
\param[in] width     Length of the side, without ghost cells
\param[in] wghost    Width of the ghost cell layer outside the sector
\param[in] height    Hight of the block, without ghost cells
\param[in] hghost    Number of ghost shells outside the slab
\param[in] construct Set to true when called from a constructor
*/
template <int verts_per_face>
void GridBlock<verts_per_face>::SetDimensions(int width, int wghost, int height, int hghost, bool construct)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Setting dimensions " << width << " by " << height << " for a GridBlock\n";
#endif
#endif

// Call base methods
   if (!construct) {
      GeodesicSector<verts_per_face>::SetDimensions(width, wghost, false);
      SphericalSlab::SetDimensions(height, hghost, false);
   };

// Free up storage (not that of the base class) because this could be a repeat call
   FreeStorage();

// Allocate duplicate element lists (used for singular corners). We use contiguous storage so that the corners are stored in sequence and can be addressed as a single array.
   dup_vert[0] = Create2D<int>(verts_per_face * (ghost_width + 1), 2);
   for (auto corner = 1; corner < verts_per_face; corner++) dup_vert[corner] = dup_vert[corner - 1] + ghost_width + 1;

// Allocate coordinates
   block_vert_cart = new GeoVector[n_verts_withghost];

#ifdef USE_SILO

// Allocate SILO index arrays. The sizes of inverse SILO maps are not known at this time.
   vert_to_silo = new int[n_verts_withghost];
   face_to_silo = new int[n_faces_withghost];

#endif

// Allocate radial grid
   xi_in = new double[n_shells_withghost + 1];
   r_in = new double[n_shells_withghost + 1];
};

/*!
\author Vladimir Florinski
\date 02/17/2020
*/
template <int verts_per_face>
void GridBlock<verts_per_face>::FreeStorage()
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Freeing storage for a GridBlock\n";
#endif
#endif

// Free up radial arrays
   delete[] xi_in;
   delete[] r_in;

#ifdef USE_SILO

// Free up SILO index arrays
   delete[] vert_to_silo;
   delete[] face_to_silo;
   delete[] silo_to_vert;
   delete[] silo_to_face;

#endif

// Free up coordinates
   delete[] block_vert_cart;

// Free up duplicate element lists
   Delete2D(dup_vert[0]);
};

/*!
\author Vladimir Florinski
\date 01/09/2025
\param[in] ximin       Smallest reference distance of the block (without ghost)
\param[in] ximax       Largest reference distance of the block (without ghost)
\param[in] corners     Corner type, true for singular corners
\param[in] borders     Radial boundary type, true for external
\param[in] vcart       Vertex coordinate array in TAS/QAS
\param[in] dist_map_in Radial map function 
\param[in] construct   Set to true when called from a constructor
*/
template <int verts_per_face>
void GridBlock<verts_per_face>::AssociateMesh(double ximin, double ximax, const bool* corners, const bool* borders,
                                              const GeoVector* vcart, std::shared_ptr<DistanceBase> dist_map_in, bool construct)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Associating mesh for a GridBlock\n";
#endif
#endif

// The dimensions were not set, cannot proceed.
   if (side_length == -1) {
      PrintError(__FILE__, __LINE__, "Cannot associate mesh: block is not initialized", true);
      return;
   };

// Copy border information and correct connectivity at singular corners
   memcpy(corner_type, corners, verts_per_face * sizeof(bool));
   n_singular = 0;
   for (auto iv = 0; iv < verts_per_face; iv++) {
      if (corners[iv]) n_singular++;
   };
   memcpy(border_type, borders, 2 * sizeof(bool));
   FixSingularCorners();

#ifdef USE_SILO
// Generate the TAS/QAS to SILO mappings
   GenerateSiloIndexing();
#endif

   dist_map = dist_map_in;
   dxi = (ximax - ximin) / n_shells;

// Compute the interface coordinates using the supplied distance map
   for (auto k = 0; k <= n_shells_withghost; k++) {
      xi_in[k] = ximin + (k - ghost_height) * dxi;
      r_in[k] = dist_map->GetPhysical(xi_in[k]);
   };

// Set up the radial mesh in EC. The domain extents are retrieved from the radial map (they correspond to xi=0 and xi=1).
   Rmin = dist_map->GetPhysical(0.0);
   Rmax = dist_map->GetPhysical(1.0);

// Copy vertex coordinates
   memcpy(block_vert_cart, vcart, n_verts_withghost * sizeof(GeoVector));
   mesh_associated = true;
};

/*!
\author Vladimir Florinski
\date 03/14/2025
\param[in] index New block index
*/
template <int verts_per_face>
void GridBlock<verts_per_face>::SetIndex(int index)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Setting index " << index << " for a GridBlock\n";
#endif
#endif

   block_index = index;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// GridBlock protected methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/08/2024
*/
template <>
constexpr void GridBlock<3>::StaticSetup(void)
{
/*
            2          2---------1          1          1---------0          0          0---------2
           / \          \       /          / \          \       /          / \          \       /
          /   \          \     /          /   \          \     /          /   \          \     /
         /     \          \   /          /     \          \   /          /     \          \   /
        /       \          \ /          /       \          \ /          /       \          \ /
       0---------1          0          2---------0          2          1---------2          1
          rot=0           rot=1           rot=2           rot=3           rot=4           rot=5
*/

// Rotation of the sub-blocks located _diagonally_ from each corner vertex
   corner_rotation[0] = 3; corner_rotation[1] = 5; corner_rotation[2] = 1;

// The storage pattern is i(ir, jr), j(ir, jr). Use the table to calculate global i and j.
//       |     i=      |     j=
//  ---------------------------------
//  ir*  | [rot][0][0] | [rot][1][0]
//   +
//  jr*  | [rot][0][1] | [rot][1][1]
//
   rotated_verts[0][0][0] =  1; rotated_verts[0][0][1] =  0; rotated_verts[0][1][0] =  0; rotated_verts[0][1][1] =  1;
   rotated_verts[1][0][0] =  1; rotated_verts[1][0][1] = -1; rotated_verts[1][1][0] =  1; rotated_verts[1][1][1] =  0;
   rotated_verts[2][0][0] =  0; rotated_verts[2][0][1] = -1; rotated_verts[2][1][0] =  1; rotated_verts[2][1][1] = -1;
   rotated_verts[3][0][0] = -1; rotated_verts[3][0][1] =  0; rotated_verts[3][1][0] =  0; rotated_verts[3][1][1] = -1;
   rotated_verts[4][0][0] = -1; rotated_verts[4][0][1] =  1; rotated_verts[4][1][0] = -1; rotated_verts[4][1][1] =  0;
   rotated_verts[5][0][0] =  0; rotated_verts[5][0][1] =  1; rotated_verts[5][1][0] = -1; rotated_verts[5][1][1] =  1;

   rotated_shift[0][0][0] =  0; rotated_shift[0][0][1] =  0;
   rotated_shift[0][1][0] =  0; rotated_shift[0][1][1] =  0;
   rotated_shift[0][2][0] =  0; rotated_shift[0][2][1] =  0;
   rotated_shift[1][0][0] =  0; rotated_shift[1][0][1] =  0;
   rotated_shift[1][1][0] = -1; rotated_shift[1][1][1] =  0;
   rotated_shift[1][2][0] =  0; rotated_shift[1][2][1] =  0;
   rotated_shift[2][0][0] =  0; rotated_shift[2][0][1] =  0;
   rotated_shift[2][1][0] = -1; rotated_shift[2][1][1] = -1;
   rotated_shift[2][2][0] = -1; rotated_shift[2][2][1] =  0;
   rotated_shift[3][0][0] = -1; rotated_shift[3][0][1] =  0;
   rotated_shift[3][1][0] =  0; rotated_shift[3][1][1] = -1;
   rotated_shift[3][2][0] = -1; rotated_shift[3][2][1] = -1;
   rotated_shift[4][0][0] = -1; rotated_shift[4][0][1] = -1;
   rotated_shift[4][1][0] =  0; rotated_shift[4][1][1] =  0;
   rotated_shift[4][2][0] =  0; rotated_shift[4][2][1] = -1;
   rotated_shift[5][0][0] =  0; rotated_shift[5][0][1] = -1;
   rotated_shift[5][1][0] =  0; rotated_shift[5][1][1] =  0;
   rotated_shift[5][2][0] =  0; rotated_shift[5][2][1] =  0;

// The storage pattern is i(ir, jr_even, jr_odd), j(ir, jr_even, jr_odd). Use the table to calculate global i and j.
//       |              i=              |              j=
//  -------------------------------------------------------------------
//  ir*  | [rot][0][0]                  | [rot][1][0]
//   +
//  jr*  | [rot][0][1+jr%square_fill]   | [rot][1][1+jr%square_fill]

   rotated_faces[0][0][0] =  1; rotated_faces[0][0][1] =  0; rotated_faces[0][0][2] =  0;
   rotated_faces[0][1][0] =  0; rotated_faces[0][1][1] =  1; rotated_faces[0][1][2] =  1;
   rotated_faces[1][0][0] =  1; rotated_faces[1][0][1] = -1; rotated_faces[1][0][2] =  0;
   rotated_faces[1][1][0] =  2; rotated_faces[1][1][1] = -1; rotated_faces[1][1][2] =  1;
   rotated_faces[2][0][0] =  0; rotated_faces[2][0][1] =  0; rotated_faces[2][0][2] = -1;
   rotated_faces[2][1][0] =  2; rotated_faces[2][1][1] = -1; rotated_faces[2][1][2] = -1;
   rotated_faces[3][0][0] = -1; rotated_faces[3][0][1] =  0; rotated_faces[3][0][2] =  0;
   rotated_faces[3][1][0] =  0; rotated_faces[3][1][1] = -1; rotated_faces[3][1][2] = -1;
   rotated_faces[4][0][0] = -1; rotated_faces[4][0][1] =  1; rotated_faces[4][0][2] =  0;
   rotated_faces[4][1][0] = -2; rotated_faces[4][1][1] =  1; rotated_faces[4][1][2] = -1;
   rotated_faces[5][0][0] =  0; rotated_faces[5][0][1] =  0; rotated_faces[5][0][2] =  1;
   rotated_faces[5][1][0] = -2; rotated_faces[5][1][1] =  1; rotated_faces[5][1][2] =  1;

/*
                   5
                 . . .
               .   .   .
             .     .     .
            1-------------2
            |      4      |
            |     . .     |
            |   .     .   |
            | .         . |
            0-------------3

*/
#ifdef USE_SILO

   silo_zonetype = DB_ZONETYPE_PRISM;
   zv_silo[0][0] = 0; zv_silo[0][1] = 0; zv_silo[1][0] = 1; zv_silo[1][1] = 0; zv_silo[2][0] = 1; zv_silo[2][1] = 1;
   zv_silo[3][0] = 0; zv_silo[3][1] = 1; zv_silo[4][0] = 0; zv_silo[4][1] = 2; zv_silo[5][0] = 1; zv_silo[5][1] = 2;

#endif

};

/*!
\author Vladimir Florinski
\date 05/08/2024
*/
template <>
constexpr void GridBlock<4>::StaticSetup(void)
{
/*
       3---------2     2---------1     1---------0     0---------3
       |         |     |         |     |         |     |         |
       |         |     |         |     |         |     |         |
       |         |     |         |     |         |     |         |
       0---------1     3---------0     2---------3     1---------2
          rot=0           rot=1           rot=2           rot=3
*/

// Rotation of the sub-blocks located _diagonally_ from each corner vertex
   corner_rotation[0] = 2; corner_rotation[1] = 3; corner_rotation[2] = 0; corner_rotation[3] = 1;

// The storage pattern is i(ir, jr), j(ir, jr). Use the table to calculate global i and j.
//       |     i=      |     j=
//  ---------------------------------
//  ir*  | [rot][0][0] | [rot][1][0]
//   +
//  jr*  | [rot][0][1] | [rot][1][1]
//
   rotated_verts[0][0][0] =  1; rotated_verts[0][0][1] =  0; rotated_verts[0][1][0] =  0; rotated_verts[0][1][1] =  1;
   rotated_verts[1][0][0] =  0; rotated_verts[1][0][1] = -1; rotated_verts[1][1][0] =  1; rotated_verts[1][1][1] =  0;
   rotated_verts[2][0][0] = -1; rotated_verts[2][0][1] =  0; rotated_verts[2][1][0] =  0; rotated_verts[2][1][1] = -1;
   rotated_verts[3][0][0] =  0; rotated_verts[3][0][1] =  1; rotated_verts[3][1][0] = -1; rotated_verts[3][1][1] =  0;

   rotated_shift[0][0][0] =  0; rotated_shift[0][0][1] =  0; rotated_shift[0][1][0] =  0; rotated_shift[0][1][1] =  0;
   rotated_shift[1][0][0] =  0; rotated_shift[1][0][1] =  0; rotated_shift[1][1][0] = -1; rotated_shift[1][1][1] =  0;
   rotated_shift[2][0][0] = -1; rotated_shift[2][0][1] =  0; rotated_shift[2][1][0] =  0; rotated_shift[2][1][1] = -1;
   rotated_shift[3][0][0] =  0; rotated_shift[3][0][1] = -1; rotated_shift[3][1][0] =  0; rotated_shift[3][1][1] =  0;

// The storage pattern is i(ir, jr_even, jr_odd), j(ir, jr_even, jr_odd). Use the table to calculate global i and j.
//       |              i=              |              j=
//  -------------------------------------------------------------------
//  ir*  | [rot][0][0]                  | [rot][1][0]
//   +
//  jr*  | [rot][0][1+jr%square_fill]   | [rot][1][1+jr%square_fill]

   rotated_faces[0][0][0] =  1; rotated_faces[0][0][1] =  0; rotated_faces[0][1][0] =  0; rotated_faces[0][1][1] =  1;
   rotated_faces[1][0][0] =  0; rotated_faces[1][0][1] = -1; rotated_faces[1][1][0] =  1; rotated_faces[1][1][1] =  0;
   rotated_faces[2][0][0] = -1; rotated_faces[2][0][1] =  0; rotated_faces[2][1][0] =  0; rotated_faces[2][1][1] = -1;
   rotated_faces[3][0][0] =  0; rotated_faces[3][0][1] =  1; rotated_faces[3][1][0] = -1; rotated_faces[3][1][1] =  0;

/*
                 7-------------6
               . .           . |
             .   .         .   |
            4-------------5    |
            |    .        |    |
            |    3 . . . .|. . 2
            |  .          |  .
            |.            |.
            0-------------1
*/

#ifdef USE_SILO

   silo_zonetype = DB_ZONETYPE_HEX;
   zv_silo[0][0] = 0; zv_silo[0][1] = 0; zv_silo[1][0] = 0; zv_silo[1][1] = 1;
   zv_silo[2][0] = 0; zv_silo[2][1] = 2; zv_silo[3][0] = 0; zv_silo[3][1] = 3;
   zv_silo[4][0] = 1; zv_silo[4][1] = 0; zv_silo[5][0] = 1; zv_silo[5][1] = 1;
   zv_silo[6][0] = 1; zv_silo[6][1] = 2; zv_silo[7][0] = 1; zv_silo[7][1] = 3;

#endif

};

/*!
\author Vladimir Florinski
\date 05/08/2024
*/
template <int verts_per_face>
void GridBlock<verts_per_face>::FixSingularCorners(void)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Fixing singular corners for a GridBlock\n";
#endif
#endif

   int i, j, i1, j1, i2, j2, etype, i_origin, j_origin, i_origin1, j_origin1, i_origin2, j_origin2;
   int iv, iv1, iv2, iv3, iv4, ie, it, it1, it2;
   int vert, vert1, vert2, edge, edge1, edge2, face, face1, face2, rot, vert_origin, face_origin, bl_side, bl_corner, n_mf;
   bool dup_found;
   std::pair base_vert = std::make_pair(0, 0);
   int** dup_edge[verts_per_face];

   dup_edge[0] = Create2D<int>(verts_per_face * ghost_width, 2);
   for (auto corner = 1; corner < verts_per_face; corner++) dup_edge[corner] = dup_edge[corner - 1] + ghost_width;

#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Fixing singular corners for a GridBlock\n";
#endif
#endif

// Reset the "dup" arrays
   for (auto corner = 0; corner < verts_per_face; corner++) {
      for (auto i_rot = 0; i_rot <= ghost_width; i_rot++) {
         dup_vert[corner][i_rot][0] = -1;
         dup_vert[corner][i_rot][1] = -1;
         if (i_rot != ghost_width) {
            dup_edge[corner][i_rot][0] = -1;
            dup_edge[corner][i_rot][1] = -1;
         };
      };
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// First pass: eliminate all interior elements of the corner blocks and build the list of duplicate vertices and edges along the cut lines.
//----------------------------------------------------------------------------------------------------------------------------------------------------

   n_mf = 0;
   for (auto corner = 0; corner < verts_per_face; corner++) {
      if (!corner_type[corner]) continue;

// Find the corner vertex and determine the amount of rotation. Calculate the starting vertex.
      CornerCoords(corner, i_origin, j_origin);
      vert_origin = vert_index_sector[i_origin][j_origin];
      rot = corner_rotation[corner];

// Reset VV, VE, and VF of vertices inside the cut and build the lists of duplicate vertices
      for (auto i_rot = 0; i_rot <= ghost_width; i_rot++) {
         for (auto j_rot = 0; j_rot <= MaxVertJ(ghost_width, i_rot); j_rot++) {

// Calculate global index of the vertex
            i = i_origin + rotated_verts[rot][0][0] * i_rot + rotated_verts[rot][0][1] * j_rot;
            j = j_origin + rotated_verts[rot][1][0] * i_rot + rotated_verts[rot][1][1] * j_rot;
            vert = vert_index_sector[i][j];

// Check if it is a duplicate vertex (lies on the cut line). In rotated coordinates, the cut line is side "0" and side "verts_per_face-1". The main block corner enters into both "dup_vert[corner][i_rot][0]" and "dup_vert[corner][j_rot][1]".
            bl_side = BoundaryVert(base_vert, ghost_width, i_rot, j_rot);
            bl_corner = CornerVert(base_vert, ghost_width, i_rot, j_rot);
            dup_found = false;
            if ((bl_side == 0) || (bl_corner == 0) || (bl_corner == 1)) {
               dup_vert[corner][i_rot][0] = vert;
               dup_found = true;
            };
            if ((bl_side == verts_per_face - 1) || (bl_corner == verts_per_face - 1) || (bl_corner == 0)) {
               dup_vert[corner][j_rot][1] = vert;
               dup_found = true;
            }; 

// These are missing vertices
            if (!dup_found)  {
               for (iv = 0; iv < edges_per_vert; iv++) vv_local[vert][iv] = -1;
               for (ie = 0; ie < edges_per_vert; ie++) ve_local[vert][ie] = -1;
               for (it = 0; it < edges_per_vert; it++) vf_local[vert][it] = -1;
               RAISE_BITS(vert_mask[vert], GEOELM_NEXI);
            }
            else RAISE_BITS(vert_mask[vert], GEOELM_CUTL);
         };
      };

// Reset EV and EF of edges inside the cut (interior edges only).
      for (auto etype_rot = 0; etype_rot < cardinal_directions; etype_rot++) {
         for (auto i_rot = 0; i_rot <= ghost_width + edge_dimax[etype_rot]; i_rot++) {
            for (auto j_rot = 0; j_rot <= MaxVertJ(ghost_width, i_rot) + edge_djmax[etype_rot]; j_rot++) {
               
// Calculate global index of the edge
               etype = (etype_rot + 2 * cardinal_directions - rot) % cardinal_directions;
               i = i_origin + rotated_shift[rot][etype_rot][0] + rotated_verts[rot][0][0] * i_rot + rotated_verts[rot][0][1] * j_rot;
               j = j_origin + rotated_shift[rot][etype_rot][1] + rotated_verts[rot][1][0] * i_rot + rotated_verts[rot][1][1] * j_rot;
               edge = edge_index_sector[etype][i][j];

// Check if it is a duplicate edge (lies on the cut line). In rotated coordinates, the cut line is side "0" and side "verts_per_face-1".
               bl_side = BoundaryEdge(base_vert, ghost_width, etype_rot, i_rot, j_rot);
               dup_found = false;
               if (bl_side == 0) {
                  dup_edge[corner][i_rot][0] = edge;
                  dup_found = true;
               };
               if (bl_side == verts_per_face - 1) {
                  dup_edge[corner][j_rot][1] = edge;
                  dup_found = true;
               };

// These are missing edges
               if (!dup_found)  {
                  for (iv = 0; iv < 2; iv++) ev_local[edge][iv] = -1;
                  for (it = 0; it < 2; it++) ef_local[edge][it] = -1;
                  RAISE_BITS(edge_mask[edge], GEOELM_NEXI);
               }
               else RAISE_BITS(edge_mask[edge], GEOELM_CUTL);
            };
         };
      };

// Starting face inside the cut
      face_origin = vf_local[vert_origin][corner * square_fill];
      i_origin = face_index_i[face_origin];
      j_origin = face_index_j[face_origin];

// Starting replacement faces in the two boundary sub-blocks on each side of the cut
      it = (rot - 1 + edges_per_vert / 2) % edges_per_vert;
      face_origin = vf_local[vert_origin][it];
      i_origin1 = face_index_i[face_origin];
      j_origin1 = face_index_j[face_origin];
      it = (rot + 1 + edges_per_vert / 2) % edges_per_vert;
      face_origin = vf_local[vert_origin][it];
      i_origin2 = face_index_i[face_origin];
      j_origin2 = face_index_j[face_origin];

// Reset FV, FE, and FF of faces inside the cut.
      for (auto i_rot = 0; i_rot <= ghost_width - 1; i_rot++) {
         i = i_origin + rotated_faces[rot][0][0] * i_rot;
         j = j_origin + rotated_faces[rot][1][0] * i_rot;

         i1 = i_origin1 + rotated_faces[(rot + edges_per_vert - 1) % edges_per_vert][0][0] * i_rot;
         j1 = j_origin1 + rotated_faces[(rot + edges_per_vert - 1) % edges_per_vert][1][0] * i_rot;
         i2 = i_origin2 + rotated_faces[(rot + 1) % edges_per_vert][0][0] * i_rot;
         j2 = j_origin2 + rotated_faces[(rot + 1) % edges_per_vert][1][0] * i_rot;

         for (auto j_rot = 0; j_rot <= MaxFaceJ(ghost_width, i_rot); j_rot++) {
            face = face_index_sector[i][j];
            RAISE_BITS(face_mask[face], GEOELM_NEXI);

            for (iv = 0; iv < verts_per_face; iv++) fv_local[face][iv] = -1;
            for (ie = 0; ie < verts_per_face; ie++) fe_local[face][ie] = -1;
            for (it = 0; it < verts_per_face; it++) ff_local[face][it] = -1;
            n_mf++;

// The steps in i and j are different in TAS for odd and even j_rot, but the same in QAS.
            i += rotated_faces[rot][0][1 + j_rot % square_fill];
            j += rotated_faces[rot][1][1 + j_rot % square_fill];

            i1 += rotated_faces[(rot + edges_per_vert - 1) % edges_per_vert][0][1 + j_rot % square_fill];
            j1 += rotated_faces[(rot + edges_per_vert - 1) % edges_per_vert][1][1 + j_rot % square_fill];
            i2 += rotated_faces[(rot + 1) % edges_per_vert][0][1 + j_rot % square_fill];
            j2 += rotated_faces[(rot + 1) % edges_per_vert][1][1 + j_rot % square_fill];
         };
      };
   };

// Check the number of missing faces
   if (n_mf != n_singular * FaceCount(ghost_width)) {
      PrintError(__FILE__, __LINE__, "Error occurred while computing the missing faces", true);
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Second pass: stitch together the cut line, i.e., fix the connectivity for the mesh elements that are on or adjacent to the cut line. Note that there are two addresses for the elements on the cut line ("right" and "left") and the variables defined there should be synchronized to ensure consistency.
//----------------------------------------------------------------------------------------------------------------------------------------------------

   for (auto corner = 0; corner < verts_per_face; corner++) {
      if (!corner_type[corner]) continue;
      rot = corner_rotation[corner];

// Determine which vertex neighbors to replace. First compute "iv1", which is the first vertex neighbor of the right vetex to replace. Its complimentary on the left side is iv2=iv1-1.
      iv1 = (rot - corner_rotation[0] + edges_per_vert) % edges_per_vert;
      iv2 = (iv1 + edges_per_vert - 1) % edges_per_vert;

// Fix VV, VE, VF of vertices on the cut line. The limits on this loop are wider by 1 in each direction to acommodate VF.
      for (auto i_rot = 1; i_rot <= ghost_width; i_rot++) {
         vert1 = dup_vert[corner][i_rot][0];
         vert2 = dup_vert[corner][i_rot][1];

// For QAS 1 VV/VE and 2 VF neighbors are replaced and for TAS 2 VV/VE and 3 VF neighbors are replaced. The limits on this loop are wider by 1 in each direction to acommodate VF.
         for (auto iiv = -1; iiv < edges_per_vert / 2; iiv++) {

// Right vertex - index increases
            iv3 = (iv1 + iiv + edges_per_vert) % edges_per_vert;
            iv4 = (iv3 + 1 + edges_per_vert) % edges_per_vert;
            if (iiv != -1) {
               if (iiv != edges_per_vert / 2 - 1) {
                  vv_local[vert1][iv3] = vv_local[vert2][iv4];
                  ve_local[vert1][iv3] = ve_local[vert2][iv4];
               };
               vf_local[vert1][iv3] = vf_local[vert2][iv4];
            };

// Left vertex - index decreases
            iv3 = (iv2 - iiv + edges_per_vert) % edges_per_vert;
            iv4 = (iv3 - 1 + edges_per_vert) % edges_per_vert;
            if (iiv != edges_per_vert / 2 - 1) {
               if (iiv != -1) {
                  vv_local[vert2][iv3] = vv_local[vert1][iv4];
                  ve_local[vert2][iv3] = ve_local[vert1][iv4];
               };
               vf_local[vert2][iv3] = vf_local[vert1][iv4];
            };
         };
      };

// The corner has one less neighbor. We can eliminate either neighbor "iv1" or "iv2".
      vert1 = dup_vert[corner][0][0];
      vv_local[vert1][iv1] = -1;
      ve_local[vert1][iv1] = -1;
      vf_local[vert1][iv1] = -1;

// Determine which edge neighbors to replace. This depends on the edge orientation (inward or outward along the cut), so we test the first edges (right and left) in "dup_edge" to see where "vert1" is located in their EV lists.
      edge1 = dup_edge[corner][0][0];
      edge2 = dup_edge[corner][0][1];

// Right side: replace neighbor 0 if outward, neighbor 1 if inward. Left side: repalce neighbor 1 if outward, neighbor 2 if inward.
      it1 = (vert1 == ev_local[edge1][0] ? 1 : 0);
      it2 = (vert1 == ev_local[edge2][0] ? 0 : 1);

// Fix EF, and FF of faces adjacent to the cut line - one element only.
      for (auto i_rot = 0; i_rot <= ghost_width - 1; i_rot++) {
         edge1 = dup_edge[corner][i_rot][0];
         edge2 = dup_edge[corner][i_rot][1];
         face1 = ef_local[edge1][it1];
         face2 = ef_local[edge2][it2];

         ef_local[edge1][(it1 + 1) % 2] = face2;
         ef_local[edge2][(it2 + 1) % 2] = face1;
         ff_local[face1][InList(verts_per_face, fe_local[face1], edge1)] = face2;
         ff_local[face2][InList(verts_per_face, fe_local[face2], edge2)] = face1;
      };
   };

   Delete2D(dup_edge[0]);
};

#ifdef USE_SILO

/*!
\author Vladimir Florinski
\date 05/08/2024

Grid blocks are saved as individual mesh objects. A visualizer (VisIt) must be able to interpolate the zone based variables to the vertices to make a smoother looking rendering of the data. This requires one layer of ghost zones on each side of the block that abuts another block, and no ghost zones at physical boundaries. This function generates the face and vertex mappings to facilitate translation between the TAS/QAS and the SILO mesh indices.
*/
template <int verts_per_face>
void GridBlock<verts_per_face>::GenerateSiloIndexing(void)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Generating SILO indexing for a GridBlock\n";
#endif
#endif

   int iv, it, vert, vert_silo, face, face_silo, idx_dup, imax, jmax;
   bool include_ghost;
   std::pair<int, int> base_vert = std::make_pair(square_fill * (ghost_width - 1), ghost_width - 1);

// Generate the TAS/QAS to SILO vertex map that includes interior and one layer of ghost vertices
   for (vert = 0; vert < n_verts_withghost; vert++) vert_to_silo[vert] = -1;
   n_verts_silo = 0;
   imax = total_length - ghost_width + 1;

   for (auto i = square_fill * ghost_width - 1; i <= imax; i++) {
      jmax = MaxVertJ(base_vert, total_length - square_fill * ghost_width + 1, i);
      for (auto j = ghost_width - 1; j <= jmax; j++) {
         vert = vert_index_sector[i][j];

// Interior vertices always included. This includes the corners.
         if (IsInteriorVertOfSector(vert)) vert_to_silo[vert] = n_verts_silo++;

// If the vertex belongs to a face that has at least one vertex in the interior, it is included in the list.
         else {
            it = 0;
            include_ghost = false;
            while (!include_ghost && (it < edges_per_vert)) {
               face = vf_local[vert][it];
               if (face != -1) {
                  iv = 0;
                  while (iv < verts_per_face && !IsInteriorVertOfSector(fv_local[face][iv])) iv++;
                  if (iv != verts_per_face) include_ghost = true;
               };
               it++;
            };

// The vertex is part of the SILO ghost layer and eligible to be included in the list. If it is a duplicate vertex (lies on the cut line), it must be only included once.
            if (include_ghost) {
               if (BITS_RAISED(vert_mask[vert], GEOELM_CUTL)) {
                  idx_dup = InList(verts_per_face * (ghost_width + 1) * 2, dup_vert[0][0], vert);
                  if (is_odd(idx_dup)) continue;
               };
               vert_to_silo[vert] = n_verts_silo++;
            };
         };
      };
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Generate the TAS/QAS to SILO face map that includes interior and one layer of ghost faces
   for (face = 0; face < n_faces_withghost; face++) face_to_silo[face] = -1;
   n_faces_silo = 0;
   imax = total_length - ghost_width;

   for (auto i = square_fill * ghost_width - 1; i <= imax; i++) {
      jmax = MaxFaceJ(base_vert, total_length - square_fill * ghost_width + 1, i);
      for (auto j = square_fill * (ghost_width - 1); j <= jmax; j++) {
         face = face_index_sector[i][j];

// Interior faces always included
         if (IsInteriorFaceOfSector(face)) face_to_silo[face] = n_faces_silo++;

// If the face has at least one vertex in the interior, it is included in the list.
         else {
            iv = 0;
            while (iv < verts_per_face && !IsInteriorVertOfSector(fv_local[face][iv])) iv++;
            if (iv != verts_per_face) face_to_silo[face] = n_faces_silo++;
         };
      };
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Allocate and generate the inverse maps (SILO to TAS/QAS vertex and face)
   silo_to_vert = new int[n_verts_silo];
   silo_to_face = new int[n_faces_silo];

   vert_silo = 0;
   for (vert = 0; vert < n_verts_withghost; vert++) {
      if (vert_to_silo[vert] != -1) silo_to_vert[vert_silo++] = vert;
   };

   face_silo = 0;
   for (face = 0; face < n_faces_withghost; face++) {
      if (face_to_silo[face] != -1) silo_to_face[face_silo++] = face;
   };
};

/*!
\author Vladimir Florinski
\date 05/08/2024
\param[in] silofile   SILO database file object
\param[in] phys_units Use physical units for output
\return Error code from SILO file I/O (zero if no error)
*/
template <int verts_per_face>
int GridBlock<verts_per_face>::WriteSiloMesh(DBfile* silofile, bool phys_units) const
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 2
   std::cerr << "Writing mesh for block " << block_index << "\n";
#endif
#endif

   int k, k1, k2, vert, vert_silo, face, face_silo, idx1, idx2, idx_dup, ivz, n_nodes, n_zones, n_intzones, lznodelist, xyz;

// Total number of elements in the block. Notice that ghost cells are not added at the physical domain boundaries (innermost and outermost slabs) as required by SILO and VisIt.
   n_nodes = n_verts_silo * (n_shells + 1);
   n_zones = n_faces_silo * n_shells;
   if (!border_type[0]) {
      n_nodes += n_verts_silo;
      n_zones += n_faces_silo;
   };
   if (!border_type[1]) {
      n_nodes += n_verts_silo;
      n_zones += n_faces_silo;
   };
   n_intzones = n_faces * n_shells;

// Arrays of node (vertex) coordinates.
   double* coords[3];
   for (xyz = 0; xyz < 3; xyz++) coords[xyz] = new double[n_nodes];

// Calculate vertex Cartesian coordinates (with ghost vertices)
   idx1 = 0;
   k1 = (border_type[0] ? ghost_height : ghost_height - 1);
   k2 = (border_type[1] ? n_shells_withghost - ghost_height : n_shells_withghost - ghost_height + 1);
   for (k = k1; k <= k2; k++) {
      for (vert_silo = 0; vert_silo < n_verts_silo; vert_silo++) {
         vert = silo_to_vert[vert_silo];
         for (xyz = 0; xyz < 3; xyz++) {
            coords[xyz][idx1] = r_in[k] * block_vert_cart[vert][xyz] * (phys_units ? unit_length_fluid : 1.0);
         };
         idx1++;
      };
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// The array "znodelist" stores the list of vertices of each zone.
   lznodelist = 2 * verts_per_face * n_zones;
   int* znodelist = new int[lznodelist];

// Number and type of zone shapes (we have only one)
   int zshapecnt[1], zshapesize[1], zshapetype[1];
   zshapecnt[0] = n_zones;
   zshapesize[0] = 2 * verts_per_face;
   zshapetype[0] = silo_zonetype;

// The first index tracks the interior zones, while the second enumerates the ghost zones.
   idx1 = 0;
   idx2 = 2 * verts_per_face * n_intzones;
   k1 = (border_type[0] ? ghost_height : ghost_height - 1);
   k2 = (border_type[1] ? n_shells_withghost - ghost_height - 1 : n_shells_withghost - ghost_height);

// Build vertex lists for each SILO zone
   for (k = k1; k <= k2; k++) {
      for (face_silo = 0; face_silo < n_faces_silo; face_silo++) {      
         face = silo_to_face[face_silo];

// Add to the interior zone list.
         if (IsInteriorShellOfSlab(k) && IsInteriorFaceOfSector(face)) {
            for (ivz = 0; ivz < 2 * verts_per_face; ivz++) {
               vert = fv_local[face][zv_silo[ivz][1]];
               znodelist[idx1++] = (k - k1 + zv_silo[ivz][0]) * n_verts_silo + vert_to_silo[vert];
            };
         }

// Add to the ghost zone list.
         else {
            for (ivz = 0; ivz < 2 * verts_per_face; ivz++) {
               vert = fv_local[face][zv_silo[ivz][1]];

// Check if the vertex is a duplicate (odd idx_dup), in which case it must be replaced with the primary (idx_dup-1).
               if (vert_mask[vert] & GEOELM_CUTL) {
                  idx_dup = InList(verts_per_face * (ghost_width + 1) * 2, dup_vert[0][0], vert);
                  if (is_odd(idx_dup)) vert = dup_vert[0][0][idx_dup - 1];
               };
               znodelist[idx2++] = (k - k1 + zv_silo[ivz][0]) * n_verts_silo + vert_to_silo[vert];
            };
         };
      };
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Unique indices for the mesh objects
   char blck_num[blck_length + 1];
   sprintf(blck_num, blck_format.c_str(), block_index);
   std::string facelistname = fl_base + "_" + blck_num;
   std::string zonelistname = zl_base + "_" + blck_num;
   std::string meshname     = um_base + "_" + blck_num;

// Build the list of external faces of the block if it has a physical boundary.
   bool ext_faces = border_type[0] || border_type[1];
   int err = 0;
   DBfacelist* fl;
   if (ext_faces) {
      fl = DBCalcExternalFacelist2(znodelist, n_nodes, 0, n_zones - n_intzones, 0, zshapetype, zshapesize, zshapecnt, 1, NULL, 0);
      err = DBPutFacelist(silofile, facelistname.c_str(), fl->nfaces, 3, fl->nodelist,
                          fl->lnodelist, 0, fl->zoneno, fl->shapesize, fl->shapecnt, fl->nshapes, NULL, NULL, 0);
   };
   err |= DBPutZonelist2(silofile, zonelistname.c_str(), n_zones, 3, znodelist,
                         lznodelist, 0, 0, n_zones - n_intzones, zshapetype, zshapesize, zshapecnt, 1, NULL);

// Write the unstructured mesh object to the SILO database
   err |= DBPutUcdmesh(silofile, meshname.c_str(), 3, NULL, coords, n_nodes, n_zones, zonelistname.c_str(),
                       ext_faces ? facelistname.c_str() : NULL, DB_DOUBLE, NULL);
   if (ext_faces) DBFreeFacelist(fl);

// Clean up
   delete[] znodelist;
   for (xyz = 0; xyz < 3; xyz++) delete[] coords[xyz];
   return err;
};

/*!
\author Vladimir Florinski
\date 05/11/2018
\param[in] silofile   SILO database file object
\param[in] phys_units Use physical units for output
\return Error code from file I/O (zero if no error)
*/
template <int verts_per_face>
int GridBlock<verts_per_face>::WriteSilo(DBfile* silofile, bool phys_units) const
{
// A GridBlock object can only write the mesh information. Derived classes will extend this functionality.
   return WriteSiloMesh(silofile, phys_units);
};

#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// GridBlock debug/testing methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

#ifdef GEO_DEBUG

/*!
\author Vladimir Florinski
\date 07/23/2019
*/
template <int verts_per_face>
void GridBlock<verts_per_face>::PrintStats(void) const
{
   std::cerr << std::endl;
   std::cerr << "--------------------------------------------------------------------------------\n";
   std::cerr << "Printing stats for block " << block_index << std::endl;
   std::cerr << "--------------------------------------------------------------------------------\n";
   std::cerr << "Height:          " << std::setw(9) << n_shells_withghost << " | "
             << "Interior height: " << std::setw(9) << n_shells  << std::endl;
   std::cerr << "Width:           " << std::setw(9) << total_length << " | "
             << "Interior width:  " << std::setw(9) << side_length << std::endl;
   std::cerr << "Faces:           " << std::setw(9) << n_faces_withghost << " | "
             << "Interior faces:  " << std::setw(9) << n_faces << std::endl;
   std::cerr << "Zones:           " << std::setw(9) << n_shells_withghost * n_faces_withghost << " | "
             << "Interior zones:  " << std::setw(9) << n_shells * n_faces << std::endl;
   std::cerr << "--------------------------------------------------------------------------------\n";
   std::cerr << "Corners: ";
   for (auto corner = 0; corner < verts_per_face; corner++) {
      std::cerr << std::setw(5) << corner_type[corner];
   };
   std::cerr << std::endl;
   std::cerr << "Borders: ";
   for (auto border = 0; border < 2; border++) {
      std::cerr << std::setw(5) << border_type[border];
   };
   std::cerr << std::endl;
   std::cerr << "--------------------------------------------------------------------------------\n";
   std::setprecision(5);
   std::cerr << "Rmin:      " << std::setw(12) << r_in[ghost_height] << " | "
             << "Rmax:      " << std::setw(12) << r_in[n_shells_withghost - ghost_height] << std::endl;
   std::cerr << "--------------------------------------------------------------------------------\n";
   std::cerr << std::endl;
};

/*!
\author Vladimir Florinski
\date 05/08/2024
\param[in] k     Shell index
\param[in] face  Face index
\param[in] rot_z The first rotation angle about the z-axis
\param[in] rot_x The second rotation angle about the x-axis
*/
template <int verts_per_face>
void GridBlock<verts_per_face>::DrawZone(int k, int face, double rot_z, double rot_x) const
{
   int ie, iv, ipt;
   double ddist;

// Precompute trig functions of the rotation angles
   double snrz = sin(rot_z);
   double csrz = cos(rot_z);
   double snrx = sin(rot_x);
   double csrx = cos(rot_x);

// Number of segments used to draw each arc
   const int nint = 100;
   GeoVector v, v3, vv[verts_per_face];

// Obtain the coordinates of the vertices for this t-face
   for (iv = 0; iv < verts_per_face; iv++) vv[iv] = block_vert_cart[fv_local[face][iv]];

   std::ofstream zdfile;
   std::string zdname = "zone_drawing_" + std::to_string(block_index) + "_" + std::to_string(k) + "_" + std::to_string(face) + ".out";
   zdfile.open(zdname.c_str(), std::ofstream::out);

// Draw the top t-edges
   for (ie = 0; ie < verts_per_face; ie++) {
   
// Compute the angle between the points and divide it into equal segments
      ddist = acos(vv[ie] * vv[(ie + 1) % verts_per_face]) / (double)nint;
      v3 = (vv[ie] ^ vv[(ie + 1) % verts_per_face]) ^ vv[ie];
      v3.Normalize();

// Draw the points to form an edge
      for (ipt = 0; ipt <= nint; ipt++) {

// Interior point vectors are linear combinations of the two vertices
         v = cos(ipt * ddist) * vv[ie] + sin(ipt * ddist) * v3;

// Rotate the view
         v.Rotate(gv_nz, snrz, csrz);
         v.Rotate(gv_nx, snrx, csrx);

// Print the projection on the xz plane
         zdfile << std::setw(12) << std::setprecision(5) << v[0] * r_in[k + 1]
                << std::setw(12) << std::setprecision(5) << v[2] * r_in[k + 1]
                << std::endl;
      };
      zdfile << std::endl;
   };

// Draw the bottom t-edges
   for (ie = 0; ie < verts_per_face; ie++) {
   
// Compute the angle between the points and divide it into equal segments
      ddist = acos(vv[ie] * vv[(ie + 1) % verts_per_face]) / (double)nint;
      v3 = (vv[ie] ^ vv[(ie + 1) % verts_per_face]) ^ vv[ie];
      v3.Normalize();

// Draw the points to form an edge
      for (ipt = 0; ipt <= nint; ipt++) {

// Interior point vectors are linear combinations of the two vertices
         v = cos(ipt * ddist) * vv[ie] + sin(ipt * ddist) * v3;

// Rotate the view
         v.Rotate(gv_nz, snrz, csrz);
         v.Rotate(gv_nx, snrx, csrx);

// Print the projection on the xz plane
         zdfile << std::setw(12) << std::setprecision(5) << v[0] * r_in[k]
                << std::setw(12) << std::setprecision(5) << v[2] * r_in[k]
                << std::endl;
      };
      zdfile << std::endl;
   };

// Draw the side edges (two points per edge is enough for straight edges)
   for (iv = 0; iv < verts_per_face; iv++) {
      v = vv[iv];

// Rotate the view
      v.Rotate(gv_nz, snrz, csrz);
      v.Rotate(gv_nx, snrx, csrx);

// Print the projection on the xz plane
      zdfile << std::setw(12) << std::setprecision(5) << v[0] * r_in[k + 1]
             << std::setw(12) << std::setprecision(5) << v[2] * r_in[k + 1]
             << std::endl
             << std::setw(12) << std::setprecision(5) << v[0] * r_in[k]
             << std::setw(12) << std::setprecision(5) << v[2] * r_in[k]
             << std::endl << std::endl;
   };

   zdfile.close();
};

#endif

template class GridBlock<3>;
//template class GridBlock<4>;

};
