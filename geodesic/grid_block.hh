/*!
\file grid_block.hh
\brief Declares grid block class, the smallest computational unit in the mesh
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_GRID_BLOCK
#define SPECTRUM_GRID_BLOCK

#ifdef USE_SILO
#include <silo.h>
#endif
#include "common/vectors.hh"
#include "geometry/distance_map.hh"
#include "geodesic/geodesic_sector.hh"
#include "geodesic/spherical_slab.hh"

namespace Spectrum {

//! Mask for the cut line near singular corners (EV)
#define GEOELM_CUTL 0x0010

#ifdef USE_SILO

//! Pattern for facelist naming
const std::string fl_base = "geofl";

//! Pattern for zonelist naming
const std::string zl_base = "geozl";

//! Pattern for mesh naming
const std::string um_base = "geoum";

//! The length, in characters, of the block index
const int blck_length = 6;

//! The format for the block index
const std::string blck_format = "%0" + std::to_string(blck_length) + 'i';

#endif

/*!
\brief A class describing a sector/slab intersection
\author Vladimir Florinski

The grid block is the smallest computational unit that must reside on a single core. Zones are arranged in radial layers (shells), where each shell is a sector.
*/
template <int verts_per_face>
class GridBlock : public GeodesicSector<verts_per_face>, public SphericalSlab
{
protected:

   using PolygonalAddressing<verts_per_face>::edges_per_vert;
   using PolygonalAddressing<verts_per_face>::cardinal_directions;
   using PolygonalAddressing<verts_per_face>::square_fill;
   using PolygonalAddressing<verts_per_face>::edge_dimax;
   using PolygonalAddressing<verts_per_face>::edge_djmax;
   using PolygonalAddressing<verts_per_face>::FaceCount;
   using GeodesicSector<verts_per_face>::side_length;
   using GeodesicSector<verts_per_face>::total_length;
   using GeodesicSector<verts_per_face>::ghost_width;
   using GeodesicSector<verts_per_face>::n_verts;
   using GeodesicSector<verts_per_face>::n_faces;
   using GeodesicSector<verts_per_face>::n_verts_withghost;
   using GeodesicSector<verts_per_face>::n_faces_withghost;
   using GeodesicSector<verts_per_face>::vert_index_sector;
   using GeodesicSector<verts_per_face>::edge_index_sector;
   using GeodesicSector<verts_per_face>::face_index_sector;
   using GeodesicSector<verts_per_face>::vert_index_i;
   using GeodesicSector<verts_per_face>::vert_index_j;
   using GeodesicSector<verts_per_face>::face_index_i;
   using GeodesicSector<verts_per_face>::face_index_j;
   using GeodesicSector<verts_per_face>::vv_local;
   using GeodesicSector<verts_per_face>::ve_local;
   using GeodesicSector<verts_per_face>::vf_local;
   using GeodesicSector<verts_per_face>::ev_local;
   using GeodesicSector<verts_per_face>::ef_local;
   using GeodesicSector<verts_per_face>::fv_local;
   using GeodesicSector<verts_per_face>::fe_local;
   using GeodesicSector<verts_per_face>::ff_local;
   using GeodesicSector<verts_per_face>::vert_mask;
   using GeodesicSector<verts_per_face>::edge_mask;
   using GeodesicSector<verts_per_face>::face_mask;
   using GeodesicSector<verts_per_face>::MaxVertJ;
   using GeodesicSector<verts_per_face>::MaxFaceJ;
   using GeodesicSector<verts_per_face>::BoundaryVert;
   using GeodesicSector<verts_per_face>::BoundaryEdge;
   using GeodesicSector<verts_per_face>::CornerVert;

//! Rotations of the corner ghost blocks
   int corner_rotation[verts_per_face];

//! Translation between rotated and normal TAS/QAS for faces
   int rotated_verts[edges_per_vert][2][2];

//! Constant shifts for rotated edge indices
   int rotated_shift[edges_per_vert][cardinal_directions][2];

//! Translation between rotated and normal TAS/QAS for faces
   int rotated_faces[edges_per_vert][2][1 + square_fill];

//! Unique numerical index of the block
   int block_index;

//! Flag telling whether this block is in the innermost or outermost slab
   bool border_type[2];

//! Number of singular corners
   int n_singular;

//! Flag telling whether corners are singular
   bool corner_type[verts_per_face];

//! List of duplicate vertices at cut lines
   int** dup_vert[verts_per_face] = {nullptr};

//! List of duplicate edges at cut lines
   int** dup_edge[verts_per_face] = {nullptr};

//! Mapping of the existing faces into the missing block at singular corners
   int** missing_faces = nullptr;

//! Radial distance of the lower boundary of the entire domain
   double Rmin;

//! Radial distance of the upper boundary of the entire domain
   double Rmax;

// FIXME - check if this is needed
//! Logarithmic ratio of radial boundaries
   double LogRmax_Rmin;

//! Constant radio of EC shell width to EC radius of bottom of zone
   double drp_ratio;

//! LUC shell width
   double dxi;

//! Reference radial coordinates of shell interfaces
   double* xi_in = nullptr;

//! Radial coordinates of shell interfaces
   double* r_in = nullptr;

//! Square of the radial coordinates of shell interfaces
   double* r2_in = nullptr;

//! Cube of the radial coordinates of shell interfaces
   double* r3_in = nullptr;

//! Exponential radial coordinates of shell interfaces
   double* rp_in = nullptr;

//! Shell widths
   double* dr = nullptr;

//! Shell midpoints, \f$(r_1+r_2)/2\f$
   double* r_mp = nullptr;

//! Shell widths in exponential coordinates
   double* drp = nullptr;

//! Distance map object
   std::shared_ptr<DistanceBase> dist_map = nullptr;

//! Cartesian coordinates of vertices on the US
   GeoVector* block_vert_cart = nullptr;

#ifdef USE_SILO

//! Zone type
   int silo_zonetype;

//! Vertex numbering in a zone
   int zv_silo[2 * verts_per_face][2];

//! Number of SILO vertices
   int n_verts_silo;

//! Number of SILO t-faces
   int n_faces_silo;

//! TAS/QAS to SILO index vertex map
   int* vert_to_silo = nullptr;

//! TAS/QAS to SILO index face map
   int* face_to_silo = nullptr;

//! SILO to TAS/QAS index vertex map
   int* silo_to_vert = nullptr;

//! SILO to TAS/QAS index face map
   int* silo_to_face = nullptr;

//! Compute zone index translation between TAS/QAS and SILO
   void GenerateSiloIndexing(void);

//! Write the mesh description to a SILO database
   int WriteSiloMesh(DBfile* silofile, bool phys_units) const;

#endif

//! Return (i,j) coordinates of a corner vertex
   SPECTRUM_DEVICE_FUNC void CornerCoords(int corner, int& i, int& j) const;

//! Determine whether a vertex belongs to the sector's interior (including boundary) - vert version
   SPECTRUM_DEVICE_FUNC bool IsInteriorVert_Int(int vert) const;

//! Determine whether a face is in the sector's interior - face version
   SPECTRUM_DEVICE_FUNC bool IsInteriorFace_Int(int face) const;

//! Calculate the increments for each object type
   SPECTRUM_DEVICE_FUNC void Setup(void);

//! Correct the connectivity at the singular corners
   SPECTRUM_DEVICE_FUNC void FixSingularCorners(void);

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC GridBlock(void);

//! Constructor with arguments
   SPECTRUM_DEVICE_FUNC GridBlock(int width, int wghost, int height, int hghost);

//! Destructor
   SPECTRUM_DEVICE_FUNC ~GridBlock();

//! Allocate memory
   SPECTRUM_DEVICE_FUNC void SetDimensions(int width, int wghost, int height, int hghost, bool construct);

//! Free all dynamically allocated memory
   SPECTRUM_DEVICE_FUNC void FreeStorage(void);

//! Set up the dimensions and geometry of the mesh
   SPECTRUM_DEVICE_FUNC void AssociateMesh(int index, double ximin, double ximax, const bool* corners, const bool* borders, const GeoVector* vcart,
                                           std::shared_ptr<DistanceBase> dist_map_in);

#ifdef USE_SILO
//! Write the entire block to a SILO database
   int WriteSilo(DBfile* silofile, bool phys_units) const;
#endif

#ifdef GEO_DEBUG

//! Print the block information
   void PrintStats(void) const;

//! Draw a picture of a single zone
   void DrawZone(int k, int face, double rot_z, double rot_x) const;

#endif

};

/*!
\author Vladimir Florinski
\date 05/08/2024
\param[in]  corner Corner index
\param[out] i      First index of corner vertex
\param[out] j      Second index of corner vertex
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC inline void GridBlock<verts_per_face>::CornerCoords(int corner, int& i, int& j) const
{
   switch(corner) {

   case 0:
      i = square_fill * ghost_width;
      j = ghost_width;
      break;

   case 1:
      i = total_length - ghost_width;
      j = ghost_width;
      break;

   case 2:
      i = total_length - ghost_width;
      j = MaxVertJ(std::make_pair(square_fill * ghost_width, ghost_width), side_length, total_length - ghost_width);
      break;

   case 3:
      i = square_fill * ghost_width;
      j = MaxVertJ(std::make_pair(square_fill * ghost_width, ghost_width), side_length, total_length - ghost_width);
      break;
   };
};

/*!
\author Vladimir Florinski
\date 05/08/2024
\param[in] vert The vertex index
\return True if the vertex is interior to the sector
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC inline bool GridBlock<verts_per_face>::IsInteriorVert_Int(int vert) const
{
   if((vert < 0) || (vert >= n_verts_withghost)) return false;
   else return GeodesicSector<verts_per_face>::IsInteriorVert_Int(vert_index_i[vert], vert_index_j[vert]);
};

/*!
\author Vladimir Florinski
\date 05/08/2024
\param[in] face The face index
\return True if the face is interior to the sector
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC inline bool GridBlock<verts_per_face>::IsInteriorFace_Int(int face) const
{
   if((face < 0) || (face >= n_faces_withghost)) return false;
   else return GeodesicSector<verts_per_face>::IsInteriorFace_Int(face_index_i[face], face_index_j[face]);
};

/*!
\author Vladimir Florinski
\date 05/08/2024
*/
template <>
SPECTRUM_DEVICE_FUNC inline void GridBlock<3>::Setup(void)
{
//          2          2---------1          1          1---------0          0          0---------2
//         / \          \       /          / \          \       /          / \          \       / 
//        /   \          \     /          /   \          \     /          /   \          \     /  
//       /     \          \   /          /     \          \   /          /     \          \   /   
//      /       \          \ /          /       \          \ /          /       \          \ /    
//     0---------1          0          2---------0          2          1---------2          1     
//
//        rot=0           rot=1           rot=2           rot=3           rot=4           rot=5

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

//                 5
//               . . .     
//             .   .   .   
//           .     .     . 
//          1-------------2
//          |      4      |
//          |     . .     |
//          |   .     .   |
//          | .         . |
//          0-------------3

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
SPECTRUM_DEVICE_FUNC inline void GridBlock<4>::Setup(void)
{
//     3---------2     2---------1     1---------0     0---------3
//     |         |     |         |     |         |     |         |
//     |         |     |         |     |         |     |         |
//     |         |     |         |     |         |     |         |
//     0---------1     3---------0     2---------3     1---------2
//
//        rot=0           rot=1           rot=2           rot=3

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

   rotated_faces[0][0][0] =  1; rotated_faces[0][0][1] =  0; rotated_faces[0][1][0] =  0;   rotated_faces[0][1][1] =  1;
   rotated_faces[1][0][0] =  0; rotated_faces[1][0][1] = -1; rotated_faces[1][1][0] =  1;   rotated_faces[1][1][1] =  0;
   rotated_faces[2][0][0] = -1; rotated_faces[2][0][1] =  0; rotated_faces[2][1][0] =  0;   rotated_faces[2][1][1] = -1;
   rotated_faces[3][0][0] =  0; rotated_faces[3][0][1] =  1; rotated_faces[3][1][0] = -1;   rotated_faces[3][1][1] =  0;

//               7-------------6
//             . .           . |
//           .   .         .   |
//          4-------------5    |
//          |    .        |    |
//          |    3 . . . .|. . 2
//          |  .          |  .
//          |.            |.
//          0-------------1

#ifdef USE_SILO

   silo_zonetype = DB_ZONETYPE_HEX;
   zv_silo[0][0] = 0; zv_silo[0][1] = 0; zv_silo[1][0] = 0; zv_silo[1][1] = 1;
   zv_silo[2][0] = 0; zv_silo[2][1] = 2; zv_silo[3][0] = 0; zv_silo[3][1] = 3;
   zv_silo[4][0] = 1; zv_silo[4][1] = 0; zv_silo[5][0] = 1; zv_silo[5][1] = 1;
   zv_silo[6][0] = 1; zv_silo[6][1] = 2; zv_silo[7][0] = 1; zv_silo[7][1] = 3;

#endif

};

};

#endif
