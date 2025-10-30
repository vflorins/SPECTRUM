/*!
\file polygonal_addressing.hh
\brief A system of numbering vertices, edges, and faces in a polygonal sector
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_POLYGONAL_ADDRESSING_HH
#define SPECTRUM_POLYGONAL_ADDRESSING_HH

namespace Spectrum {

/*!
\brief A class describing the Triangular Addressing Scheme (TAS) and Quad Addressing Scheme (QAS)
\author Vladimir Florinski

This class's role is to fill in the bootstrap arrays ("vert_vert", "vert_edge", etc.) that are used later in the "GeodesicSector" class to generate the connectivity arrays ("vv_local", "ve_local", etc.) The bootstrap arrays have minimum utility beyond that, and it is recommended that they are not used in classes derived from "GeodesicSector".
*/
template <int verts_per_face>
class PolygonalAddressing
{
private:

//! Calculate the increments for each object type
   static void StaticSetup(void);

//! Structure to enable static initialization before "main()"
   struct StaticInitializer
   {
      StaticInitializer(void) {StaticSetup();};
   };
   
//! Static initializer
   static inline StaticInitializer si{};

protected:

//! Number of edges meeting at a vertex
   static constexpr int edges_per_vert = EdgesAtVert<verts_per_face>();

//! Rotational symmetry
   static constexpr int cardinal_directions = edges_per_vert / 2;

//! How many zones fill a quadrilateral
   static constexpr int square_fill = edges_per_vert / verts_per_face;

//! Use vertex macros to loop over edge i-index. Edges have the lower limit as vertices, but the upper limit could be off by one.
   static inline int edge_dimax[cardinal_directions];

//! Use vertex macros to loop over edge j-index. Edges have the lower limit as vertices, but the upper limit could be off by one.
   static inline int edge_djmax[cardinal_directions];

//! Vertex neighbor offsets for vertices
   static inline int vert_vert[edges_per_vert][2];

//! Edge neighbor offsets for vertices (i, j, etype)
   static inline int vert_edge[edges_per_vert][3];

//! Face neighbor offsets for vertices (i, j)
   static inline int vert_face[edges_per_vert][2];

//! Vertex neighbor offsets for edges of each type (i, j)
   static inline int edge_vert[cardinal_directions][2][2];

//! Face neighbor offsets for edges of each type (i, j)
   static inline int edge_face[cardinal_directions][2][2];

//! Vertex neighbor offsets for faces (i, j)
   static inline int face_vert[square_fill][verts_per_face][2];

//! Edge neighbor offsets for faces (i, j, etype)
   static inline int face_edge[square_fill][verts_per_face][3];

//! Face neighbor offsets for faces
   static inline int face_face[square_fill][verts_per_face][2];

//! Default constructor (protected)
   PolygonalAddressing(void) = default;

public:

//! Edge count for a given side length
   static constexpr int EdgeCount(int len);

//! Face count for a given side length
   static constexpr int FaceCount(int len);

//! Vertex count for a given side length - Euler's formula for planar graphs
   static constexpr int VertCount(int len);

#ifdef GEO_DEBUG
//! Test the static iniitalization
   static constexpr int GetVertVert(int ie, int ij) {return vert_vert[ie][ij];};
#endif

};

/*!
\author Vladimir Florinski
\date 01/08/20256
\param[in] other Object to move into this
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC inline PolygonalAddressing<verts_per_face>::PolygonalAddressing(PolygonalAddressing<verts_per_face>&& other)
                                                               : PolygonalAddressing<verts_per_face>()
{
};

/*!
\author Vladimir Florinski
\date 05/06/2024
\param[in] len Length of a side
\return Number of edges in a sector whose side is "len"
*/
template <int verts_per_face>
inline constexpr int PolygonalAddressing<verts_per_face>::EdgeCount(int len)
{
   return cardinal_directions * len * (len + 1) / square_fill;
};

/*!
\author Vladimir Florinski
\date 05/06/2024
\param[in] len Length of a side
\return Number of faces in a sector whose side is "len"
*/
template <int verts_per_face>
inline constexpr int PolygonalAddressing<verts_per_face>::FaceCount(int len)
{
   return len * len;
};

/*!
\author Vladimir Florinski
\date 05/06/2024
\param[in] len Length of a side
\return Number of vertices in a sector whose side is "len"
*/
template <int verts_per_face>
inline constexpr int PolygonalAddressing<verts_per_face>::VertCount(int len)
{
   return EdgeCount(len) - FaceCount(len) + 1;
};

/*!
\author Vladimir Florinski
\date 05/07/2024
*/
template <>
inline constexpr void PolygonalAddressing<3>::StaticSetup(void)
{
// Indexing scheme: vertex "V" and edges "T0", "T1", and "T2" have the same i and j coordinates. Face "F" (upward) has coordinates (i,2j). Face "S" (downward) has coordinates (i,2j+1). When obtaining edge and vertex coorindates from face coordinates, divide the second coordinate by two, so both "F" and "S" will point to (i,j).

/*          -----------
             \       / \
              T1 S  T2  \
               \   /  F  \
                \ /       \
                 *----T0----
                 V
*/
   edge_dimax[0] = -1; edge_dimax[1] =  0; edge_dimax[2] = -1;
   edge_djmax[0] =  0; edge_djmax[1] = -1; edge_djmax[2] =  0;

/*              vertex-vertex                   vertex-edge                    vertex-face
  
                 4         3                    -         -     	       -----------
                                                 \       /      	      / \	/ \
                                                  4     3       	     /   \  4  /   \
                                                   \   /        	    /  5  \   /  3  \
                                                    \ /         	   /	   \ /       \
            5         X         2          ----5-----X-----2----	  ----------X----------
                                                    / \         	   \	   / \       /
                                                   /   \        	    \  0  /   \  2  /
                                                  0     1       	     \   /  1  \   /
                                                 /       \      	      \ /	\ /
                 0         1                    -         -     	       -----------
*/
   vert_vert[0][0] = -1; vert_vert[0][1] = -1;
   vert_vert[1][0] =  0; vert_vert[1][1] = -1;
   vert_vert[2][0] =  1; vert_vert[2][1] =  0;
   vert_vert[3][0] =  1; vert_vert[3][1] =  1;
   vert_vert[4][0] =  0; vert_vert[4][1] =  1;
   vert_vert[5][0] = -1; vert_vert[5][1] =  0;

   vert_edge[0][0] = -1; vert_edge[0][1] = -1; vert_edge[0][2] =  2;
   vert_edge[1][0] =  0; vert_edge[1][1] = -1; vert_edge[1][2] =  1;
   vert_edge[2][0] =  0; vert_edge[2][1] =  0; vert_edge[2][2] =  0;
   vert_edge[3][0] =  0; vert_edge[3][1] =  0; vert_edge[3][2] =  2;
   vert_edge[4][0] =  0; vert_edge[4][1] =  0; vert_edge[4][2] =  1;
   vert_edge[5][0] = -1; vert_edge[5][1] =  0; vert_edge[5][2] =  0;

   vert_face[0][0] = -1; vert_face[0][1] = -1;
   vert_face[1][0] = -1; vert_face[1][1] = -2;
   vert_face[2][0] =  0; vert_face[2][1] = -1;
   vert_face[3][0] =  0; vert_face[3][1] =  0;
   vert_face[4][0] =  0; vert_face[4][1] =  1;
   vert_face[5][0] = -1; vert_face[5][1] =  0;

/*                           edge-vertex                                                 edge-face 
  
               type 0            type 1          type 2             type 0                 type 1                    type 2
  
                                 1                    0               -                    -----------          -----------
                                  \                  /               / \                  / \       /            \       / \
                                   \                /               /   \                /   \  1  /              \  1  /   \
                                    \              /               /  0  \              /  0  \   /                \   /  0  \
                                     \            /               /       \            /       \ /                  \ /       \
            0---------1               0          1               -----------          -----------                    -----------
                                                                  \       /
                                                                   \  1  /
                                                                    \   /
                                                                     \ /
                                                                      -
*/
   edge_vert[0][0][0] =  0; edge_vert[0][0][1] =  0; edge_vert[0][1][0] =  1; edge_vert[0][1][1] =  0;
   edge_vert[1][0][0] =  0; edge_vert[1][0][1] =  0; edge_vert[1][1][0] =  0; edge_vert[1][1][1] =  1;
   edge_vert[2][0][0] =  1; edge_vert[2][0][1] =  1; edge_vert[2][1][0] =  0; edge_vert[2][1][1] =  0;

   edge_face[0][0][0] =  0; edge_face[0][0][1] =  0; edge_face[0][1][0] =  0; edge_face[0][1][1] = -1;
   edge_face[1][0][0] = -1; edge_face[1][0][1] =  0; edge_face[1][1][0] =  0; edge_face[1][1][1] =  1;
   edge_face[2][0][0] =  0; edge_face[2][0][1] =  0; edge_face[2][1][0] =  0; edge_face[2][1][1] =  1;

/*               face-vertex                      face-edge                                  face-face
  
                 2     1---------0               -     -----0-----          ---------------------          -
                / \     \       /               / \     \       /            \       / \       /          / \
               /   \     \     /               /   \     1     2              \  2  /   \  1  /          /   \
              /     \     \   /               2     1     \   /                \   /     \   /          /  0  \
             /       \     \ /               /       \     \ /                  \ /       \ /          /       \
            0---------1     2               -----0-----     -                    -----------          -----------
                                                                                  \   0   /          / \       / \
                                                                                   \     /          /   \     /   \
                                                                                    \   /          /  1  \   /  2  \
                                                                                     \ /          /       \ /       \
                                                                                      -          ---------------------
*/
// Upward faces are element [0] and downward faces are element [1]. The correct set is chosen based on (j % 2).

   face_vert[0][0][0] =  0; face_vert[0][0][1] =  0;
   face_vert[0][1][0] =  1; face_vert[0][1][1] =  0;
   face_vert[0][2][0] =  1; face_vert[0][2][1] =  1;
   face_vert[1][0][0] =  1; face_vert[1][0][1] =  1;
   face_vert[1][1][0] =  0; face_vert[1][1][1] =  1;
   face_vert[1][2][0] =  0; face_vert[1][2][1] =  0;

   face_edge[0][0][0] =  0; face_edge[0][0][1] =  0; face_edge[0][0][2] =  0;
   face_edge[0][1][0] =  1; face_edge[0][1][1] =  0; face_edge[0][1][2] =  1;
   face_edge[0][2][0] =  0; face_edge[0][2][1] =  0; face_edge[0][2][2] =  2;
   face_edge[1][0][0] =  0; face_edge[1][0][1] =  1; face_edge[1][0][2] =  0;
   face_edge[1][1][0] =  0; face_edge[1][1][1] =  0; face_edge[1][1][2] =  1;
   face_edge[1][2][0] =  0; face_edge[1][2][1] =  0; face_edge[1][2][2] =  2;

   face_face[0][0][0] =  0; face_face[0][0][1] = -1;
   face_face[0][1][0] =  1; face_face[0][1][1] =  1;
   face_face[0][2][0] =  0; face_face[0][2][1] =  1;
   face_face[1][0][0] =  0; face_face[1][0][1] =  1;
   face_face[1][1][0] = -1; face_face[1][1][1] = -1;
   face_face[1][2][0] =  0; face_face[1][2][1] = -1;
};

/*!
\author Vladimir Florinski
\date 05/07/2024
*/
template <>
inline constexpr void PolygonalAddressing<4>::StaticSetup(void)
{
// Indexing scheme: vertex "V", edges "T0" and "T1", and face "F" have the same i and j coordinates.

/*          -----------
            |         |
            T1   F    |
            |         |
            *----T0----
            V
*/
   edge_dimax[0] = -1; edge_dimax[1] =  0;
   edge_djmax[0] =  0; edge_djmax[1] = -1;

/*              vertex-vertex                   vertex-edge                    vertex-face
  
                      2                              -                   ---------------------
                                                     |                   |         |         |
                                                     2                   |    3    |    2    |
                                                     |                   |         |         |
            3         +         1          -----3----+----1----          |---------+---------|
                                                     |                   |         |         |
                                                     0                   |    0    |    1    |
                                                     |                   |         |         |
                      0                              -                   ---------------------
*/
   vert_vert[0][0] =  0; vert_vert[0][1] = -1;
   vert_vert[1][0] =  1; vert_vert[1][1] =  0;
   vert_vert[2][0] =  0; vert_vert[2][1] =  1;
   vert_vert[3][0] = -1; vert_vert[3][1] =  0;

   vert_edge[0][0] =  0; vert_edge[0][1] = -1; vert_edge[0][2] =  1;
   vert_edge[1][0] =  0; vert_edge[1][1] =  0; vert_edge[1][2] =  0;
   vert_edge[2][0] =  0; vert_edge[2][1] =  0; vert_edge[2][2] =  1;
   vert_edge[3][0] = -1; vert_edge[3][1] =  0; vert_edge[3][2] =  0;

   vert_face[0][0] = -1; vert_face[0][1] = -1;
   vert_face[1][0] =  0; vert_face[1][1] = -1;
   vert_face[2][0] =  0; vert_face[2][1] =  0;
   vert_face[3][0] = -1; vert_face[3][1] =  0;

/*                edge-vertex                               edge-face 
  
                                 1          -----------          ---------------------
                                 |          |         |          |         |         |
                                 |          |    0    |          |    0    |    1    |
                                 |          |         |          |         |         |
            0---------1          0          -----------          ---------------------
                                            |         |
                                            |    1    |
                                            |         |
                                            -----------
*/
   edge_vert[0][0][0] =  0; edge_vert[0][0][1] =  0; edge_vert[0][1][0] =  1; edge_vert[0][1][1] =  0;
   edge_vert[1][0][0] =  0; edge_vert[1][0][1] =  0; edge_vert[1][1][0] =  0; edge_vert[1][1][1] =  1;

   edge_face[0][0][0] =  0; edge_face[0][0][1] =  0; edge_face[0][1][0] =  0; edge_face[0][1][1] = -1;
   edge_face[1][0][0] = -1; edge_face[1][0][1] =  0; edge_face[1][1][0] =  0; edge_face[1][1][1] =  0;

/*          face-vertex           face-edge                      face-face
  
                                                                -----------
                                                                |         |
                                                                |    2    |
                                                                |         |
            3---------2          -----2-----          ----------+---------+----------
            |         |          |         |          |         |         |         |
            |         |          3         1          |    3    |         |    1    |
            |         |          |         |          |         |         |         |
            0---------1          -----0-----          ----------+---------+----------
                                                                |         |
                                                                |    0    |
                                                                |         |
                                                                -----------
*/
   face_vert[0][0][0] =  0; face_vert[0][0][1] =  0;
   face_vert[0][1][0] =  1; face_vert[0][1][1] =  0;
   face_vert[0][2][0] =  1; face_vert[0][2][1] =  1;
   face_vert[0][3][0] =  0; face_vert[0][3][1] =  1;

   face_edge[0][0][0] =  0; face_edge[0][0][1] =  0; face_edge[0][0][2] =  0;
   face_edge[0][1][0] =  1; face_edge[0][1][1] =  0; face_edge[0][1][2] =  1;
   face_edge[0][2][0] =  0; face_edge[0][2][1] =  1; face_edge[0][2][2] =  0;
   face_edge[0][3][0] =  0; face_edge[0][3][1] =  0; face_edge[0][3][2] =  1;

   face_face[0][0][0] =  0; face_face[0][0][1] = -1;
   face_face[0][1][0] =  1; face_face[0][1][1] =  0;
   face_face[0][2][0] =  0; face_face[0][2][1] =  1;
   face_face[0][3][0] = -1; face_face[0][3][1] =  0;
};

template class PolygonalAddressing<3>;
template class PolygonalAddressing<4>;

};

#endif