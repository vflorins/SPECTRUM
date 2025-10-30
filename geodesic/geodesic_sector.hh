/*!
\file geodesic_sector.hh
\brief Declares sector class, a spherical polygon divided by grid lines
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_GEODESIC_SECTOR_HH
#define SPECTRUM_GEODESIC_SECTOR_HH

#include <cstdint>

#include "geodesic/polygonal_addressing.hh"

namespace Spectrum {

//! Mask for nonexisting element (FEV)
#define GEOELM_NEXI 0x0001

//! Mask for interior sector elements with boundary (FEV)
#define GEOELM_INTR 0x0002

//! Mask for ghost elements (FEV)
#define GEOELM_GHST 0x0004

//! Mask for the sector boundary (EV)
#define GEOELM_BNDR 0x0008

//! Smallest number of ghost t-faces
#define min_ghost_width 1

//! Largest number of ghost t-faces
#define max_ghost_width 4

//! Smallest ratio between block width and ghost width
#define min_width_to_ghost 4

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A class describing sectors in the mesh
\author Vladimir Florinski

This is a base class for grid blocks that provides access to the sector numbering schemes for the faces, edges, and vertices and the mesh connectivity lists. The numbering schemes are called Triangular Addressing Scheme (TAS) and Quad Addressing Scheme (QAS) corresponding to "verts_per_face" equal to 3 and 4, respectively. This class has no concept of coordinates and no possibility of singular corners (those are added in the "GridBlock" class).
*/
template <int verts_per_face>
class GeodesicSector : public PolygonalAddressing<verts_per_face>
{
protected:

   using PolygonalAddressing<verts_per_face>::edges_per_vert;
   using PolygonalAddressing<verts_per_face>::cardinal_directions;
   using PolygonalAddressing<verts_per_face>::square_fill;
   using PolygonalAddressing<verts_per_face>::edge_dimax;
   using PolygonalAddressing<verts_per_face>::edge_djmax;
   using PolygonalAddressing<verts_per_face>::vert_vert;
   using PolygonalAddressing<verts_per_face>::vert_edge;
   using PolygonalAddressing<verts_per_face>::vert_face;
   using PolygonalAddressing<verts_per_face>::edge_vert;
   using PolygonalAddressing<verts_per_face>::edge_face;
   using PolygonalAddressing<verts_per_face>::face_vert;
   using PolygonalAddressing<verts_per_face>::face_edge;
   using PolygonalAddressing<verts_per_face>::face_face;
   using PolygonalAddressing<verts_per_face>::EdgeCount;
   using PolygonalAddressing<verts_per_face>::FaceCount;
   using PolygonalAddressing<verts_per_face>::VertCount;

//! Length of a side of the sector (in edges)
   int side_length = -1;

//! Length of a side with ghost cells
   int total_length;

//! Width of the ghost cell layer outside the sector
   int ghost_width;

//! Number of vertices in the sector
   int n_verts;

//! Number of edges in the sector
   int n_edges;

//! Number of t-faces in the sector
   int n_faces;

//! Total number of vertices, including ghost vertices
   int n_verts_withghost;

//! Total number of edges, including ghost edges
   int n_edges_withghost;

//! Total number of t-faces, including ghost faces
   int n_faces_withghost;

//! Sector index vertex map
   int** vert_index_sector = nullptr;

//! Sector index edge map
   int** edge_index_sector[cardinal_directions] = {nullptr};
   
//! Sector index face map
   int** face_index_sector = nullptr;

//! Sector reverse vertex lookup - first  index (i)
   int* vert_index_i = nullptr;

//! Sector reverse vertex lookup - second index (j)
   int* vert_index_j = nullptr;

//! Sector reverse edge lookup - first  index (i)
   int* edge_index_i = nullptr;

//! Sector reverse edge lookup - second index (j)
   int* edge_index_j = nullptr;

//! Sector reverse face lookup - first  index (i)
   int* face_index_i = nullptr;

//! Sector reverse face lookup - second index (j)
   int* face_index_j = nullptr;

//! Vertex-vertex connectivity array
   int** vv_local = nullptr;

//! Vertex-edge connectivity array
   int** ve_local = nullptr;

//! Vertex-face connectivity array
   int** vf_local = nullptr;

//! Edge-vertex connectivity array
   int** ev_local = nullptr;

//! Edge-face connectivity array
   int** ef_local = nullptr;

//! Face-vertex connectivity array
   int** fv_local = nullptr;

//! Face-edge connectivity array
   int** fe_local = nullptr;

//! Face-face connectivity array
   int** ff_local = nullptr;

//! A mask for vertices
   uint16_t* vert_mask = nullptr;

//! A mask for edges
   uint16_t* edge_mask = nullptr;

//! A mask for faces
   uint16_t* face_mask = nullptr;

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Maximum value of the vertex j-index for the given vertex i-index
   int MaxVertJ(int len, int i) const;

//! Maximum value of the vertex j-index for the given vertex i-index for a sub-block
   int MaxVertJ(std::pair<int, int> base_vertex, int len, int i) const;

//! Maximum value of the face j-index for the given face i-index
   int MaxFaceJ(int len, int i) const;

//! Maximum value of the face j-index for the given face i-index for a sub-block
   int MaxFaceJ(std::pair<int, int> base_vertex, int len, int i) const;

//! Determine whether a vertex is a corner of a sub-block
   int CornerVert(std::pair<int, int> base_vertex, int len, int i, int j) const;

//! Determine whether a vertex is on the line extending a side of a sub-block
   int SideLine(std::pair<int, int> base_vertex, int len, int i, int j) const;

//! Determine whether a vertex is on the line extending a side of the sector's interior
   int SideLineOfSector(int i, int j) const;

//! Determine whether a vertex is on the boundary of a sub-block
   int BoundaryVert(std::pair<int, int> base_vertex, int len, int i, int j) const;

//! Determine whether a vertex is on the boundary of the sector's interior
   int BoundaryVertOfSector(int i, int j) const;

//! Determine whether an edge is on the boundary of a sub-block
   int BoundaryEdge(std::pair<int, int> base_vertex, int len, int etype, int i, int j) const;

//! Determine whether an edge is on the boundary of the sector's interior
   int BoundaryEdgeOfSector(int etype, int i, int j) const;

//! Determine whether a vertex belongs to a sub-block
   bool IsInteriorVert(std::pair<int, int> base_vertex, int len, int i, int j) const;

//! Determine whether a vertex belongs to the sector's interior (possibly on the boundary)
   bool IsInteriorVertOfSector(int i, int j) const;

//! Determine whether an edge belongs to a sub-block
   bool IsInteriorEdge(std::pair<int, int> base_vertex, int len, int etype, int i, int j) const;

//! Determine whether an edge is in the sector interior's (possibly on the boundary)
   bool IsInteriorEdgeOfSector(int etype, int i, int j) const;

//! Determine whether a face belongs to a sub-block
   bool IsInteriorFace(std::pair<int, int> base_vertex, int len, int i, int j) const;

//! Determine whether a face is in the sector's interior
   bool IsInteriorFaceOfSector(int i, int j) const;

//! Check if the vertex is outside the mesh (TAS only)
   bool IsClippedVert(int i, int j) const;

//! Check if the edge is outside the mesh (TAS only)
   bool IsClippedEdge(int etype, int i, int j) const;

//! Check if the face is outside the mesh (TAS only)
   bool IsClippedFace(int i, int j) const;

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Compute the vertex, edge, and face indices in the TAS/QAS
   void ComputeIndices(void);

//! Compute the local VV connectivity array
   void VertVertConn(void);

//! Compute the local VE connectivity array
   void VertEdgeConn(void);

//! Compute the local VF connectivity array
   void VertFaceConn(void);

//! Compute the local EV connectivity array
   void EdgeVertConn(void);

//! Compute the local EF connectivity array
   void EdgeFaceConn(void);

//! Compute the local FV connectivity array
   void FaceVertConn(void);

//! Compute the local FE connectivity array
   void FaceEdgeConn(void);

//! Compute the local FF connectivity array
   void FaceFaceConn(void);

public:

//! Default constructor
   GeodesicSector(void);

//! Copy constructor
   GeodesicSector(const GeodesicSector& other);

//! Move constructor
   GeodesicSector(GeodesicSector&& other) noexcept;

//! Move constructor
   SPECTRUM_DEVICE_FUNC GeodesicSector(GeodesicSector&& other);

//! Constructor with arguments
   GeodesicSector(int width, int wghost);

//! Destructor
   ~GeodesicSector(void);

//! Allocate memory and set up connectivity between mesh elements
   void SetDimensions(int width, int wghost, bool construct);

//! Free all dynamically allocated memory
   void FreeStorage(void);

//! Return the number of t-faces excluding ghost faces
   int InteriorFaces(void) const {return n_faces;};

//! Return the number of t-faces including ghost faces
   int TotalFaces(void) const {return n_faces_withghost;};

//! Return the number of vertices excluding ghost vertices
   int InteriorVerts(void) const {return n_verts;};

//! Return the number of vertices including ghost vertices
   int TotalVerts(void) const {return n_verts_withghost;};

#ifdef GEO_DEBUG

//! Prints the TAS/QAS address table
   void PrintAddresses(int type) const;

//! Prints a mesh connectivity array
   void PrintConn(int type) const;

//! Prints the mask
   void PrintMask(int type) const;

#endif

};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// GeodesicSector inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/08/2024
\param[in] len Side length
\param[in] i   i-index
\return Largest j-index for this i
*/
template <>
inline int GeodesicSector<3>::MaxVertJ(int len, int i) const
{
   return i;
};

/*!
\author Vladimir Florinski
\date 05/08/2024
\param[in] len Side length
\param[in] i   i-index
\return Largest j-index for this i
*/
template <>
inline int GeodesicSector<4>::MaxVertJ(int len, int i) const
{
   return len;
};

/*!
\author Vladimir Florinski
\date 05/21/2024
\param[in] base_vertex Base vertex of the sub-block
\param[in] len         Sub-block side length
\param[in] i           i-index
\return Largest j-index for this i
*/
template <>
inline int GeodesicSector<3>::MaxVertJ(std::pair<int, int> base_vertex, int len, int i) const
{
   return base_vertex.second - base_vertex.first + i;
};

/*!
\author Vladimir Florinski
\date 05/21/2024
\param[in] base_vertex Base vertex of the sub-block
\param[in] len         Length of the side of the sub-block
\param[in] i           i-index
\return Largest j-index for this i
*/
template <>
inline int GeodesicSector<4>::MaxVertJ(std::pair<int, int> base_vertex, int len, int i) const
{
   return base_vertex.second + len;
};

/*!
\author Vladimir Florinski
\date 05/08/2024
\param[in] len Side length
\param[in] i   i-index
\return Largest j-index for this i
*/
template <>
inline int GeodesicSector<3>::MaxFaceJ(int len, int i) const
{
   return square_fill * i;
};

/*!
\author Vladimir Florinski
\date 05/08/2024
\param[in] len Side length
\param[in] i   i-index
\return Largest j-index for this i
*/
template <>
inline int GeodesicSector<4>::MaxFaceJ(int len, int i) const
{
   return len - 1;
};

/*!
\author Vladimir Florinski
\date 05/21/2024
\param[in] base_vertex Base vertex of the sub-block
\param[in] len         Sub-block side length
\param[in] i           i-index
\return Largest j-index for this i
*/
template <>
inline int GeodesicSector<3>::MaxFaceJ(std::pair<int, int> base_vertex, int len, int i) const
{
   return square_fill * (base_vertex.second - base_vertex.first + i);
};

/*!
\author Vladimir Florinski
\date 05/21/2024
\param[in] base_vertex Base vertex of the sub-block
\param[in] len         Length of the side of the sub-block
\param[in] i           i-index
\return Largest j-index for this i
*/
template <>
inline int GeodesicSector<4>::MaxFaceJ(std::pair<int, int> base_vertex, int len, int i) const
{
   return base_vertex.second + len - 1;
};

/*!
\author Vladimir Florinski
\date 05/21/2024
\param[in] base_vertex Base vertex of the sub-block
\param[in] len         Length of the side of the sub-block
\param[in] i           First index
\param[in] j           Second index
\return Corner number (-1 if not a corner vertex)
*/
template <int verts_per_face>
inline int GeodesicSector<verts_per_face>::CornerVert(std::pair<int, int> base_vertex, int len, int i, int j) const
{
   if (i == base_vertex.first) {
      if (j == base_vertex.second) return 0;
      else if (j == MaxVertJ(base_vertex, len, i)) return 3;
      else return -1;
   }
   else if (i == base_vertex.first + len) {
      if (j == base_vertex.second) return 1;
      else if (j == MaxVertJ(base_vertex, len, i)) return 2;
      else return -1;
   }
   else return -1;
};

/*!
\author Vladimir Florinski
\date 05/21/2024
\param[in] base_vertex Base vertex of the sub-block
\param[in] len         Length of the side of the sub-block
\param[in] i           First index
\param[in] j           Second index
\return Side number (-1 if not on the side line)
*/
template <int verts_per_face>
inline int GeodesicSector<verts_per_face>::SideLine(std::pair<int, int> base_vertex, int len, int i, int j) const
{
   int j0 = base_vertex.second;
   int i1 = base_vertex.first + len;
   int j2 = MaxVertJ(base_vertex, len, i);
   int i3 = base_vertex.first;

// Corner vertices belong to two side lines. In that case the side lying in the CC direction from the vertex is returned.

//                 1   2
//                  1 2                      3         1
//                   2                       3         1
//                  2 1                  2 2 3 2 2 2 2 2 2 2
//                 2   1                     3         1
//                2     1                    3         1
//               2       1                   3         1
//          0 0 0 0 0 0 0 1 0 0          0 0 0 0 0 0 0 1 0 0
//             2           1                 3         1
//            2             1                3         1

   if      ((j == j0) && (i != i1)) return 0;
   else if ((i == i1) && (j != j2)) return 1;
   else if ((j == j2) && (i != i3)) return 2;
   else if ((verts_per_face == 4) && (i == i3)) return 3;
   else return -1;
};

/*!
\author Vladimir Florinski
\date 05/07/2024
\param[in] i First index
\param[in] j Second index
\return Side number (-1 if not on the side line)
*/
template <int verts_per_face>
inline int GeodesicSector<verts_per_face>::SideLineOfSector(int i, int j) const
{
   return SideLine(std::make_pair(square_fill * ghost_width, ghost_width), side_length, i, j);
};

/*!
\author Vladimir Florinski
\date 05/21/2024
\param[in] base_vertex Base vertex of the sub-block
\param[in] len         Length of the side of the sub-block
\param[in] i           First index
\param[in] j           Second index
\return Side number (-1 if not on the boundary)
*/
template <int verts_per_face>
inline int GeodesicSector<verts_per_face>::BoundaryVert(std::pair<int, int> base_vertex, int len, int i, int j) const
{
   int side = SideLine(base_vertex, len, i, j);

   if (side % 2) {
   if ((j >= base_vertex.second) && (j <= MaxVertJ(base_vertex, len, i))) return side;
      else return -1;
   }
   else {
      if ((i >= base_vertex.first) && (i <= base_vertex.first + len)) return side;
      else return -1;
   };
};

/*!
\author Vladimir Florinski
\date 05/07/2024
\param[in] i First index
\param[in] j Second index
\return Side number (-1 if not on the boundary)
*/
template <int verts_per_face>
inline int GeodesicSector<verts_per_face>::BoundaryVertOfSector(int i, int j) const
{
   return BoundaryVert(std::make_pair(square_fill * ghost_width, ghost_width), side_length, i, j);
};

/*
\author Vladimir Florinski
\date 05/21/2024
\param[in] base_vertex Base vertex of the sub-block
\param[in] len         Length of the side of the sub-block
\param[in] etype       Edge type
\param[in] i           First index
\param[in] j           Second index
\return Side number (-1 if not on the boundary)
*/
template <int verts_per_face>
inline int GeodesicSector<verts_per_face>::BoundaryEdge(std::pair<int, int> base_vertex, int len, int etype, int i, int j) const
{
   int ii[2] = {i + edge_vert[etype][0][0], i + edge_vert[etype][1][0]};
   int jj[2] = {j + edge_vert[etype][0][1], j + edge_vert[etype][1][1]};

   int side[2] = {BoundaryVert(base_vertex, len, ii[0], jj[0]), BoundaryVert(base_vertex, len, ii[1], jj[1])};
   int corner[2];

// Both vertices lie on the same side of the sub-block. This also covers the case when neither vertex does.
   if (side[0] == side[1]) return side[0];

// If the vertices are on different sides, check the corner status. Return the side of the vertex that is not in the corner.
   else if ((side[0] != -1) && (side[1] != -1)) {
      corner[0] = CornerVert(base_vertex, len, ii[0], jj[0]);
      corner[1] = CornerVert(base_vertex, len, ii[1], jj[1]);
      if (corner[0] != -1) return side[1];
      else if (corner[1] != -1) return side[0];
      else return -1;
   }         
   else return -1;
};

/*
\author Vladimir Florinski
\date 05/07/2024
\param[in] etype Edge type
\param[in] i     First index
\param[in] j     Second index
\return Side number (-1 if not on the boundary)
*/
template <int verts_per_face>
inline int GeodesicSector<verts_per_face>::BoundaryEdgeOfSector(int etype, int i, int j) const
{
   return BoundaryEdge(std::make_pair(square_fill * ghost_width, ghost_width), side_length, etype, i, j);
};

/*!
\author Vladimir Florinski
\date 05/21/2024
\param[in] base_vertex Base vertex of the sub-block
\param[in] len         Length of the side of the sub-block
\param[in] i           First index
\param[in] j           Second index
\return True if the vertex belongs to the sub-block
*/
template <int verts_per_face>
inline bool GeodesicSector<verts_per_face>::IsInteriorVert(std::pair<int, int> base_vertex, int len, int i, int j) const
{
   return (i >= base_vertex.first) && (i <= base_vertex.first + len)
       && (j >= base_vertex.second) && (j <= MaxVertJ(base_vertex, len, i));
};

/*!
\author Vladimir Florinski
\date 05/07/2024
\param[in] i First index
\param[in] j Second index
\return True if the vertex belongs to the sector's interior
*/
template <int verts_per_face>
inline bool GeodesicSector<verts_per_face>::IsInteriorVertOfSector(int i, int j) const
{
   return IsInteriorVert(std::make_pair(square_fill * ghost_width, ghost_width), side_length, i, j);
};

/*!
\author Vladimir Florinski
\date 05/21/2024
\param[in] base_vertex Base vertex of the sub-block
\param[in] len         Length of the side of the sub-block
\param[in] etype       Edge type
\param[in] i           First index
\param[in] j           Second index
\return True if the edge belongs to the sub-block
*/
template <int verts_per_face>
inline bool GeodesicSector<verts_per_face>::IsInteriorEdge(std::pair<int, int> base_vertex, int len, int etype, int i, int j) const
   {
      return (i >= base_vertex.first) && (i <= base_vertex.first + len + edge_dimax[etype])
          && (j >= base_vertex.second) && (j <= MaxVertJ(base_vertex, len, i) + edge_djmax[etype]);
   };

/*
\author Vladimir Florinski
\date 05/07/2024
\param[in] etype Edge type
\param[in] i     First index
\param[in] j     Second index
\return True if the edge is interior to the sector
*/
template <int verts_per_face>
inline bool GeodesicSector<verts_per_face>::IsInteriorEdgeOfSector(int etype, int i, int j) const
{
   return IsInteriorEdge(std::make_pair(square_fill * ghost_width, ghost_width), side_length, etype, i, j);
};

/*!
\author Vladimir Florinski
\date 05/21/2024
\param[in] base_vertex Base vertex of the sub-block
\param[in] len         Length of the side of the sub-block
\param[in] i           First index
\param[in] j           Second index
\return True if the face belongs to the sub-block
*/
template <int verts_per_face>
inline bool GeodesicSector<verts_per_face>::IsInteriorFace(std::pair<int, int> base_vertex, int len, int i, int j) const
{
   return (i >= base_vertex.first) && (i <= base_vertex.first + len - 1)
       && (j >= square_fill * base_vertex.second) && (j <= MaxFaceJ(base_vertex, len, i));
};

/*!
\author Vladimir Florinski
\date 05/07/2024
\param[in] i First index
\param[in] j Second index
\return True if the face is interior to the sector
*/
template <int verts_per_face>
inline bool GeodesicSector<verts_per_face>::IsInteriorFaceOfSector(int i, int j) const
{
   return IsInteriorFace(std::make_pair(square_fill * ghost_width, ghost_width), side_length, i, j);
};

/*!
\author Vladimir Florinski
\date 05/07/2024
\param[in] i First index
\param[in] j Second index
\return True if the vertex lies outside the mesh (does not exist)
*/
template <int verts_per_face>
inline bool GeodesicSector<verts_per_face>::IsClippedVert(int i, int j) const
{
   if ((i < 0) || (i > total_length) || (j < 0) || (j > MaxVertJ(total_length, i))) return true;
   if (verts_per_face != 3) return false;

   if (IsInteriorVert(std::make_pair(0, 0), ghost_width - 1, i, j)) return true;
   if (IsInteriorVert(std::make_pair(total_length - ghost_width + 1, 0), ghost_width - 1, i, j)) return true;
   if (IsInteriorVert(std::make_pair(total_length - ghost_width + 1, MaxVertJ(total_length, total_length - ghost_width + 1)), ghost_width - 1, i, j)) return true;

   return false;
};

/*!
\author Vladimir Florinski
\date 05/21/2024
\param[in] etype Edge type
\param[in] i     First index
\param[in] j     Second index
\return True if the edge lies outside the mesh (does not exist)
*/
template <int verts_per_face>
inline bool GeodesicSector<verts_per_face>::IsClippedEdge(int etype, int i, int j) const
{
   if ((i < 0) || (i > total_length + edge_dimax[etype]) || (j < 0) || (j > MaxVertJ(total_length, i) + edge_djmax[etype])) return true;
   if (verts_per_face != 3) return false;

   if (etype == 0) {
      if (IsInteriorEdge(std::make_pair(0, 0), ghost_width, etype, i, j)) return true;
      if (IsInteriorEdge(std::make_pair(total_length - ghost_width, 0), ghost_width, etype, i, j)) return true;
      if (IsInteriorEdge(std::make_pair(total_length - ghost_width + 1, MaxVertJ(total_length, total_length - ghost_width + 1)), ghost_width - 1, etype, i, j)) return true;
   }
   else if (etype == 1) {
      if (IsInteriorEdge(std::make_pair(0, 0), ghost_width - 1, etype, i, j)) return true;
      if (IsInteriorEdge(std::make_pair(total_length - ghost_width, 0), ghost_width, etype, i, j)) return true;
      if (IsInteriorEdge(std::make_pair(total_length - ghost_width, MaxVertJ(total_length, total_length - ghost_width)), ghost_width, etype, i, j)) return true;
   }
   else if (etype == 2) {
      if (IsInteriorEdge(std::make_pair(0, 0), ghost_width, etype, i, j)) return true;
      if (IsInteriorEdge(std::make_pair(total_length - ghost_width + 1, 0), ghost_width - 1, etype, i, j)) return true;
      if (IsInteriorEdge(std::make_pair(total_length - ghost_width, MaxVertJ(total_length, total_length - ghost_width)), ghost_width, etype, i, j)) return true;
   };

   return false;
};

/*!
\author Vladimir Florinski
\date 05/21/2024
\param[in] i First index
\param[in] j Second index
\return True if the face lies outside the mesh (does not exist)
*/
template <int verts_per_face>
inline bool GeodesicSector<verts_per_face>::IsClippedFace(int i, int j) const
{
   if ((i < 0) || (i > total_length - 1) || (j < 0) || (j > MaxFaceJ(total_length, i))) return true;
   if (verts_per_face != 3) return false;

   if (IsInteriorFace(std::make_pair(0, 0), ghost_width, i, j)) return true;
   if (IsInteriorFace(std::make_pair(total_length - ghost_width, 0), ghost_width, i, j)) return true;
   if (IsInteriorFace(std::make_pair(total_length - ghost_width, MaxVertJ(total_length, total_length - ghost_width)), ghost_width, i, j)) return true;

   return false;
};

};

#endif