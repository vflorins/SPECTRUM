/*!
\file geodesic_sector.cc
\brief Implements the GeodesicSector class
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include <cstring>
#include "common/definitions.hh"
#include "common/print_warn.hh"
#include "geodesic/geodesic_sector.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// GeodesicSector public methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/28/2024
\param[in] other Object to initialize from
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC GeodesicSector<verts_per_face>::GeodesicSector(const GeodesicSector& other)
                                                   : PolygonalAddressing<verts_per_face>()
{
   if (other.side_length != -1) {
      SetDimensions(other.side_length, other.ghost_width, true);
   };
};

/*!
\author Vladimir Florinski
\date 07/23/2019
\param[in] width  Length of the side, without ghost cells
\param[in] wgohst Width of the ghost cell layer outside the sector
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC GeodesicSector<verts_per_face>::GeodesicSector(int width, int wghost)
                                                   : PolygonalAddressing<verts_per_face>()
{
   SetDimensions(width, wghost, true);
};

/*!
\author Vladimir Florinski
\date 07/23/2019
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC GeodesicSector<verts_per_face>::~GeodesicSector()
{
   FreeStorage();
};

/*!
\author Vladimir Florinski
\date 06/28/2024
\param[in] width     Length of the side, without ghost cells
\param[in] wgohst    Width of the ghost cell layer outside the sector
\param[in] construct Unused, but included for consistency with derived classes
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC void GeodesicSector<verts_per_face>::SetDimensions(int width, int wghost, bool construct)
{
   if ((width < min_block_width) || (width > max_block_width) || (wghost < min_ghost_width) || (wghost > max_ghost_width)) {
      PrintError(__FILE__, __LINE__, "Cannot allocate a sector with these dimensions", true);
      return;
   };

// Free up storage because this could be a repeat call
   FreeStorage();

// Set up the grid dimensions
   side_length = width;
   ghost_width = wghost;
   total_length = side_length + (1 + square_fill) * ghost_width;

   n_verts = VertCount(side_length);
   n_verts_withghost = VertCount(total_length);
   n_edges = EdgeCount(side_length);
   n_edges_withghost = EdgeCount(total_length);
   n_faces = FaceCount(side_length);
   n_faces_withghost = FaceCount(total_length);

// Index arrays are used to quickly access parts of the block for exchanges across boundaries and input/output.
   vert_index_sector = Create2D<int>(total_length + 1, total_length + 1);
   for (auto etype = 0; etype < cardinal_directions; etype++) {
      edge_index_sector[etype] = Create2D<int>(total_length + 1, total_length + 1);
   };
   face_index_sector = Create2D<int>(total_length, square_fill * total_length);

// Reverse lookup arrays
   vert_index_i = new int[n_verts_withghost];
   vert_index_j = new int[n_verts_withghost];
   edge_index_i = new int[n_edges_withghost];
   edge_index_j = new int[n_edges_withghost];
   face_index_i = new int[n_faces_withghost];
   face_index_j = new int[n_faces_withghost];

// Allocate memory for connectivity arrays
   vv_local = Create2D<int>(n_verts_withghost, edges_per_vert);
   ve_local = Create2D<int>(n_verts_withghost, edges_per_vert);
   vf_local = Create2D<int>(n_verts_withghost, edges_per_vert);
   ev_local = Create2D<int>(n_edges_withghost, 2);
   ef_local = Create2D<int>(n_edges_withghost, 2);
   fv_local = Create2D<int>(n_faces_withghost, verts_per_face);
   fe_local = Create2D<int>(n_faces_withghost, verts_per_face);
   ff_local = Create2D<int>(n_faces_withghost, verts_per_face);

// Allocate and zero off the masks
   vert_mask = new uint16_t[n_verts_withghost];
   edge_mask = new uint16_t[n_edges_withghost];
   face_mask = new uint16_t[n_faces_withghost];
   memset(vert_mask, 0x0, n_verts_withghost * sizeof(uint16_t));
   memset(edge_mask, 0x0, n_edges_withghost * sizeof(uint16_t));
   memset(face_mask, 0x0, n_faces_withghost * sizeof(uint16_t));

// Compute addresses and build connectivity lists
   ComputeIndices();

   VertVertConn();
   VertEdgeConn();
   VertFaceConn();
   EdgeVertConn();
   EdgeFaceConn();
   FaceVertConn();
   FaceEdgeConn();
   FaceFaceConn();
};

/*!
\author Vladimir Florinski
\date 07/23/2019
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC void GeodesicSector<verts_per_face>::FreeStorage(void)
{
// Free up masks
   delete[] vert_mask;
   delete[] edge_mask;
   delete[] face_mask;

// Free up connectivity memory
   Delete2D(vv_local);
   Delete2D(ve_local);
   Delete2D(vf_local);
   Delete2D(ev_local);
   Delete2D(ef_local);
   Delete2D(fv_local);
   Delete2D(fe_local);
   Delete2D(ff_local);

// Free up index arrays
   delete[] vert_index_i;
   delete[] vert_index_j;
   delete[] edge_index_i;
   delete[] edge_index_j;
   delete[] face_index_i;
   delete[] face_index_j;

   Delete2D(face_index_sector);
   for (auto etype = 0; etype < cardinal_directions; etype++) Delete2D(edge_index_sector[etype]);
   Delete2D(vert_index_sector);

// This marks this sector as un-allocated
   side_length = -1;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// GeodesicSector protected methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/21/2024
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC void GeodesicSector<verts_per_face>::ComputeIndices(void)
{
   int i, j, etype, vert, edge, face;
   std::pair base_vert = std::make_pair(0, 0);

// Count all vertices (even non-existing)
//
//                    5 
//                   / \                   2---------5---------8
//                  /   \                  |         |         |
//                 /     \                 |         |         |
//                /       \                |         |         |
//               2---------4               1---------4---------7
//              / \       / \              |         |         |
//             /   \     /   \             |         |         |
//            /     \   /     \            |         |         |
//           /       \ /       \           0---------3---------6
//          0---------1---------3

   vert = 0;
   for (i = 0; i <= total_length; i++) {
      for (j = 0; j <= total_length; j++) {
         if (IsInteriorVert(base_vert, total_length, i, j)) {
            vert_index_i[vert] = i;
            vert_index_j[vert] = j;
            vert_index_sector[i][j] = vert++;
         }
         else vert_index_sector[i][j] = -1;
      };
   };

// Count all edges (even non-existing). Edges have a "type" whose meaning is:
// TAS: S (type 0), NE (type 1), NW (type 2)
// QAS: Horizontal (type 0), Vertical (type 1)
//
//                    -
//                   / \                   -----2---------5-----
//                  /   \                  |         |         |
//                 8     5                 7         9        11
//                /       \                |         |         |
//               -----2-----               |----1----+----4----|
//              / \       / \              |         |         |
//             /   \     /   \             6         8        10
//            6     3   7     4            |         |         |
//           /       \ /       \           -----0---------3-----
//          -----0---------1-----

   edge = 0;
   for (etype = 0; etype < cardinal_directions; etype++) {
      for (i = 0; i <= total_length; i++) {
         for (j = 0; j <= total_length; j++) {
            if (IsInteriorEdge(base_vert, total_length, etype, i, j)) {
               edge_index_i[edge] = i;
               edge_index_j[edge] = j;
               edge_index_sector[etype][i][j] = edge++;
            }
            else edge_index_sector[etype][i][j] = -1;
         };
      };
   };

// Count all faces (even non-existing).
//
//                    - 
//                   / \                   ---------------------
//                  /   \                  |         |         |
//                 /  3  \                 |    1    |    3    |
//                /       \                |         |         |
//               -----------               |---------+---------|
//              / \       / \              |         |         |
//             /   \  2  /   \             |    0    |    2    |
//            /  0  \   /  1  \            |         |         |
//           /       \ /       \           ---------------------
//          ---------------------

   face = 0;
   for (i = 0; i < total_length; i++) {
      for (j = 0; j <= square_fill * (total_length - 1); j++) {
         if (IsInteriorFace(base_vert, total_length, i, j)) {
            face_index_i[face] = i;
            face_index_j[face] = j;
            face_index_sector[i][j] = face++;
         }
         else face_index_sector[i][j] = -1;
      };
   };
};

/*!
\author Vladimir Florinski
\date 05/07/2024
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC void GeodesicSector<verts_per_face>::VertVertConn(void)
{
   int i, j, i1, j1, vert, iv;

   for (i = 0; i <= total_length; i++) {
      for (j = 0; j <= MaxVertJ(total_length, i); j++) {

// Clipped - add the NEXI mask
         vert = vert_index_sector[i][j];
         if (IsClippedVert(i, j)) {
            for (iv = 0; iv < edges_per_vert; iv++) vv_local[vert][iv] = -1;
            RAISE_BITS(vert_mask[vert], GEOELM_NEXI);
         }

// Neighbor list
         else {

// Set INTR, GHST, and BNDR masks
            if (BoundaryVertOfSector(i, j) != -1) RAISE_BITS(vert_mask[vert], GEOELM_BNDR);
            if (IsInteriorVertOfSector(i, j)) RAISE_BITS(vert_mask[vert], GEOELM_INTR);
            else RAISE_BITS(vert_mask[vert], GEOELM_GHST);

            for (iv = 0; iv < edges_per_vert; iv++) {
               i1 = i + vert_vert[iv][0];
               j1 = j + vert_vert[iv][1];
               if (IsClippedVert(i1, j1)) vv_local[vert][iv] = -1;
               else vv_local[vert][iv] = vert_index_sector[i1][j1];
            };
         };
      };
   };
};

/*!
\author Vladimir Florinski
\date 05/07/2024
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC void GeodesicSector<verts_per_face>::VertEdgeConn(void)
{
   int i, j, i1, j1, vert, ie, etype;

   for (i = 0; i <= total_length; i++) {
      for (j = 0; j <= MaxVertJ(total_length, i); j++) {

// Clipped
         vert = vert_index_sector[i][j];
         if (IsClippedVert(i, j)) {
            for (ie = 0; ie < edges_per_vert; ie++) ve_local[vert][ie] = -1;
         }

// Neighbor list
         else {
            for (ie = 0; ie < edges_per_vert; ie++) {
               i1 = i + vert_edge[ie][0];
               j1 = j + vert_edge[ie][1];
               etype  = vert_edge[ie][2];
               if (IsClippedEdge(etype, i1, j1)) ve_local[vert][ie] = -1;
               else ve_local[vert][ie] = edge_index_sector[etype][i1][j1];
            };
         };
      };
   };
};

/*!
\author Vladimir Florinski
\date 05/07/2024
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC void GeodesicSector<verts_per_face>::VertFaceConn(void)
{
   int i, j, i1, j1, vert, it;

   for (i = 0; i <= total_length; i++) {
      for (j = 0; j <= MaxVertJ(total_length, i); j++) {

// Clipped
         vert = vert_index_sector[i][j];
         if (IsClippedVert(i, j)) {
            for (it = 0; it < edges_per_vert; it++) vf_local[vert][it] = -1;
         }

// Neighbor list
         else {
            for (it = 0; it < edges_per_vert; it++) {
               i1 = i + vert_face[it][0];
               j1 = square_fill * j + vert_face[it][1];
               if (IsClippedFace(i1, j1)) vf_local[vert][it] = -1;
               else vf_local[vert][it] = face_index_sector[i1][j1];
            };
         };
      };
   };
};

/*!
\author Vladimir Florinski
\date 05/07/2024
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC void GeodesicSector<verts_per_face>::EdgeVertConn(void)
{
   int i, j, i1, j1, edge, iv, etype;

   for (etype = 0; etype < cardinal_directions; etype++) {
      for (i = 0; i <= total_length + edge_dimax[etype]; i++) {
         for (j = 0; j <= MaxVertJ(total_length, i) + edge_djmax[etype]; j++) {

// Clipped - add the NEXI mask
            edge = edge_index_sector[etype][i][j];
            if (IsClippedEdge(etype, i, j)) {
               for (iv = 0; iv < 2; iv++) ev_local[edge][iv] = -1;
               RAISE_BITS(edge_mask[edge], GEOELM_NEXI);
            }

// Neighbor list
            else {

// Set INTR, GHST, and BNDR masks
               if (BoundaryEdgeOfSector(etype, i, j) != -1) RAISE_BITS(edge_mask[edge], GEOELM_BNDR);
               if (IsInteriorEdgeOfSector(etype, i, j)) RAISE_BITS(edge_mask[edge], GEOELM_INTR);
               else RAISE_BITS(edge_mask[edge], GEOELM_GHST);

               for (iv = 0; iv < 2; iv++) {
                  i1 = i + edge_vert[etype][iv][0];
                  j1 = j + edge_vert[etype][iv][1];
                  ev_local[edge][iv] = vert_index_sector[i1][j1];
               };
            };
         };
      };
   };
};

/*!
\author Vladimir Florinski
\date 05/07/2024
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC void GeodesicSector<verts_per_face>::EdgeFaceConn(void)
{
   int i, j, i1, j1, edge, it, etype;

   for (etype = 0; etype < cardinal_directions; etype++) {
      for (i = 0; i <= total_length + edge_dimax[etype]; i++) {
         for (j = 0; j <= MaxVertJ(total_length, i) + edge_djmax[etype]; j++) {

// Clipped
            edge = edge_index_sector[etype][i][j];
            if (IsClippedEdge(etype, i, j)) {
               for (it = 0; it < 2; it++) ef_local[edge][it] = -1;
            }

// Neighbor list
            else {
               for (it = 0; it < 2; it++) {
                  i1 = i + edge_face[etype][it][0];
                  j1 = square_fill * j + edge_face[etype][it][1];
                  if (IsClippedFace(i1, j1)) ef_local[edge][it] = -1;
                  else ef_local[edge][it] = face_index_sector[i1][j1];
               };
            };
         };
      };
   };
};

/*!
\author Vladimir Florinski
\date 05/07/2024
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC void GeodesicSector<verts_per_face>::FaceVertConn(void)
{
   int i, j, i1, j1, face, iv;
   
   for (i = 0; i <= total_length - 1; i++) {
      for (j = 0; j <= MaxFaceJ(total_length, i); j++) {

// Clipped - add the NEXI mask
         face = face_index_sector[i][j];
         if (IsClippedFace(i, j)) {
            for (iv = 0; iv < verts_per_face; iv++) fv_local[face][iv] = -1;
            RAISE_BITS(face_mask[face], GEOELM_NEXI);
         }

// Neighbor list
         else {

// Set INTR and GHST masks
            if (IsInteriorFaceOfSector(i, j)) RAISE_BITS(face_mask[face], GEOELM_INTR);
            else RAISE_BITS(face_mask[face], GEOELM_GHST);

            for (iv = 0; iv < verts_per_face; iv++) {
               i1 = i + face_vert[j % square_fill][iv][0];
               j1 = j / square_fill + face_vert[j % square_fill][iv][1];
               fv_local[face][iv] = vert_index_sector[i1][j1];
            };
         };
      };
   };
};

/*!
\author Vladimir Florinski
\date 05/07/2024
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC void GeodesicSector<verts_per_face>::FaceEdgeConn(void)
{
   int i, j, i1, j1, face, ie, etype;

   for (i = 0; i <= total_length - 1; i++) {
      for (j = 0; j <= MaxFaceJ(total_length, i); j++) {

// Clipped
         face = face_index_sector[i][j];
         if (IsClippedFace(i, j)) {
            for (ie = 0; ie < verts_per_face; ie++) fe_local[face][ie] = -1;
         }

// Neighbor list
         else {
            for (ie = 0; ie < verts_per_face; ie++) {
               i1 = i + face_edge[j % square_fill][ie][0];
               j1 = j / square_fill + face_edge[j % square_fill][ie][1];
               etype = face_edge[j % square_fill][ie][2];
               fe_local[face][ie] = edge_index_sector[etype][i1][j1];
            };
         };
      };
   };
};

/*!
\author Vladimir Florinski
\date 05/07/2024
*/
template <int verts_per_face>
SPECTRUM_DEVICE_FUNC void GeodesicSector<verts_per_face>::FaceFaceConn(void)
{
   int i, j, i1, j1, face, it;
   
   for (i = 0; i <= total_length - 1; i++) {
      for (j = 0; j <= MaxFaceJ(total_length, i); j++) {

// Clipped
         face = face_index_sector[i][j];
         if (IsClippedFace(i, j)) {
            for (it = 0; it < verts_per_face; it++) ff_local[face][it] = -1;
         }

// Neighbor list
         else {
            for (it = 0; it < verts_per_face; it++) {
               i1 = i + face_face[j % square_fill][it][0];
               j1 = j + face_face[j % square_fill][it][1];
               if (IsClippedFace(i1, j1)) ff_local[face][it] = -1;
               else ff_local[face][it] = face_index_sector[i1][j1];
            };
         };
      };
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// GeodesicSector debug methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

#ifdef GEO_DEBUG

/*!
\author Vladimir Florinski
\date 01/10/2020
\param[in] type vertex (1), edge (2), or face (3)
*/
template <int verts_per_face>
void GeodesicSector<verts_per_face>::PrintAddresses(int type) const
{
   switch (type) {

   case 1:
      std::cout << "Printing vertex addresses\n";
      PrintConnectivity(total_length + 1, total_length + 1, 0, vert_index_sector);
      break;

   case 2:
      std::cout << "Printing edge addresses\n";
      for(auto etype = 0; etype < cardinal_directions; etype++) {
         PrintConnectivity(total_length + 1, total_length + 1, 0, edge_index_sector[etype]);
      };
      break;

   case 3:
      std::cout << "Printing face addresses\n";
      PrintConnectivity(total_length, square_fill * total_length - 1, 0, face_index_sector);
      break;
   };
};

/*!
\author Vladimir Florinski
\date 01/10/2020
\param[in] type The connectivity array to print

The type is one of the folloing: 1 is vertex-vertex, 2 is vertex-edge, 3 is vertex-face, 4 is edge-vertex, 6 is edge-face, 7 is face-vertex, 8 is face-edge, 9 is face-face.
*/
template <int verts_per_face>
void GeodesicSector<verts_per_face>::PrintConn(int type) const
{
   switch (type) {

   case 1:
      std::cout << "Printing vert-vert connectivity\n";
      PrintConnectivity(n_verts_withghost, edges_per_vert, 0, vv_local);
      break;

   case 2:
      std::cout << "Printing vert-edge connectivity\n";
      PrintConnectivity(n_verts_withghost, edges_per_vert, 0, ve_local);
      break;

   case 3:
      std::cout << "Printing vert-face connectivity\n";
      PrintConnectivity(n_verts_withghost, edges_per_vert, 0, vf_local);
      break;

   case 4:
      std::cout << "Printing edge-vert connectivity for block\n";
      PrintConnectivity(n_edges_withghost, 2, 0, ev_local);
      break;

   case 6:
      std::cout << "Printing edge-face connectivity for block\n";
      PrintConnectivity(n_edges_withghost, 2, 0, ef_local);
      break;

   case 7:
      std::cout << "Printing face-vert connectivity for block\n";
      PrintConnectivity(n_faces_withghost, verts_per_face, 0, fv_local);
      break;

   case 8:
      std::cout << "Printing face-edge connectivity for block\n";
      PrintConnectivity(n_faces_withghost, verts_per_face, 0, fe_local);
      break;

   case 9:
      std::cout << "Printing face-face connectivity for block\n";
      PrintConnectivity(n_faces_withghost, verts_per_face, 0, ff_local);
      break;
   };
};

/*!
\author Vladimir Florinski
\date 01/10/2020
\param[in] type Type of request: 1 is print vertex mask, 2 is print edge mask, 3 is print face mask
*/
template <int verts_per_face>
void GeodesicSector<verts_per_face>::PrintMask(int type) const
{
   uint16_t bitselect;

   switch (type) {

   case 1:
      std::cout << "Printing vertex mask\n";
      for (auto vert = 0; vert < n_verts_withghost; vert++) {
         std::cout << std::setw(6) << vert << " | ";
         bitselect = 1;
         while (bitselect) {
            std::cout << std::setw(2) << bool(vert_mask[vert] & bitselect);
            bitselect <<= 1;
         };
         std::cout << std::endl;
      };
      break;

   case 2:
      std::cout << "Printing edge mask\n";
      for (auto edge = 0; edge < n_edges_withghost; edge++) {
         std::cout << std::setw(6) << edge << " | ";
         bitselect = 1;
         while (bitselect) {
            std::cout << std::setw(2) << bool(edge_mask[edge] & bitselect);
            bitselect <<= 1;
         };
         std::cout << std::endl;
      };
      break;

   case 3:
      std::cout << "Printing face mask\n";
      for (auto face = 0; face < n_faces_withghost; face++) {
         std::cout << std::setw(6) << face << " | ";
         bitselect = 1;
         while (bitselect) {
            std::cout << std::setw(2) << bool(face_mask[face] & bitselect);
            bitselect <<= 1;
         };
         std::cout << std::endl;
      };
      break;
   };
};

#endif

template class GeodesicSector<3>;
template class GeodesicSector<4>;

};
