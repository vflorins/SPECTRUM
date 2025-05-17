/*!
\file spherical_tesselation.cc
\brief Implements SphericalTesselation class, a recursive polygonal partitioning of a sphere
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include <geodesic/spherical_tesselation.hh>

namespace Spectrum {

/*!
\author Vladimir Florinski
\date 04/15/2020
\param[in]  n_nodes1    Total nodes of type 1
\param[in]  n_nbrs1     Number of type 2 neighbors per type 1 node
\param[in]  n_sing1     Number of singular nodes of type 1
\param[in]  n_nbrs1s    Number of type 2 neighbors per type 1 singular node
\param[in]  start_sing1 Starting index of singular nodes of type 1
\param[in]  conn_12     Forward connectivity list
\param[in]  n_nodes2    Total nodes of type 2
\param[in]  n_nbrs2     Number of type 1 neighbors per type 2 node
\param[in]  n_sing2     Number of singular nodes of type 2
\param[in]  n_nbrs2s    Number of type 1 neighbors per type 2 singular node
\param[in]  start_sing2 Starting index of singular nodes of type 2
\param[out] conn_21     Reverse connectivity list
\return Error code described in "tesselate.hh"
*/
TERR_TYPE BuildReverse(int n_nodes1, int n_nbrs1, int n_sing1, int n_nbrs1s, int start_sing1, const int* const* conn_12,
                       int n_nodes2, int n_nbrs2, int n_sing2, int n_nbrs2s, int start_sing2,       int* const* conn_21)
{
   int n_nbrs1_actual, n_nbrs2_actual, node1, node2, nbr;
   TERR_TYPE err = TESERR_NOERR;

// Create an array of entry points into the "conn_21" table.
   int* entry_point = new int[n_nodes2];
   memset(entry_point, 0x0, n_nodes2 * SZINT);

// Loop over type 1 nodes. Each type 1 node enters once in "conn_21" for each of its type 2 neighbor. The array "ventry_point" points to the first unfilled position.
   for (node1 = 0; node1 < n_nodes1; node1++) {
      n_nbrs1_actual = ((node1 >= start_sing1) && (node1 < start_sing1 + n_sing1) ? n_nbrs1s : n_nbrs1);
   
      for (nbr = 0; nbr < n_nbrs1_actual; nbr++) {
         node2 = conn_12[node1][nbr];

// Bad entry in the "conn_12" array encountered - skip to the next entry to avoid a memory error.
         if (node2 < 0 || node2 >= n_nodes2) {
            RAISE_BITS(err, TESERR_INDEX);
            continue;
         };

// Enter node1 and increment the entry point by one. An error flag is raised if some node2 has all its connections filled already.
         n_nbrs2_actual = ((node2 >= start_sing2) && (node2 < start_sing2 + n_sing2) ? n_nbrs2s : n_nbrs2);
         if (entry_point[node2] == n_nbrs2_actual) RAISE_BITS(err, TESERR_OVERF);
         else conn_21[node2][entry_point[node2]++] = node1;
      };
   };

// Check if some connections were left unfilled.
   for (node2 = 0; node2 < n_nodes2; node2++) {
      n_nbrs2_actual = ((node2 >= start_sing2) && (node2 < start_sing2 + n_sing2) ? n_nbrs2s : n_nbrs2);
      if (entry_point[node2] < n_nbrs2_actual) {
         RAISE_BITS(err, TESERR_UNDER);
         break;
      };
   };

// Complete the lists for singular "node2".
   for (node2 = start_sing2; node2 < start_sing2 + n_sing2; node2++) {
      for (nbr = n_nbrs2s; nbr < n_nbrs2; nbr++) conn_21[node2][nbr] = -1;
   };

   delete[] entry_point;
   return err;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Initialization SphericalTesselation methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/01/2024
*/
template <PolyType poly_type, int max_division>
SphericalTesselation<poly_type, max_division>::SphericalTesselation(void)
                                             : Polyhedron<poly_type>()
{
   AllocateStorage();
   ComputeAll();
};

/*!
\author Vladimir Florinski
\date 05/01/2024
*/
template <PolyType poly_type, int max_division>
SphericalTesselation<poly_type, max_division>::~SphericalTesselation()
{
   FreeStorage();
};

/*!
\author Vladimir Florinski
\date 05/01/2024
*/
template <PolyType poly_type, int max_division>
void SphericalTesselation<poly_type, max_division>::AllocateStorage(void)
{
// Number of elements at division 0
   verts_per_face[0] = poly_p[poly_type];
   edges_per_vert[0] = poly_q[poly_type];
   nverts[0] = Polyhedron<poly_type>::Nv;
   nedges[0] = Polyhedron<poly_type>::Ne;
   nfaces[0] = Polyhedron<poly_type>::Nf;
   children_per_face[0] = (poly_type == POLY_DODECAHEDRON ? 5 : 4);
   newverts_at_edge[0] = (poly_type == POLY_DODECAHEDRON ? false : true);
   newverts_at_face[0] = ((poly_type == POLY_HEXAHEDRON) || (poly_type == POLY_DODECAHEDRON) ? true : false);

// Number of elements at divisions > 0
   for (auto div = 1; div <= max_division; div++) {
      verts_per_face[div] = (poly_type == POLY_HEXAHEDRON ? 4 : 3);
      edges_per_vert[div] = 2 * verts_per_face[div] / (verts_per_face[div] - 2); // singular vertices have fewer

      children_per_face[div] = 4;
      newverts_at_edge[div] = true;
      newverts_at_face[div] = (poly_type == POLY_HEXAHEDRON ? true : false);

      nfaces[div] = nfaces[div - 1] * children_per_face[div - 1];
      nverts[div] = nverts[div - 1];
      if (newverts_at_edge[div - 1]) nverts[div] += nedges[div - 1];
      if (newverts_at_face[div - 1]) nverts[div] += nfaces[div - 1];
      nedges[div] = nverts[div] + nfaces[div] - 2;
   };

// Allocate memory for vertex coordinates. Vertices at lower divisions are subsets of vertices at higher divisions, so a single array is sufficient.
   vert_cart = new GeoVector[nverts[max_division]];

// Allocate memory for the connectivity arrays
   for (auto div = 0; div <= max_division; div++) {
      vv_con[div] = Create2D<int>(nverts[div], edges_per_vert[div]);
      ve_con[div] = Create2D<int>(nverts[div], edges_per_vert[div]);
      vf_con[div] = Create2D<int>(nverts[div], edges_per_vert[div]);

      ev_con[div] = Create2D<int>(nedges[div], 2);
      ef_con[div] = Create2D<int>(nedges[div], 2);

      fv_con[div] = Create2D<int>(nfaces[div], verts_per_face[div]);
      fe_con[div] = Create2D<int>(nfaces[div], verts_per_face[div]);
      ff_con[div] = Create2D<int>(nfaces[div], verts_per_face[div]);
   };

// Calculate Cartesian vertex coordinates at division 0
   for (auto vert = 0; vert < nverts[0]; vert++) {
      vert_cart[vert] = gv_nr;
      vert_cart[vert].ToCartesian(cos(Polyhedron<poly_type>::vlat[vert]), sin(Polyhedron<poly_type>::vlat[vert]),
                                  sin(Polyhedron<poly_type>::vlon[vert]), cos(Polyhedron<poly_type>::vlon[vert]));
   };

// Copy base polyhedron connectivity
   memcpy(vv_con[0][0], Polyhedron<poly_type>::vert_vert[0], nverts[0] * edges_per_vert[0] * sizeof(int));
   memcpy(fv_con[0][0], Polyhedron<poly_type>::face_vert[0], nfaces[0] * verts_per_face[0] * sizeof(int));
};

/*!
\author Vladimir Florinski
\date 08/30/2019
*/
template <PolyType poly_type, int max_division>
void SphericalTesselation<poly_type, max_division>::ComputeAll(void)
try {
   for (auto div = 0; div <= max_division; div++) {
      RefineVert(div);   // VV, div>0
      AddEdges(div);     // EV
      RefineFace(div);   // FV, div>0
      VertEdgeConn(div); // VE
      VertFaceConn(div); // VF
      EdgeFaceConn(div); // EF
      FaceEdgeConn(div); // FE
      FaceFaceConn(div); // FF
   };
}

catch(const TessError& err) {
   err.PrintErrors();
   return;
};

/*!
\author Vladimir Florinski
\date 07/23/2019
*/
template <PolyType poly_type, int max_division>
void SphericalTesselation<poly_type, max_division>::FreeStorage(void)
{
// Free memory used for vertex coordinates
   delete[] vert_cart;

// Release connectivity array memory
   for (auto div = 0; div <= max_division; div++) {
      Delete2D(vv_con[div]);
      Delete2D(ve_con[div]);
      Delete2D(vf_con[div]);
      Delete2D(ev_con[div]);
      Delete2D(ef_con[div]);
      Delete2D(fv_con[div]);
      Delete2D(fe_con[div]);
      Delete2D(ff_con[div]);
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Connectivity SphericalTesselation methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 04/13/2020
\param[in] div Division (daughter)
*/
template <PolyType poly_type, int max_division>
void SphericalTesselation<poly_type, max_division>::RefineVert(int div)
{
// May only be called for division 1 or higher.
   if (div == 0) return;

   static const std::string callerID = "RefineVert";
   if ((div < 1) || (div > max_division)) throw TessError(callerID, div, TESERR_INPUT, __LINE__);
   int vert0, vert1, vert2, vert3, vert_old, vert_new, edge, edge0, edge1, edge2, face, face0, face1;
   int iv, ie, ie0, ie1, ie2;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Compute the coordinates of the new vertices.
//----------------------------------------------------------------------------------------------------------------------------------------------------

// This variable will keep the index of the next vertex to be inserted
   vert_new = nverts[div - 1];

// Insert vertices at edge midpoints (triangle and quad parent face). The number of new vertices added is thus equal to the number of edges in the parent division.
   for (edge = 0; newverts_at_edge[div - 1] && (edge < nedges[div - 1]); edge++) {
      vert0 = ev_con[div - 1][edge][0];
      vert1 = ev_con[div - 1][edge][1];
      if ((vert0 < 0) || (vert0 >= nverts[div - 1]) || (vert1 < 0) || (vert1 >= nverts[div - 1])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

// Insert a new vertex between "vert0" and "vert1".
      vert_cart[vert_new] = Bisect(vert_cart[vert0], vert_cart[vert1]);
      vert_new++;
   };

// Insert vertices at face centers (quad and penta parent face). The number of new vertices added is thus equal to the number of faces in the parent division.
   for (face = 0; newverts_at_face[div - 1] && (face < nfaces[div - 1]); face++) {

// To get to this point a face must have either four of five vertices. When dividing a division 0 face we use the circumcenter defined with vertices 0,1,2. When dividing a division 1 and higher quad face we use the intersection of great circles through edge midpoints; great circles through pairs of vertices would produce poor results.
      if (div == 1) {
         vert0 = fv_con[div - 1][face][0];
         vert1 = fv_con[div - 1][face][1];
         vert2 = fv_con[div - 1][face][2];
         vert3 = 0;
      }
      else {
         vert0 = nverts[div - 1] + fe_con[div - 1][face][0];
         vert1 = nverts[div - 1] + fe_con[div - 1][face][1];
         vert2 = nverts[div - 1] + fe_con[div - 1][face][2];
         vert3 = nverts[div - 1] + fe_con[div - 1][face][3];
      };
      if ((vert0 < 0) || (vert0 >= nverts[div]) || (vert1 < 0) || (vert1 >= nverts[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);
      if ((vert2 < 0) || (vert2 >= nverts[div]) || (vert3 < 0) || (vert3 >= nverts[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);
      if (div == 1) vert_cart[vert_new] = PlaneNormal(vert_cart[vert0], vert_cart[vert1], vert_cart[vert2]);
      else vert_cart[vert_new] = GreatCircleInt(vert_cart[vert0], vert_cart[vert2], vert_cart[vert1], vert_cart[vert3]);
      vert_new++;
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Compute vertex neighbors of the new vertices.
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Compute neighbors of the new edge vertices (triangle and quad parent face).
   if (newverts_at_edge[div - 1]) {
      for (vert_new = nverts[div - 1]; vert_new < nverts[div - 1] + nedges[div - 1]; vert_new++) {

// Some handy indices for VV connectivity
         ie0 = 0;
         ie1 = NVertNbrs(div, vert_new);
         ie2 = ie1 / 2;

// Compute the parent's "edge0" that the vertex belongs to. That edge has two lower division vertices that are neighbors of the new vertex.
         edge0 = vert_new - nverts[div - 1];
         vert0 = ev_con[div - 1][edge0][0];
         vert1 = ev_con[div - 1][edge0][1];
         if((vert0 < 0) || (vert0 >= nverts[div - 1]) || (vert1 < 0) || (vert1 >= nverts[div - 1])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

         vv_con[div][vert_new][ie0] = vert0;
         vv_con[div][vert_new][ie2] = vert1;

// Find two lower division faces sharing "edge0".
         face0 = ef_con[div - 1][edge0][0];
         face1 = ef_con[div - 1][edge0][1];
         if ((face0 < 0) || (face0 >= nfaces[div - 1]) || (face1 < 0) || (face1 > nfaces[div - 1])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

// For quad meshes neighbors 1 and 3 are vertices at the centers of the adjacent parent faces.
         if (newverts_at_face[div - 1]) {
            vv_con[div][vert_new][    ie2 / 2] = nverts[div - 1] + nedges[div - 1] + face1;
            vv_con[div][vert_new][3 * ie2 / 2] = nverts[div - 1] + nedges[div - 1] + face0;
         }

// For triangle meshes we obtain two pairs of parent edges belonging to "face0" and "face1" that are C and CC from "edge0". The indices of these two edges will give us the indices of the new vertices at their midpoints.
         else {

// Neighbors 1 and 2 are at midpoints of edges of "face1".
            edge1 = EdgeCC(div - 1, face1, edge0, 1);
            edge2 = EdgeCC(div - 1, face1, edge0, 2);
            if ((edge1 < 0) || (edge1 >= nedges[div - 1]) || (edge2 < 0) || (edge2 >= nedges[div - 1])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

            vv_con[div][vert_new][ie0 + 1] = nverts[div - 1] + edge1;
            vv_con[div][vert_new][ie2 - 1] = nverts[div - 1] + edge2;

// Neighbors 4 and 5 are at midpoints of edges of "face0".
            edge1 = EdgeCC(div - 1, face0, edge0, 1);
            edge2 = EdgeCC(div - 1, face0, edge0, 2);
            if ((edge1 < 0) || (edge1 >= nedges[div - 1]) || (edge2 < 0) || (edge2 >= nedges[div - 1])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

            vv_con[div][vert_new][ie2 + 1] = nverts[div - 1] + edge1;
            vv_con[div][vert_new][ie1 - 1] = nverts[div - 1] + edge2;
         };
      };
   };

// Compute the neighbors of the new center vertices (quad and penta parent face). The following code only executes for quad parent faces since penta faces have no new vertices at edges.
   if (newverts_at_face[div - 1] && newverts_at_edge[div - 1]) {
      for (vert_new = nverts[div - 1] + nedges[div - 1]; vert_new < nverts[div]; vert_new++) {
         face0 = vert_new - nverts[div - 1] - nedges[div - 1];

         for (ie = 0; ie < verts_per_face[div - 1]; ie++) {
            edge0 = fe_con[div - 1][face0][ie];
            if ((edge0 < 0) || (edge0 >= nedges[div - 1])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

// The new vertex has the index of its parent edge.
            vv_con[div][vert_new][ie] = nverts[div - 1] + edge0;
         };
      };
   }

// This code only executes for penta parent faces.
   else if (newverts_at_face[div - 1]) {
      for (vert_new = nverts[div - 1]; vert_new < nverts[div]; vert_new++) {
         face0 = vert_new - nverts[div - 1];

         for (iv = 0; iv < verts_per_face[div - 1]; iv++) {
            vert0 = fv_con[div - 1][face0][iv];
            if ((vert0 < 0) || (vert0 >= nverts[div - 1])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

            vv_con[div][vert_new][iv] = vert0;
         };
      };
   };

// Pad the list for singular vertices.
   for (vert_new = nverts[div - 1]; vert_new < nverts[div]; vert_new++) {
      for (iv = NVertNbrs(div, vert_new); iv < edges_per_vert[div]; iv++) vv_con[div][vert_new][iv] = -1;
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Compute vertex neighbors of the parent vertices.
//----------------------------------------------------------------------------------------------------------------------------------------------------

   for (vert_old = 0; vert_old < nverts[div - 1]; vert_old++) {

// Triangle and quad parent faces
      if (newverts_at_edge[div - 1]) {
         for (ie = 0; ie < NVertNbrs(div - 1, vert_old); ie++) {
            edge0 = ve_con[div - 1][vert_old][ie];
            if ((edge0 < 0) || (edge0 >= nedges[div - 1])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

// The new vertex has the index of its parent edge.
            vv_con[div][vert_old][ie] = nverts[div - 1] + edge0;
         };
      }

// Penta parent faces - new vertex neighbors are inserted between the old vertex neighbors.
      else {
         for (iv = 0; iv < NVertNbrs(div, vert_old); iv++) {

// Old vertex neighbor
            if (!(iv % 2)) vv_con[div][vert_old][iv] = vv_con[div - 1][vert_old][iv / 2];

// Old face neighbor gives the index of the new vertex.
            else vv_con[div][vert_old][iv] = nverts[div - 1] + vf_con[div - 1][vert_old][(iv / 2 + 1) % NVertNbrs(div - 1, vert_old)];
         };
      };

// Pad the list for singular vertices.
      for (iv = NVertNbrs(div, vert_old); iv < edges_per_vert[div]; iv++) vv_con[div][vert_old][iv] = -1;
   };
};

/*!
\author Vladimir Florinski
\date 04/13/2020
\param[in] div Division
*/
template <PolyType poly_type, int max_division>
void SphericalTesselation<poly_type, max_division>::AddEdges(int div)
{
   static const std::string callerID = "AddEdges";
   if ((div < 0) || (div > max_division)) throw TessError(callerID, div, TESERR_INPUT, __LINE__);
   int vert1, vert2, edge = 0, iv;

// Loop over vertices. An edge is created for each _unique_ connection in VV. Because there are two connections between each pair of vertices ("vert1"->"vert2" and "vert2"->"vert1"), an edge is only created if "vert1" is greater than "vert2".
   for (vert1 = 0; vert1 < nverts[div]; vert1++) {
      for (iv = 0; iv < NVertNbrs(div, vert1); iv++) {
         vert2 = vv_con[div][vert1][iv];
         if ((vert2 < 0) || (vert2 >= nverts[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

// Insert an edge connecting "vert1" and "vert2" (once per vertex pair).
         if (vert2 > vert1) {

// Do we have too many edges?
            if (edge >= nedges[div]) throw TessError(callerID, div, TESERR_OVERF, __LINE__);
            ev_con[div][edge][0] = vert1;
            ev_con[div][edge][1] = vert2;
            edge++;
         };
      };
   };

// Not enough edges were created because the number of vertices or edges was wrong.
   if (edge < nedges[div]) throw TessError(callerID, div, TESERR_UNDER, __LINE__);
};

/*!
\author Vladimir Florinski
\date 04/13/2020
\param[in] div Division (daughter)
*/
template <PolyType poly_type, int max_division>
void SphericalTesselation<poly_type, max_division>::RefineFace(int div)
{
// May only be called for division 1 or higher.
   if (div == 0) return;

   static const std::string callerID = "RefineFace";
   if ((div < 1) || (div > max_division)) throw TessError(callerID, div, TESERR_INPUT, __LINE__);
   int vert1, vert2, edge1, edge2, face, face_new, ie;

// Each parent face is divided into four. One of the daughter faces (diagram) inherits the index of the parent.

/*
                      2
                     /0\                   3----------2----------2
                    /   \                  |0        3|1        0|
                   /  3  \                 |    3     |    2     |
                  /1     2\                |1        2|2        3|
                 2---------1               3----------+----------1
                /2\2     1/1\              |3        2|2        1|
               /   \0=par/   \             |  0=par   |    1     |
              /  1  \   /  2  \            |0        1|3        0|
             /0     1\0/2     0\           0----------0----------1
            0---------0---------1
*/

   for (face = 0; face < nfaces[div - 1]; face++) {

// Center face of a triangle parent face. New vertices are numbered after the parent edges they bisect. The CC ordering of the parent's FE array ensures CC ordering of the daughter's FV.
      if (newverts_at_edge[div - 1] && !newverts_at_face[div - 1]) {
         for (ie = 0; ie < verts_per_face[div]; ie++) {
            edge1 = fe_con[div - 1][face][ie];
            if ((edge1 < 0) || (edge1 >= nedges[div - 1])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

// The new vertex's index is computed from the parent edge index.
            fv_con[div][face][ie] = nverts[div - 1] + edge1;
         };
      };

// Corner faces - work with pairs of edges or vertices of the parent division.
      for (ie = 0; ie < verts_per_face[div - 1]; ie++) {
         edge1 = fe_con[div - 1][face][ie];
         edge2 = EdgeCC(div - 1, face, edge1, -1);
         if ((edge1 < 0) || (edge1 >= nedges[div - 1]) || (edge2 < 0) || (edge2 >= nedges[div - 1])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

// Corner vertex between "edge1" and "edge2"
         vert1 = fv_con[div - 1][face][ie];
         vert2 = VertCC(div - 1, face, vert1, 1);
         if ((vert1 < 0) || (vert1 >= nverts[div - 1]) || (vert2 < 0) || (vert2 >= nverts[div - 1])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

// Triangle parent faces
         if (newverts_at_edge[div - 1] && !newverts_at_face[div - 1]) {
            face_new = nfaces[div - 1] + (children_per_face[div - 1] - 1) * face + ie;

// Corner, edge, edge
            fv_con[div][face_new][0] = vert1;
            fv_con[div][face_new][1] = nverts[div - 1] + edge1;
            fv_con[div][face_new][2] = nverts[div - 1] + edge2;
         }

// Quad parent faces
         else if (newverts_at_edge[div - 1]) {
            if (!ie) face_new = face;
            else face_new = nfaces[div - 1] + (children_per_face[div - 1] - 1) * face + ie - 1;

// Corner, edge, center, edge
            fv_con[div][face_new][0] = vert1;
            fv_con[div][face_new][1] = nverts[div - 1] + edge1;
            fv_con[div][face_new][2] = nverts[div - 1] + nedges[div - 1] + face;
            fv_con[div][face_new][3] = nverts[div - 1] + edge2;
         }

// Penta parent faces
         else if (newverts_at_face[div - 1]) {
            if (!ie) face_new = face;
            else face_new = nfaces[div - 1] + (children_per_face[div - 1] - 1) * face + ie - 1;
         
// Corner, corner, center
            fv_con[div][face_new][0] = vert1;
            fv_con[div][face_new][1] = vert2;
            fv_con[div][face_new][2] = nverts[div - 1] + face;
         };
      };
   };
};

/*!
\author Vladimir Florinski
\date 05/02/2024
\param[in] div Division
*/
template <PolyType poly_type, int max_division>
void SphericalTesselation<poly_type, max_division>::VertEdgeConn(int div)
{
   static const std::string callerID = "VertEdgeConn";
   if ((div < 0) || (div > max_division)) throw TessError(callerID, div, TESERR_INPUT, __LINE__);
   int sing_vert, sing_nbrs, offset;
   TERR_TYPE err;

// Build VE table, which is the reverse of EV table (EV must be available).
   if ((poly_type == POLY_DODECAHEDRON) && div) {
      sing_vert = nfaces[0];
      sing_nbrs = verts_per_face[0];
      offset = nverts[0];
   }
   else {
      sing_vert = nverts[0];
      sing_nbrs = edges_per_vert[0];
      offset = 0;
   };
   err = BuildReverse(nedges[div], 2, 0, 2, 0, ev_con[div], nverts[div], edges_per_vert[div], sing_vert, sing_nbrs, offset, ve_con[div]);
   if (err != TESERR_NOERR) throw TessError(callerID, div, err, __LINE__);

// Make VE CC-ordered and synchronized with VV.
   int vert, vert0, vert1, edge, iv, ie;
   int vert_unsorted[edges_per_vert[div]];

   for (vert = 0; vert < nverts[div]; vert++) {

// For convenience, build the list of other vertices (different from "vert") of every neighbor edge.
      for (ie = 0; ie < NVertNbrs(div, vert); ie++) {
         edge = ve_con[div][vert][ie];
         if ((edge < 0) || (edge >= nedges[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

         vert0 = ev_con[div][edge][0];
         vert1 = ev_con[div][edge][1];
         if ((vert0 < 0) || (vert0 >= nverts[div]) || (vert1 < 0) || (vert1 >= nverts[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);
         if ((vert0 != vert) && (vert1 != vert)) throw TessError(callerID, div, TESERR_MISMT, __LINE__);
         
         vert_unsorted[ie] = (vert0 == vert ? vert1 : vert0);
      };

// Sort in ascending order. At the end of this loop VV and VE tables for "vert" are fully synchronized - see the diagram in "VertFaceConn()".
      iv = 0;
      while (iv < NVertNbrs(div, vert) - 1) {
         vert1 = vv_con[div][vert][iv];
         if ((vert1 < 0) || (vert1 >= nverts[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

// Search the unsorted portion of the list for "vert1".
         ie = InList(NVertNbrs(div, vert) - iv, vert_unsorted + iv, vert1);
         if (ie == -1) throw TessError(callerID, div, TESERR_MISMT, __LINE__);
         ie += iv;

// Swap the two entries.
         std::swap(vert_unsorted[iv], vert_unsorted[ie]);
         std::swap(ve_con[div][vert][iv], ve_con[div][vert][ie]);
         iv++;
      };
   };
};

/*!
\author Vladimir Florinski
\date 05/02/2024
\param[in] div Division
*/
template <PolyType poly_type, int max_division>
void SphericalTesselation<poly_type, max_division>::VertFaceConn(int div)
{
   static const std::string callerID = "VertFaceConn";
   if ((div < 0) || (div > max_division)) throw TessError(callerID, div, TESERR_INPUT, __LINE__);
   int sing_vert, sing_nbrs, offset;
   TERR_TYPE err;

// Build VF table, which is the reverse of FV table (FV must be available).
   if ((poly_type == POLY_DODECAHEDRON) && div) {
      sing_vert = nfaces[0];
      sing_nbrs = verts_per_face[0];
      offset = nverts[0];
   }
   else {
      sing_vert = nverts[0];
      sing_nbrs = edges_per_vert[0];
      offset = 0;
   };
   err = BuildReverse(nfaces[div], verts_per_face[div], 0, verts_per_face[div], 0, fv_con[div],
                      nverts[div], edges_per_vert[div], sing_vert, sing_nbrs, offset, vf_con[div]);
   if (err != TESERR_NOERR) throw TessError(callerID, div, err, __LINE__);

// Make VF CC-ordered and synchronized with VV, VE.
   int vert, vert1, face, iv, it;
   int vert_unsorted[edges_per_vert[div]];

   for (vert = 0; vert < nverts[div]; vert++) {

// For convenience, build the list of vertices _clockwise_ from "vert" in every neighbor face's FV. We rely on CC ordering of the FV table.
      for (it = 0; it < NVertNbrs(div, vert); it++) {
         face = vf_con[div][vert][it];
         if ((face < 0) || (face >= nfaces[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

         vert1 = VertCC(div, face, vert, -1);
         if ((vert1 < 0) || (vert1 >= nverts[div])) throw TessError(callerID, div, TESERR_MISMT, __LINE__);

         vert_unsorted[it] = vert1;
      };

// Sort in ascending order. At the end of this loop VV, VE, and VF tables for "vert" are fully synchronized as shown in the diagram below.

/*
                 4---------3
                / \       / \              ----------2----------
               /   \  4  /   \             |         |         |
              /  5  4   3  3  \            |    3    2    2    |
             /       \ /       \           |         |         |
            5----5----X----2----2          3----3----+----1----1
             \       / \       /           |         |         |
              \  0  0   1  2  /            |    0    0    1    |
               \   /  1  \   /             |         |         |
                \ /       \ /              ----------0----------
                 0---------1
*/

      iv = 0;
      while (iv < NVertNbrs(div, vert) - 1) {
         vert1 = vv_con[div][vert][iv];
         if ((vert1 < 0) || (vert1 >= nverts[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

// Search the unsorted portion of the list for "vert1".
         it = InList(NVertNbrs(div, vert) - iv, vert_unsorted + iv, vert1);
         if (it == -1) throw TessError(callerID, div, TESERR_MISMT, __LINE__);
         it += iv;

// Swap the two entries
         std::swap(vert_unsorted[iv], vert_unsorted[it]);
         std::swap(vf_con[div][vert][iv], vf_con[div][vert][it]);
         iv++;
      };
   };
};

/*!
\author Vladimir Florinski
\date 08/30/2019
\param[in] div The division
*/
template <PolyType poly_type, int max_division>
void SphericalTesselation<poly_type, max_division>::EdgeFaceConn(int div)
{
   static const std::string callerID = "EdgeFaceConn";
   if ((div < 0) || (div > max_division)) throw TessError(callerID, div, TESERR_INPUT, __LINE__);
   int vert0, vert1, edge, face1, face2;

// The neighbor face search is based on finding two faces that share both this edge's vertices.
   for (edge = 0; edge < nedges[div]; edge++) {
      vert0 = ev_con[div][edge][0];
      vert1 = ev_con[div][edge][1];
      if ((vert0 < 0) || (vert0 >= nverts[div]) || (vert1 < 0) || (vert1 >= nverts[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);
      
      Match2vf(div, vert0, vert1, face1, face2);
      if ((face1 < 0) || (face1 >= nfaces[div]) || (face2 < 0) || (face2 >= nfaces[div])) throw TessError(callerID, div, TESERR_MISMT, __LINE__);

// Synchronize EV and EF such that EF[0] is on the left and EF[1] on the right of the vector from EV[0] to EV[1] as shown in the diagram below. We rely on CC ordering of the FV table.
//
/*
                 -
                / \              -----------
               /   \             |         |
              /  0  \            |    0    |
             /       \           |         |
            0---------1          0---------1
             \       /           |         |
              \  1  /            |    1    |
               \   /             |         |
                \ /              -----------
                 -
*/

      if (VertCC(div, face1, vert0, 1) == vert1) {
         ef_con[div][edge][0] = face1;
         ef_con[div][edge][1] = face2;
      }
      else {
         ef_con[div][edge][0] = face2;
         ef_con[div][edge][1] = face1;
      };
   };
};

/*!
\author Vladimir Florinski
\date 04/13/2020
\param[in] div Division
*/
template <PolyType poly_type, int max_division>
void SphericalTesselation<poly_type, max_division>::FaceEdgeConn(int div)
{
   static const std::string callerID = "EdgeFaceConn";
   if ((div < 0) || (div > max_division)) throw TessError(callerID, div, TESERR_INPUT, __LINE__);
   int vert1, vert2, edge, face, ie;

// Find the edge connecting each pair of vertices. We rely on CC ordering of FV.
   for (face = 0; face < nfaces[div]; face++) {
      for (ie = 0; ie < verts_per_face[div]; ie++) {
         vert1 = fv_con[div][face][ie];
         vert2 = VertCC(div, face, vert1, 1);
         if ((vert1 < 0) || (vert1 >= nverts[div]) || (vert2 < 0) || (vert2 >= nverts[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

         edge = Match2ve(div, vert1, vert2);
         if ((edge < 0) || (edge >= nedges[div])) throw TessError(callerID, div, TESERR_MISMT, __LINE__);

         fe_con[div][face][ie] = edge;
      };
   };
};

/*!
\author Vladimir Florinski
\date 08/30/2019
\param[in] div Division
*/
template <PolyType poly_type, int max_division>
void SphericalTesselation<poly_type, max_division>::FaceFaceConn(int div)
{
   static const std::string callerID = "FaceFaceConn";
   if ((div < 0) || (div > max_division)) throw TessError(callerID, div, TESERR_INPUT, __LINE__);
   int edge, face, face0, face1, it;

// Use FE to find the edge, then pick the appropriate neighbor face. By the end of this loop FV, FE, anf FF tables are all synchronized as shown in the diagram below.
//
/*
                                                     -----------
            ----------2----------                    |         |
             \       / \       /                     |    2    |
              \  2  /   \  1  /                      |         |
               \   2     1   /             ----------3----2----2----------
                \ /       \ /              |         |         |         |
                 0----0----1               |    3    3         1    1    |
                  \       /                |         |         |         |
                   \  0  /                 ----------0----0----1----------
                    \   /                            |         |
                     \ /                             |    0    |
                      -                              |         |
                                                     -----------
*/

   for (face = 0; face < nfaces[div]; face++) {
      for (it = 0; it < verts_per_face[div]; it++) {
         edge = fe_con[div][face][it];
         if ((edge < 0) || (edge >= nedges[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);

         face0 = ef_con[div][edge][0];
         face1 = ef_con[div][edge][1];
         if ((face0 < 0) || (face0 >= nfaces[div]) || (face1 < 0) || (face1 >= nfaces[div])) throw TessError(callerID, div, TESERR_INDEX, __LINE__);
         if (face0 == face) ff_con[div][face][it] = face1;
         else if (face1 == face) ff_con[div][face][it] = face0;
         else throw TessError(callerID, div, TESERR_MISMT, __LINE__);
      };
   };
};

template class SphericalTesselation<POLY_TETRAHEDRON, 5>;
template class SphericalTesselation<POLY_HEXAHEDRON, 5>;
template class SphericalTesselation<POLY_OCTAHEDRON, 5>;
template class SphericalTesselation<POLY_DODECAHEDRON, 5>;
template class SphericalTesselation<POLY_ICOSAHEDRON, 5>;

template class SphericalTesselation<POLY_TETRAHEDRON, 6>;
template class SphericalTesselation<POLY_HEXAHEDRON, 6>;
template class SphericalTesselation<POLY_OCTAHEDRON, 6>;
template class SphericalTesselation<POLY_DODECAHEDRON, 6>;
template class SphericalTesselation<POLY_ICOSAHEDRON, 6>;

template class SphericalTesselation<POLY_TETRAHEDRON, 7>;
template class SphericalTesselation<POLY_HEXAHEDRON, 7>;
template class SphericalTesselation<POLY_OCTAHEDRON, 7>;
template class SphericalTesselation<POLY_DODECAHEDRON, 7>;
template class SphericalTesselation<POLY_ICOSAHEDRON, 7>;

};
