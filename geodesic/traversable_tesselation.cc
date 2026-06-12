/*!
\file traversable_tesselation.cc
\brief Implements a tesselation with sector traversing capabilities
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include <common/print_warn.hh>
#include <geodesic/traversable_tesselation.hh>

namespace Spectrum {

/*!
\author Vladimir Florinski
\date 05/03/2024
\param[in]  div   Division
\param[in]  face1 Starting t-face
\param[in]  vert1 Base vertex
\param[in]  dir   Direction of the move, with "vert1" in the lower left corner
\param[out] face2 Destination t-face (-1 for bad input)
\param[out] vert2 Base vertex of "face2" (-1 for bad input)

This function performs a basic move that allow visiting every face in the tesselation. The base vertex ("vert1") determines the base edge (CC from the base vertex). The moves are defined for the t-face resting on its base edge, so that the base vertex is in the SW corner. They are S (dir=1), NE (dir=2), and NW (dir=3).
*/
template <PolyType poly_type, int max_division>
SPECTRUM_DEVICE_FUNC void TraversableTesselation<poly_type, max_division>::Step(int div, int face1, int vert1, int dir, int& face2, int& vert2) const
{
   int iv1, iv2;
   if ((dir < 0) || (dir >= verts_per_face[div])) return;

// Find the order of "vert1" in FV of "face1".
   iv1 = InList(verts_per_face[div], fv_con[div][face1], vert1);

#ifdef GEO_DEBUG
   if (iv1 < 0) PrintError(__FILE__, __LINE__, "Incorrect base vertex", true);
#endif

/*
                 -
                / \
               /   \                  ----------2          ----------2
              /  1  \	             / \       /            \       / \
             /       \	            /   \  2  /              \  2  /   \
            1---------2            /  1  \   /                \   /  1  \
             \       /            /       \ /                  \ /       \
              \  2  /            1---------*                    1----------
               \   /
                \ /
                 -
*/

// Common to all cases "face2" is neighbor "iv1+dir" of "face1"
   face2 = ff_con[div][face1][(iv1 + dir) % verts_per_face[div]];
   
   switch (dir) {

   case 0:
      vert2 = fv_con[div][face1][(iv1 + 1) % verts_per_face[div]];
      break;

   case 1:
      vert2 = fv_con[div][face1][(iv1 + 1) % verts_per_face[div]];

// "vert2" is CC from "vert1". For triangular mesh (*) it is not the correct "vert2" yet.
      iv2 = InList(verts_per_face[div], fv_con[div][face2], vert2);
      vert2 = fv_con[div][face2][(iv2 + 1) % verts_per_face[div]];
      break;

   case 2:
      vert2 = fv_con[div][face1][(iv1 + 2) % verts_per_face[div]];
      break;

   default:
      break;
   };
};

/*!
\author Vladimir Florinski
\date 05/03/2024
\param[in]  div   Division
\param[in]  face1 Starting t-face
\param[in]  vert1 Base vertex
\param[in]  dir   Direction of the move, with "vert1" in the lower left corner
\param[out] face2 Destination t-face (-1 for bad input)
\param[out] vert2 Base vertex of "face2" (-1 for bad input)

This function performs a basic move that allow visiting every face in the tesselation. The base vertex ("vert1") determines the base edge (CC from the base vertex). The moves are defined for the t-face resting on its base edge, so that the base vertex is in the SW corner. They are S (dir=1) and E (dir=2).
*/
template <int max_division>
SPECTRUM_DEVICE_FUNC void TraversableTesselation<POLY_HEXAHEDRON, max_division>::Step(int div, int face1, int vert1, int dir, int& face2, int& vert2) const
{
   int iv1, iv2;
   if ((dir < 0) || (dir >= verts_per_face[div])) return;

// Find the order of "vert1" in FV of "face1".
   iv1 = InList(verts_per_face[div], fv_con[div][face1], vert1);

#ifdef GEO_DEBUG
   if (iv1 < 0) PrintError(__FILE__, __LINE__, "Incorrect base vertex", true);
#endif

/*
            -----------                                        -----------
            |	      |                                         |         |
            |	 1    |          ---------------------          |    2    |          ---------------------
            |	      |          |         |         |          |         |          |         |         |
            1----------         |    1    |    2    |          2----------          |    2    |    1    |
            |	      |          |         |         |          |         |          |         |         |
            |   2    |          1---------2----------          |    1    |          2---------1----------
            |	      |                                         |         |
            2----------                                        1----------
*/

// Common to all cases "face2" is neighbor "iv1+dir" of "face1"
   face2 = ff_con[div][face1][(iv1 + dir) % verts_per_face[div]];
   
   switch (dir) {

   case 0:

// Find the order of "vert1" in FV of "face2", then add 1.
      iv2 = InList(verts_per_face[div], fv_con[div][face2], vert1);
      vert2 = fv_con[div][face2][(iv2 + 1) % verts_per_face[div]];
      break;

   case 1:
      vert2 = fv_con[div][face1][(iv1 + 1) % verts_per_face[div]];
      break;

   case 2:
      vert2 = fv_con[div][face1][(iv1 + verts_per_face[div] - 1) % verts_per_face[div]];
      break;

    case 3:

// Find the order of "vert1" in FV of "face2", then add 3.
      iv2 = InList(verts_per_face[div], fv_con[div][face2], vert1);
      vert2 = fv_con[div][face2][(iv2 + verts_per_face[div] - 1) % verts_per_face[div]];
      break;

   default:
      break;
   };
};

/*!
\author Vladimir Florinski
\date 05/03/2024
\param[in]  divs    Division of the sector
\param[in]  sect    Sector
\param[in]  divf    Division of the faces in the mesh
\param[in]  nghost  Width of the ghost cell layer
\param[out] flist   TAS/QAS array of faces in the sector+ghost
\param[out] vlist   TAS/QAS array of vertices in the sector+ghost
\param[out] corners Corner type, true for singular corners

A crawler through a triangular region consisting of a sector surrounded by a layer of ghost faces on each side. The crawl is performed in the TAS pattern starting from the base (SW) corner, and the index of each t-face visited and each vertex encountered are recorded. NB: The TAS is the addressing systems used by grid blocks.
*/
template <PolyType poly_type, int max_division>
void TraversableTesselation<poly_type, max_division>::GetAllInsideFaceNative(int divs, int sect, int divf, int nghost,
                                                                             int* flist, int* vlist, bool* corners) const
{

// Pentagonal sectors are not allowed
   if ((poly_type == POLY_DODECAHEDRON) && !divs) {
      PrintError(__FILE__, __LINE__, "Division 0 sectors not allowed for DOCECAHEDRON tesselation", true);
      return;
   };

// TODO Implement triangular sectors with tetra- and tri-corners.
   if ((poly_type == POLY_TETRAHEDRON) || (poly_type == POLY_OCTAHEDRON)) {
      PrintError(__FILE__, __LINE__, "Not implemented yet for TETRAHEDRON and OCTAHEDRON tesselations", true);
      return;
   };

   int vert, face, iv, it, i, j, findex = 0, vindex = 0;
   
// Record singular corners.
   for (iv = 0; iv < verts_per_face[divs]; iv++) corners[iv] = (fv_con[divs][sect][iv] < nverts[0]);

// Find the interior SW corner vertex and t-face.
   int vert_start = fv_con[divs][sect][0];
   it = 0;
   while (!IsInside(divs, sect, divf, vf_con[divf][vert_start][it]) && it < edges_per_vert[divs]) it++;
   int face_start = vf_con[divf][vert_start][it];
   int side_length = Pow2(divf - divs);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Position the crawler in the SW corner
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Walk SW to (nghost,0).
   for (j = nghost; j > 0; j--) {
      Step(divf, face_start, vert_start, 0, face_start, vert_start);
      Step(divf, face_start, vert_start, 1, face_start, vert_start);
   };

// Record the SW unused triangle.
   for (i = 0; i < nghost; i++) {
      for (j = 0; j < i; j++) {
         vlist[vindex++] = -1;
         flist[findex++] = -1;
         flist[findex++] = -1;
      };
      vlist[vindex++] = -1;
      flist[findex++] = -1;
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Region 1: leftmost "nghost" diagonals
//----------------------------------------------------------------------------------------------------------------------------------------------------
/*
                               /\
                              /  \
                             /    \
                            /      \
                           /--------\
                          / \      / \
                         /   \    /   \
                        /     \  /     \
                       /       \/       \
                      /--------/\--------\
                     /\       /  \       /\
                    /..\     /    \     /  \
                   /....\   /      \   /    \
                  /......\ /        \ /      \
                 /........\----------\--------\
                /\.........\        / \       /\
               /  \.........\      /   \     /  \
              /    \.........\    /     \   /    \
             /      \.........\  /       \ /      \
            ----------------------------------------
*/

// We record three things per step: (1) the bottom-left vertex of the unshaded face, (2) the unshaded face, and (3) the neighbor shaded face. At the top of each diagonal, we record one extra bottom-left vertex and unshaded face.
   for (i = nghost; i < 2 * nghost; i++) {

// Start of the diagonal
      face = face_start;
      vert = vert_start;

// Record SW-S corner block, starting from an unshaded face, moving NW.
      for (j = 0; j < i - nghost; j++) {
         vlist[vindex++] = vert;
         flist[findex++] = face;
         Step(divf, face, vert, 2, face, vert);
         flist[findex++] = face;
         Step(divf, face, vert, 0, face, vert);
      };

// The last unshaded face of the SW-S corner block
      vlist[vindex++] = vert;
      flist[findex++] = face;

// Jump the cut or position at the start of the SW block. If the cut was crossed, adjust the base vertex.
      Step(divf, face, vert, 2, face, vert);
      if (corners[0]) vert = VertCC(divf, face, vert, -1);

// Record SW corner block, starting from a shaded face, moving NW.
      for (j = 0; j < 2 * nghost - i - 1; j++) {
         if (corners[0]) {
            flist[findex++] = -1;
            vlist[vindex++] = -1;
            flist[findex++] = -1;
         }
         else {
            flist[findex++] = face;
            Step(divf, face, vert, 0, face, vert);
            vlist[vindex++] = vert;
            flist[findex++] = face;
            Step(divf, face, vert, 2, face, vert);
         };
      };

// The last shaded face in the SW corner block.
      if (corners[0]) flist[findex++] = -1;
      else flist[findex++] = face;

// For non-singular "corners[0]", move N to the start of the SW-W block.
      if (!corners[0]) Step(divf, face, vert, 0, face, vert);

// Record SW-N corner block, starting from an unshaded face, moving NW.
      for (j = 0; j < i - nghost; j++) {
         vlist[vindex++] = vert;
         flist[findex++] = face;
         Step(divf, face, vert, 2, face, vert);
         flist[findex++] = face;
         Step(divf, face, vert, 0, face, vert);
      };

// The last unshaded face of this diagonal.
      vlist[vindex++] = vert;
      flist[findex++] = face;

// Move E to the next diagonal.
      Step(divf, face_start, vert_start, 1, face_start, vert_start);
      Step(divf, face_start, vert_start, 2, face_start, vert_start);
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Region 2: middle "side_length" diagonals
//----------------------------------------------------------------------------------------------------------------------------------------------------
/*
                               /\
                              /  \
                             /    \
                            /      \
                           /--------\
                          /.\      / \
                         /...\    /   \
                        /.....\  /     \
                       /.......\/       \
                      /.........\--------\
                     /\..........\       /\
                    /  \..........\     /  \
                   /    \..........\   /    \
                  /      \..........\ /      \
                 /--------\..........\--------\
                /\       / \..........\       /\
               /  \     /   \..........\     /  \
              /    \   /     \..........\   /    \
             /      \ /       \..........\ /      \
            ----------------------------------------
*/

// We record three things per step: (1) the bottom-left vertex of the unshaded face, (2) the unshaded face, and (3) the neighbor shaded face. At the top of each diagonal, we record one extra bottom-left vertex and unshaded face. After reaching the last diagonal we record an extra diagonal of bottom-right vertices. Finally, at the top of that extra diagonal the top vertex of the unshaded face is also recorded.
   for (i = 2 * nghost; i < 2 * nghost + side_length; i++) {

// Start of the diagonal
      vert = vert_start;
      face = face_start;

// Record large central trapezoidal block, starting from an unshaded face, moving NW.
      for (j = 0; j < i; j++) {
         vlist[vindex++] = vert;
         flist[findex++] = face;
         Step(divf, face, vert, 2, face, vert);
         flist[findex++] = face;
         Step(divf, face, vert, 0, face, vert);
      };
      vlist[vindex++] = vert;
      flist[findex++] = face;

// Move E to the next diagonal. Don't make the move when we reach the last diagonal because this may be a cut line.
      if (i != 2 * nghost + side_length - 1) {
         Step(divf, face_start, vert_start, 1, face_start, vert_start);
         Step(divf, face_start, vert_start, 2, face_start, vert_start);
      };
   };

// Repeat the last diagonal to record another line of vertices, moving NW. This must be done here because parts of this line will be inaccessible from Region 3, for singular "corners[1]" or "corners[2]".
   vert = vert_start;
   face = face_start;
   for (j = 0; j < 2 * nghost + side_length - 1; j++) {
      vlist[vindex++] = VertCC(divf, face, vert, 1);
      Step(divf, face, vert, 2, face, vert);
      Step(divf, face, vert, 0, face, vert);
   };
   vlist[vindex++] = VertCC(divf, face, vert, 1);
   vlist[vindex++] = VertCC(divf, face, vert, 2);

// For singular "corners[1]", take "nghost" steps NW to prepare for Region 3 sweep.
   for (j = 0; corners[1] && (j < nghost); j++) {
      Step(divf, face_start, vert_start, 2, face_start, vert_start);
      Step(divf, face_start, vert_start, 0, face_start, vert_start);
   };

// Take one step NE. This will take us either into the SE corner block (non-singular "corners[1]") or just above the corner vertex (singular "corners[1]"), in which case one more SE step is needed.
   Step(divf, face_start, vert_start, 1, face_start, vert_start);
   if (corners[1]) Step(divf, face_start, vert_start, 2, face_start, vert_start);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Region 3: rightmost nghost diagonals
//----------------------------------------------------------------------------------------------------------------------------------------------------
/*
                               /\
                              /  \
                             /    \
                            /      \
                           /--------\
                          / \........\
                         /   \........\
                        /     \........\
                       /       \........\
                      /--------/\........\
                     /\       /  \........\
                    /  \     /    \........\
                   /    \   /      \........\
                  /      \ /        \........\
                 /--------/----------\........\
                /\       / \        / \......./\
               /  \     /   \      /   \...../  \
              /    \   /     \    /     \.../    \
             /      \ /       \  /       \./      \
            ----------------------------------------
*/

// We record three things per step: (1) the top-right vertex of the shaded face, (2) the shaded face, and (3) the neighbor unshaded face (except in unused triangles).
   for (i = 2 * nghost + side_length; i < 3 * nghost + side_length; i++) {

// Record the SE unused triangle.
      for (j = 0; j < i - 2 * nghost - side_length; j++) {
         vlist[vindex++] = -1;
         flist[findex++] = -1;
         flist[findex++] = -1;
      };
      vlist[vindex++] = -1;
      flist[findex++] = -1;

// Start of the diagonal
      vert = vert_start;
      face = face_start;

// Record the SE corner block, starting from a shaded face, moving NW.
      for (j = 0; j < 3 * nghost + side_length - i - 1; j++) {
         if (corners[1]) {
            flist[findex++] = -1;
            vlist[vindex++] = -1;
            flist[findex++] = -1;
         }
         else {
            flist[findex++] = face;
            vlist[vindex++] = vert;
            Step(divf, face, vert, 0, face, vert);
            flist[findex++] = face;
            Step(divf, face, vert, 2, face, vert);
         };
      };
      if (corners[1]) flist[findex++] = -1;
      else {
         flist[findex++] = face;
         Step(divf, face, vert, 0, face, vert);
      };

// Record the NE trapezoidal block, starting from an unshaded face, moving NW.
      for (j = 0; j < i - 2 * nghost; j++) {
         if (!j) vlist[vindex++] = VertCC(divf, face, vert, 1);
         flist[findex++] = face;
         Step(divf, face, vert, 2, face, vert);
         vlist[vindex++] = vert;
         flist[findex++] = face;
         Step(divf, face, vert, 0, face, vert);
      };

// The last unshaded face and its top vertex. The top vertex must be recorded here in case we have a singular "corners[2]", when the diagonal ends here.
      flist[findex++] = face;
      vlist[vindex++] = VertCC(divf, face, vert, 2);

// Record the S corner block, starting from an unshaded face, moving NW.
      for (j = 0; j < 3 * nghost + side_length - i - 1; j++) {
         if (corners[2]) {
            flist[findex++] = -1;
            if (j) vlist[vindex++] = -1;
            flist[findex++] = -1;
         }
         else {
            Step(divf, face, vert, 2, face, vert);
            flist[findex++] = face;
            if (j) vlist[vindex++] = vert;
            Step(divf, face, vert, 0, face, vert);
            flist[findex++] = face;
         };
      };

// The last shaded face and its base vertex. If we are in he last diagonal, that vertex has been already entered, so we skip its recording here.
      if (corners[2]) {
         flist[findex++] = -1;
         if (i != 3 * nghost + side_length - 1) vlist[vindex++] = -1;
      }
      else {
         Step(divf, face, vert, 2, face, vert);
         flist[findex++] = face;
         if (i != 3 * nghost + side_length - 1) vlist[vindex++] = vert;
      };

// Record the N unused triangle.
      for (j = 0; j < i - 2 * nghost - side_length; j++) {
         vlist[vindex++] = -1;
         flist[findex++] = -1;
         flist[findex++] = -1;
      };
      vlist[vindex++] = -1;
      flist[findex++] = -1;

// Move E to the next diagonal, unless we are at the end of the block.
      if (i != 3 * nghost + side_length - 1) {
         if (corners[1]) {
            Step(divf, face_start, vert_start, 1, face_start, vert_start);
            Step(divf, face_start, vert_start, 2, face_start, vert_start);
         }
         else {
            Step(divf, face_start, vert_start, 0, face_start, vert_start);
            Step(divf, face_start, vert_start, 1, face_start, vert_start);
         };
      };
   };
};

/*!
\author Vladimir Florinski
\date 05/03/2024
\param[in]  divs    Division of the sector
\param[in]  sect    Sector
\param[in]  divf    Division of the faces in the mesh
\param[in]  nghost  Width of the ghost cell layer
\param[out] flist   TAS/QAS array of faces in the sector+ghost
\param[out] vlist   TAS/QAS array of vertices in the sector+ghost
\param[out] corners Corner type, true for singular corners

A crawler through a square region consisting of a sector surrounded by a layer of ghost faces on each side. The crawl is performed in the QAS pattern starting from the base (SW) corner, and the index of each t-face visited and each vertex encountered are recorded. NB: The QAS is the addressing systems used by grid blocks.
*/
template <int max_division>
void TraversableTesselation<POLY_HEXAHEDRON, max_division>::GetAllInsideFaceNative(int divs, int sect, int divf, int nghost,
                                                                                   int* flist, int* vlist, bool* corners) const
{
   int vert, face, iv, it, i, j, findex = 0, vindex = 0;
   
// Record singular corners.
   for (iv = 0; iv < verts_per_face[divs]; iv++) corners[iv] = (fv_con[divs][sect][iv] < nverts[0]);

// Find the interior SW corner vertex and t-face.
   int vert_start = fv_con[divs][sect][0];
   it = 0;
   while (!IsInside(divs, sect, divf, vf_con[divf][vert_start][it]) && it < edges_per_vert[divs]) it++;
   int face_start = vf_con[divf][vert_start][it];
   int side_length = Pow2(divf - divs);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Position the walker in the SW corner
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Walk W to (0, nghost).
   for (i = nghost; i > 0; i--) {
      Step(divf, face_start, vert_start, 3, face_start, vert_start);
   };

// Walk S to (0, 0) for non-singular "corners[0]".
   for (j = nghost; !corners[0] && (j > 0); j--) {
      Step(divf, face_start, vert_start, 0, face_start, vert_start);
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Region 1: leftmost "nghost" columns
//----------------------------------------------------------------------------------------------------------------------------------------------------
/*
//          --------------------------------
//          |........|            |        |
//          |........|            |        |
//          |........|            |        |
//          |........|------------+--------|
//          |........|            |        |
//          |........|            |        |
//          |........|            |        |
//          |........|            |        |
//          |........|            |        |
//          |........|------------+--------|
//          |........|            |        |
//          |........|            |        |
//          |........|            |        |
//          --------------------------------
*/

// We record two things per step: (1) the bottom-left vertex of the face and (2) the face. At the top of each columns, we also record the top-left vertex.
   for (i = 0; i < nghost; i++) {

// Bottom of the column
      vert = vert_start;
      face = face_start;

// Record SW corner block, moving N.
      for (j = 0; j < nghost; j++) {
         if (corners[0]) {
            vlist[vindex++] = -1;
            flist[findex++] = -1;
         }
         else {
            vlist[vindex++] = vert;
            flist[findex++] = face;
            if(j != nghost - 1) Step(divf, face, vert, 2, face, vert);
         };
      };

// Move N to the start of the W side block for non-singular "corners[0]". For singular "corners[0]", "vert_start" and "face_start" are already here.
      if (!corners[0]) Step(divf, face, vert, 2, face, vert);

// Record W side block, moving N.
      for (j = 0; j < side_length; j++) {
         vlist[vindex++] = vert;
         flist[findex++] = face;
         if (j != side_length - 1) Step(divf, face, vert, 2, face, vert);
      };

// Move N to the start of the NW corner block for non-singular "corners[3]".
      if (!corners[3]) Step(divf, face, vert, 2, face, vert);

// Record NW corner block, moving N. For singular "corners[3]" we must record the existing first top-left vertex.
      for (j = 0; j < nghost; j++) {
         if (corners[3]) {
            if (j == 0) vlist[vindex++] = VertCC(divf, face, vert, -1);
            else vlist[vindex++] = -1;
            flist[findex++] = -1;
         }
         else {
            vlist[vindex++] = vert;
            flist[findex++] = face;
            if (j != nghost - 1) Step(divf, face, vert, 2, face, vert);
         };
      };

// Record topmost vertex in the column.
      if (corners[3]) vlist[vindex++] = -1;
      else vlist[vindex++] = VertCC(divf, face, vert, -1);

// Move E to the next column.
      Step(divf, face_start, vert_start, 1, face_start, vert_start);
   };

// For singular "corners[0]" we must walk S to (nghost,0).
   for (j = nghost - 1; corners[0] && (j >= 0); j--) {
      Step(divf, face_start, vert_start, 0, face_start, vert_start);
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Region 2: middle "side_length" columns
//----------------------------------------------------------------------------------------------------------------------------------------------------
/*
//          --------------------------------
//          |        |............|        |
//          |        |............|        |
//          |        |............|        |
//          |--------|............|--------|
//          |        |............|        |
//          |        |............|        |
//          |        |............|        |
//          |        |............|        |
//          |        |............|        |
//          |--------|............|--------|
//          |        |............|        |
//          |        |............|        |
//          |        |............|        |
//          --------------------------------
*/

// We record two things per step: (1) the bottom-left vertex of the face and (2) the face. At the top of each columns, we also record the top-left vertex. After reaching the last column we record an extra column of bottom-right vertices. Finally, at the top of that extra column the top-right vertex is also recorded.
   for (i = nghost; i < nghost + side_length; i++) {

// Bottom of the column
      vert = vert_start;
      face = face_start;

// Record large central rectangular block, moving N.
      for (j = 0; j < 2 * nghost + side_length; j++) {
         vlist[vindex++] = vert;
         flist[findex++] = face;
         if (j != 2 * nghost + side_length - 1) Step(divf, face, vert, 2, face, vert);
      };

// Record topmost vertex in the column.
      iv = InList(verts_per_face[divf], fv_con[divf][face], vert);
      vlist[vindex++] = fv_con[divf][face][(iv + verts_per_face[divf] - 1) % verts_per_face[divf]];

// Move E to the next column. Don't make the move when we reach the last column because this may be a cut line.
      if (i != nghost + side_length - 1) Step(divf, face_start, vert_start, 1, face_start, vert_start);
   };

// Repeat the last column to record another line of vertices, moving N. This must be done here because parts of this line will be inaccessible from Region 3, for singular "corners[1]" or "corners[2]".
   vert = vert_start;
   face = face_start;
   for (j = 0; j < 2 * nghost + side_length; j++) {
      vlist[vindex++] = VertCC(divf, face, vert, 1);
      if (j != 2 * nghost + side_length - 1) Step(divf, face, vert, 2, face, vert);
   };
   vlist[vindex++] = VertCC(divf, face, vert, 2);

// Walk N to (nghost+side_length,nghost) for singular "corners[1]".
   for (j = nghost - 1; corners[0] && (j >= 0); j--) {
      Step(divf, face_start, vert_start, 2, face_start, vert_start);
   };

// Take one step E. We are now in the first column
   Step(divf, face_start, vert_start, 1, face_start, vert_start);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Region 3: rightmost nghost columns
//----------------------------------------------------------------------------------------------------------------------------------------------------
/*
//          --------------------------------
//          |        |            |........|
//          |        |            |........|
//          |        |            |........|
//          |--------+------------|........|
//          |        |            |........|
//          |        |            |........|
//          |        |            |........|
//          |        |            |........|
//          |        |            |........|
//          |--------+------------|........|
//          |        |            |........|
//          |        |            |........|
//          |        |            |........|
//          --------------------------------
*/

// We record two things per step: (1) the bottom-right vertex of the face and (2) the face. At the top of each columns, we also record the top-right vertex.
   for (i = nghost + side_length; i < 2 * nghost + side_length; i++) {

// Bottom of the column
      vert = vert_start;
      face = face_start;

// Record SE corner block, moving N.
      for (j = 0; j < nghost; j++) {
         if (corners[1]) {
            vlist[vindex++] = -1;
            flist[findex++] = -1;
         }
         else {
            vlist[vindex++] = VertCC(divf, face, vert, 1);
            flist[findex++] = face;
            if (j != nghost - 1) Step(divf, face, vert, 2, face, vert);
         };
      };

// Move N to the start of the E side block for non-singular "corners[1]". For singular "corners[1]", "vert_start" and "face_start" are already here.
      if (!corners[1]) Step(divf, face, vert, 2, face, vert);

// Record E side block, moving N.
      for (j = 0; j < side_length; j++) {
         vlist[vindex++] = VertCC(divf, face, vert, 1);
         flist[findex++] = face;
         if (j != side_length - 1) Step(divf, face, vert, 2, face, vert);
      };

// Move N to the start of the NE corner block for non-singular corners[2].
      if (!corners[2]) Step(divf, face, vert, 2, face, vert);

// Record NE corner block, moving N. For singular "corners[2]" we must record the existing first top-right vertex.
      for (j = 0; j < nghost; j++) {
         if (corners[2]) {
            if (j == 0) vlist[vindex++] = VertCC(divf, face, vert, 2);
            else vlist[vindex++] = -1;
            flist[findex++] = -1;
         }
         else {
            vlist[vindex++] = VertCC(divf, face, vert, 1);
            flist[findex++] = face;
            if (j != nghost - 1) Step(divf, face, vert, 2, face, vert);
         };
      };

// Record topmost vertex in the column.
      if (corners[2]) vlist[vindex++] = -1;
      else vlist[vindex++] = VertCC(divf, face, vert, 2);

// Move E to the next column, unless we are at the end of the block.
      if (i != 2 * nghost + side_length - 1) Step(divf, face_start, vert_start, 1, face_start, vert_start);
   };
};

template class TraversableTesselation<POLY_HEXAHEDRON, 4>;
template class TraversableTesselation<POLY_HEXAHEDRON, 5>;
template class TraversableTesselation<POLY_HEXAHEDRON, 6>;
template class TraversableTesselation<POLY_HEXAHEDRON, 7>;
template class TraversableTesselation<POLY_ICOSAHEDRON, 4>;
template class TraversableTesselation<POLY_ICOSAHEDRON, 5>;
template class TraversableTesselation<POLY_ICOSAHEDRON, 6>;
template class TraversableTesselation<POLY_ICOSAHEDRON, 7>;

};
