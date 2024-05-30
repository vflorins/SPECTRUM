/*!
\file tesselate_request.cc
\brief Implements a tesselation with additional helper functions
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "geodesic/requestable_tesselation.hh"

namespace Spectrum {

/*!
\author Vladimir Florinski
\date 08/30/2019
\param[in] face Face index
\return The lowest division that can have this face index
*/
template <PolyType poly_type, int max_division>
SPECTRUM_DEVICE_FUNC int RequestableTesselation<poly_type, max_division>::GetMinDivision(int face) const
{
   if((face < 0) || (face >= nfaces[max_division])) return -1;

   int div = 0;
   while(face >= nfaces[div]) div++;
   return div;
};

/*!
\author Vladimir Florinski
\date 04/06/2020
\param[in] div  Division
\param[in] face Face
\return The immediate parent of this face (-1 if face is division 0)
*/
template <PolyType poly_type, int max_division>
SPECTRUM_DEVICE_FUNC int RequestableTesselation<poly_type, max_division>::GetParent(int div, int face) const
{
   int mindiv = GetMinDivision(face);
   if((mindiv == -1) || (div < mindiv) || (div > max_division)) return -1;

// Parent face has the same index.
   if(div > mindiv) return face;

// "mindiv" and "div" are the same.
   else {

// No parent if face is division 0.
      if(!div) return -1;

// General formula that takes into account that one of the children inherits the parent's index.
      else return (face - nfaces[div - 1]) / (children_per_face[div - 1] - 1);
   };
};

/*!
\author Vladimir Florinski
\date 04/06/2020
\param[in]  div      Division
\param[in]  face     Face
\param[out] children Daughter faces in an array
*/
template <PolyType poly_type, int max_division>
SPECTRUM_DEVICE_FUNC void RequestableTesselation<poly_type, max_division>::GetChildren(int div, int face, int* children) const
{
   int mindiv = GetMinDivision(face);
   if((mindiv == -1) || (div < mindiv) || (div > max_division)) {
      memset(children, -1, children_per_face[div] * SZINT);
      return;
   };

// General formula that takes into account that one of the children inherits the parent's index.
   children[0] = face;
   for(auto it = 1; it < children_per_face[div]; it++) {
      children[it] = nfaces[div] + face * (children_per_face[div] - 1) + it - 1;
   };
};

/*!
\author Vladimir Florinski
\date 08/30/2019
\param[in] divs Division of the sector
\param[in] sect Sector
\param[in] divf Division of the t-face
\param[in] face Face
\return True if "face" is inside this sector
*/
template <PolyType poly_type, int max_division>
SPECTRUM_DEVICE_FUNC bool RequestableTesselation<poly_type, max_division>::IsInside(int divs, int sect, int divf, int face) const
{
   if(face == -1) return false;

// The sector must be the (grand)parent of the face. This also covers the trivial case of divf=divs.
   int parent = face;
   int div = divf;
   while(div > divs) parent = GetParent(div--, parent);
   return (parent == sect);
};

/*!
\author Vladimir Florinski
\date 04/06/2020
\param[in] div  Division
\param[in] face Face
\param[in] v    Vector from the origin to the point
\return True if the point lies inside, false otherwise
*/
template <PolyType poly_type, int max_division>
SPECTRUM_DEVICE_FUNC bool RequestableTesselation<poly_type, max_division>::VectorInsideFace(int div, int face, const GeoVector& v) const
{
// In the interest of speed we do not check if "face" and "div" are valid.
   GeoVector edge_norm;

// Check the sign of the projection of "v" onto planes through each vector pair. This relies on counter-clockwise ordering of the vertices.
   for(auto iv = 0; iv < verts_per_face[div]; iv++) {
      edge_norm = vert_cart[fv_con[div][face][iv]] ^ vert_cart[fv_con[div][face][(iv + 1) % verts_per_face[div]]];
      if(v * edge_norm < 0.0) return false;
   };
   return true;
};

/*!
\author Vladimir Florinski
\date 04/06/2020
\param[in]     divr  Starting (root) level in the tree
\param[in]     root  Current root
\param[in]     divl  Ending (leaf) level in the tree
\param[in,out] list  An array of nodes at level "divl"
\param[in,out] index Current position in the "list"
\note This is a recursive function
*/
template <PolyType poly_type, int max_division>
SPECTRUM_DEVICE_FUNC void RequestableTesselation<poly_type, max_division>::DescendTree(int divr, int root, int divl, int* list, int& index) const
{
   int br, branches[5];

// Reached the leaf level - record face in the list and advance the index.
   if(divr == divl) list[index++] = root;

// Go to the next branch level of the tree - this function is intended to be called recursively.
   else {
      GetChildren(divr, root, branches);
      for(br = 0; br < children_per_face[divr]; br++) {
         DescendTree(divr + 1, branches[br], divl, list, index);
      };
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 08/30/2019
\param[in]  div      Division
\param[in]  face     Face
\param[out] edges    Edges of this face
\param[out] ef       EF connectivity for each edge
\param[out] vertices Vertices of this face
\param[out] vf       VF connectivity for each vertex
*/
template <PolyType poly_type, int max_division>
void RequestableTesselation<poly_type, max_division>::ExchangeSites(int div, int face, int* edges, int* const* ef, int* vertices, int* const* vf) const
{
   int edge, vert, ie, iv;

// Copy the EF array
   for(ie = 0; ie < verts_per_face[div]; ie++) {
      edge = fe_con[div][face][ie];
      edges[ie] = edge;
      memcpy(ef[ie], ef_con[div][edge], 2 * SZINT);
   };

// Copy the VF array, irrespective of singular vertices
   for(iv = 0; iv < verts_per_face[div]; iv++) {
      vert = fv_con[div][face][iv];
      vertices[iv] = vert;
      memcpy(vf[iv], vf_con[div][vert], edges_per_vert[div] * SZINT);
   };
};

/*!
\author Vladimir Florinski
\date 08/30/2019
\param[in]  div  Division
\param[in]  face Face
\param[out] v    Cartesian coordinates of the vertices
*/
template <PolyType poly_type, int max_division>
void RequestableTesselation<poly_type, max_division>::FaceVertCoords(int div, int face, GeoVector* v) const
{
   for(auto iv = 0; iv < verts_per_face[div]; iv++) v[iv] = vert_cart[fv_con[div][face][iv]];
};

/*!
\author Vladimir Florinski
\date 08/30/2019
\param[in] div Division
\param[in] v   Vector from the origin to the point
\return Index of the t-face where this point lies
*/
template <PolyType poly_type, int max_division>
int RequestableTesselation<poly_type, max_division>::Locate(int div, const GeoVector& v) const
{
   double sp, spmax = 0.0;
   int vert_near = 1, vert, face, div1, it, children[4];

// Find the largest (division 0) face. First, narrow down the search by computing the largest scalar product.
   for(vert = 0; vert < nverts[0]; vert++) {
      sp = v * vert_cart[vert];
      if(sp > spmax) {
         spmax = sp;
         vert_near = vert;
      };
   };

// Search among the neighbor faces of "vert_near" using "VectorInsideFace()"
   it = 0;
   while(!VectorInsideFace(0, vf_con[0][vert_near][it], v) && it < edges_per_vert[div] - 1) it++;
   face = vf_con[0][vert_near][it];

// Recursively search among the 4 daughter faces. For a level 0 request, this loop is skipped.
   for(div1 = 1; div1 <= div; div1++) {
      GetChildren(div1 - 1, face, children);
      it = 0;
      while(!VectorInsideFace(div1, children[it], v) && it < 4) it++;

// This is now the parent t-face for the next iteration
      face = children[it];
   };

   return face;
};

/*!
\author Vladimir Florinski
\date 08/30/2019
\param[in]  divs Division of the sector
\param[in]  sect Sector
\param[in]  divf Division of the faces in the mesh
\param[out] list An array of faces at divf inside the sector in tree format
*/
template <PolyType poly_type, int max_division>
void RequestableTesselation<poly_type, max_division>::GetAllInsideFaceTree(int divs, int sect, int divf, int* list) const
{
   int idx = 0;

// Zeroth element is the sector itself. The rest of the work is done by "DescendTree()".
   list[idx++] = sect;
   DescendTree(divs, sect, divf, list, idx);
};

/*!
\author Vladimir Florinski
\date 08/30/2019
\param[in]  length Length of the list
\param[in]  list   The list of indices
\param[out] vcart  Array of vertex coordinates
*/
template <PolyType poly_type, int max_division>
void RequestableTesselation<poly_type, max_division>::FillVertCoordArrays(int length, const int* list, GeoVector* vcart) const
{
   int vert, idx;

   for(idx = 0; idx < length; idx++) {
      vert = list[idx];
      if((vert >= 0) && (vert < nverts[max_division])) vcart[idx] = vert_cart[vert];
      else vcart[idx] = gv_zeros;
   };
};

template class RequestableTesselation<POLY_HEXAHEDRON, 4>;
template class RequestableTesselation<POLY_HEXAHEDRON, 5>;
template class RequestableTesselation<POLY_HEXAHEDRON, 6>;
template class RequestableTesselation<POLY_ICOSAHEDRON, 4>;
template class RequestableTesselation<POLY_ICOSAHEDRON, 5>;
template class RequestableTesselation<POLY_ICOSAHEDRON, 6>;

};
