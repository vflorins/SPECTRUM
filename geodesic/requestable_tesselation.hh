/*!
\file requestable_tesselation.hh
\brief Declares a tesselation with additional helper functions
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_REQUESTABLE_TESSELATION_HH
#define SPECTRUM_REQUESTABLE_TESSELATION_HH

#include <geodesic/spherical_tesselation.hh>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RequestableTesselation class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A tesselation class with added functionality to find neighbors/children/parents and locate a position
\author Vladimir Florinski
*/
template <PolyType poly_type, int max_division>
class RequestableTesselation : public SphericalTesselation<poly_type, max_division>
{
protected:

   using SphericalTesselation<poly_type, max_division>::nverts;
   using SphericalTesselation<poly_type, max_division>::nfaces;
   using SphericalTesselation<poly_type, max_division>::verts_per_face;
   using SphericalTesselation<poly_type, max_division>::edges_per_vert;
   using SphericalTesselation<poly_type, max_division>::children_per_face;
   using SphericalTesselation<poly_type, max_division>::vert_cart;
   using SphericalTesselation<poly_type, max_division>::vf_con;
   using SphericalTesselation<poly_type, max_division>::ef_con;
   using SphericalTesselation<poly_type, max_division>::fv_con;
   using SphericalTesselation<poly_type, max_division>::fe_con;
   using SphericalTesselation<poly_type, max_division>::NVertNbrs;

//! Return the smallest possible division for a given face index.
   SPECTRUM_DEVICE_FUNC int GetMinDivision(int face) const;

//! Return the immediate parent of a given face.
   SPECTRUM_DEVICE_FUNC int GetParent(int div, int face) const;

//! Return all daughter faces of a given face.
   SPECTRUM_DEVICE_FUNC void GetChildren(int div, int face, int* children) const;

//! Check whether a face lies inside a sector
   SPECTRUM_DEVICE_FUNC bool IsInside(int divs, int sect, int divf, int face) const;

//! Check whether a vector from the origin lies in the interior of a face.
   SPECTRUM_DEVICE_FUNC bool VectorInsideFace(int div, int face, const GeoVector& v) const;

//! Recursively traverse a face tree and build a list of leaf nodes.
   SPECTRUM_DEVICE_FUNC void DescendTree(int divr, int root, int divl, int* list, int& index) const;

public:

//! Default constructor
   RequestableTesselation(void) = default;

//! Return the edge neighbors of a face
   int EdgeNeighborsOfFace(int div, int face, int* edges) const;

//! Return the vertex neighbors of a face
   int VertNeighborsOfFace(int div, int face, int* verts) const;

//! Return the face neighbors of an edge
   int FaceNeighborsOfEdge(int div, int edge, int* faces) const;

//! Return the face neighbors of a vertex
   int FaceNeighborsOfVert(int div, int vert, int* faces) const;

//! Return the vertices of a single t-face.
   void FaceVertCoords(int div, int face, GeoVector* v) const;

//! Find the face where a given vector lies.
   int Locate(int div, const GeoVector& v) const;

//! Generate a list of t-faces that lie inside a sector in tree format.
   void GetAllInsideFaceTree(int divs, int sect, int divf, int* list) const;

//! Fill provided arrays with vertex polar coordinates.
   void FillVertCoordArrays(int length, const int* list, GeoVector* vcart) const;
};

};

#endif
