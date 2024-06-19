/*!
\file traversable_tesselation.hh
\brief Declares a tesselation with sector traversing capabilities
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAVERSABLE_TESSELATION_HH
#define SPECTRUM_TRAVERSABLE_TESSELATION_HH

#include "geodesic/requestable_tesselation.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// TraversableTesselation class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A SphericalTesselation derived class that can generate flat addressing for the elements suitable for sectors
\author Vladimir Florinski
*/
template <PolyType poly_type, int max_division>
class TraversableTesselation : public RequestableTesselation<poly_type, max_division>
{
   using RequestableTesselation<poly_type, max_division>::verts_per_face;
   using RequestableTesselation<poly_type, max_division>::edges_per_vert;
   using RequestableTesselation<poly_type, max_division>::vf_con;
   using RequestableTesselation<poly_type, max_division>::fv_con;
   using RequestableTesselation<poly_type, max_division>::ff_con;
   using RequestableTesselation<poly_type, max_division>::NVertNbrs;
   using RequestableTesselation<poly_type, max_division>::VertCC;
   using RequestableTesselation<poly_type, max_division>::IsInside;

protected:

//! Perform one of basic moves on the mesh.
   SPECTRUM_DEVICE_FUNC void Step(int div, int face1, int vert1, int dir, int& face2, int& vert2) const;

public:

//! Default constructor
   TraversableTesselation(void) = default;

//! Compute lists of faces and vertices in a sector with ghost cells.
   void GetAllInsideFaceNative(int divs, int sect, int divf, int nghost, int* flist, int* vlist, bool* corners) const;
};

/*!
\brief Partial specialization of TraversableTesselation for POLY_HEXAHEDRON
\author Vladimir Florinski
*/
template <int max_division>
class TraversableTesselation<POLY_HEXAHEDRON, max_division> : public RequestableTesselation<POLY_HEXAHEDRON, max_division>
{
   using RequestableTesselation<POLY_HEXAHEDRON, max_division>::verts_per_face;
   using RequestableTesselation<POLY_HEXAHEDRON, max_division>::edges_per_vert;
   using RequestableTesselation<POLY_HEXAHEDRON, max_division>::vf_con;
   using RequestableTesselation<POLY_HEXAHEDRON, max_division>::fv_con;
   using RequestableTesselation<POLY_HEXAHEDRON, max_division>::ff_con;
   using RequestableTesselation<POLY_HEXAHEDRON, max_division>::NVertNbrs;
   using RequestableTesselation<POLY_HEXAHEDRON, max_division>::VertCC;
   using RequestableTesselation<POLY_HEXAHEDRON, max_division>::IsInside;

protected:

//! Perform one of basic moves on the mesh.
   SPECTRUM_DEVICE_FUNC void Step(int div, int face1, int vert1, int dir, int& face2, int& vert2) const;

public:

//! Default constructor
   TraversableTesselation(void) = default;

//! Destructor
   ~TraversableTesselation() = default;

//! Compute lists of faces and vertices in a sector with ghost cells.
   void GetAllInsideFaceNative(int divs, int sect, int divf, int nghost, int* flist, int* vlist, bool* corners) const;
};

};

#endif
