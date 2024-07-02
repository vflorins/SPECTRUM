/*!
\file drawable_tesselation.hh
\brief Extends the SphericalTesselation class with some visual and diagnostic methods
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DRAWABLE_TESSELATION_HH
#define SPECTRUM_DRAWABLE_TESSELATION_HH

#include "geodesic/spherical_tesselation.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DrawableTesselation class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief An extension of SphericalTesselation to generate figures and print debug information
\author Vladimir Florinski
*/
template <int poly_type, int max_division>
class DrawableTesselation : public SphericalTesselation<poly_type, max_division>
{
protected:

   using SphericalTesselation<poly_type, max_division>::nverts;
   using SphericalTesselation<poly_type, max_division>::nedges;
   using SphericalTesselation<poly_type, max_division>::nfaces;
   using SphericalTesselation<poly_type, max_division>::verts_per_face;
   using SphericalTesselation<poly_type, max_division>::edges_per_vert;
   using SphericalTesselation<poly_type, max_division>::vert_cart;
   using SphericalTesselation<poly_type, max_division>::vv_con;
   using SphericalTesselation<poly_type, max_division>::ve_con;
   using SphericalTesselation<poly_type, max_division>::vf_con;
   using SphericalTesselation<poly_type, max_division>::ev_con;
   using SphericalTesselation<poly_type, max_division>::ef_con;
   using SphericalTesselation<poly_type, max_division>::fv_con;
   using SphericalTesselation<poly_type, max_division>::fe_con;
   using SphericalTesselation<poly_type, max_division>::ff_con;
   using SphericalTesselation<poly_type, max_division>::NVertNbrs;

public:

//! Default constructor
   DrawableTesselation(void) = default;

//! Draw a projection of the grid onto the xz plane.
   void DrawGridArcs(int div, bool opaque, double rot_z, double rot_x, bool smooth) const;

//! Print the summary for the tesselation at every division.
   void PrintStats(void) const;

#ifdef GEO_DEBUG

//! Test the connectivity of the mesh.
   void TestConnectivity(int div) const;

//! Print a mesh connectivity array.
   void PrintConn(int div, int type) const;

#endif

};

};

#endif
