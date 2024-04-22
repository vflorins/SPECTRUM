/*!
\file mapped_vertex.hh
\brief Class representing a single point geometric element
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_MAPPED_VERTEX
#define SPECTRUM_MAPPED_VERTEX

#ifdef GEO_DEBUG
#include <vector>
#include <iostream>
#include <iomanip>
#endif
#include "geometry/mapped_element.hh"

namespace Spectrum {

//! Number of anchor nodes per vertex
constexpr int n_vert_anchors = 1;

//! Reference coordinates of the vertex anchors
constexpr double vert_anchors[n_vert_anchors] = {0.0};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// MappedVert class
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A class describing a single point in 3D space
\author Vladimir Florinski
*/
class MappedVert : public MappedElement<n_vert_anchors>
{
public:

   using MappedElement<n_vert_anchors>::SetAnchors;

//! Default constructor
   SPECTRUM_DEVICE_FUNC MappedVert(void);

//! Constructor with arguments
   SPECTRUM_DEVICE_FUNC MappedVert(const GeoVector* new_anc);

//! Compute the position vector
   SPECTRUM_DEVICE_FUNC GeoVector Position(void) const;

#ifdef GEO_DEBUG

//! Tests if positions of the anchors are correct
   void TestAnchors(void) const;

#endif

};

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 04/02/2024
*/
SPECTRUM_DEVICE_FUNC inline MappedVert::MappedVert(void)
{
   ref_anchors = vert_anchors;
};

/*!
\author Vladimir Florinski
\date 08/08/2019
\param[in] new_anc New physical coordinates of the anchor
*/
SPECTRUM_DEVICE_FUNC inline MappedVert::MappedVert(const GeoVector* new_anc)
                                      : MappedVert()
{
   SetAnchors(new_anc);
};

/*!
\author Vladimir Florinski
\date 07/10/2019
\return Position (anchor point itself)
*/
SPECTRUM_DEVICE_FUNC inline GeoVector MappedVert::Position(void) const
{
   return anchors[0];
};

#ifdef GEO_DEBUG

/*!
\author Vladimir Florinski
\date 04/02/2024
*/
inline void MappedVert::TestAnchors(void) const
{
   GeoVector error;

   std::cerr << "-------------------------------------------------------------------------------------------\n";
   std::cerr << "Vertex anchor position test\n";
   std::cerr << "-------------------------------------------------------------------------------------------\n";
   std::cerr << " Anchor |                           Value                               |     Error\n";
   std::cerr << "-------------------------------------------------------------------------------------------\n";

   for(auto apt = 0; apt < n_vert_anchors; apt++) {
      error = Position() - anchors[apt];
      std::cerr << std::setw(5) << apt << "   |   "
                << std::setprecision(15)
                << std::fixed << std::setw(19) << anchors[apt][0]
                << std::fixed << std::setw(19) << anchors[apt][1]
                << std::fixed << std::setw(19) << anchors[apt][2] << "   |   "
                << std::fixed << std::setw(9) << std::setprecision(5) << error.Norm() / anchors[apt].Norm() * 100.0 << " %"
                << std::endl;
   };
   std::cerr << "-------------------------------------------------------------------------------------------\n";
};

#endif

};

#endif
