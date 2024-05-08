/*!
\file mapped_element.hh
\brief Base class representing a geometric element
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_MAPPED_ELEMENT
#define SPECTRUM_MAPPED_ELEMENT

#include "common/vectors.hh"

namespace Spectrum {

#define MAX_ELEMENT_ORDER 3

// This and derived files declare a set of classes to perform coordinate transformation from a reference element (point, line, or polygon) to a point, curved line, and curved polygon. The mapping is defined in terms of "anchor points", which are points on the elements where the mapping is exact (although the derivatives are not), and a set of Lagrangian basis functions associated with each point.

//----------------------------------------------------------------------------------------------------------------------------------------------------
// MappedElement class
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A class describing a generic mapped element in 3D space
\author Vladimir Florinski
*/
template <int n_anchors>
class MappedElement
{
protected:

//! Anchor points
   GeoVector anchors[n_anchors];

//! Reference coordinates of the anchors (use pointer to save on memory)
   const double* ref_anchors;

//! Default constructor
   SPECTRUM_DEVICE_FUNC MappedElement(void) = default;

public:

//! Return the number of anchors
   SPECTRUM_DEVICE_FUNC static int GetNAnchors(void);

//! Initialize the anchor set
   SPECTRUM_DEVICE_FUNC void SetAnchors(const GeoVector* new_anc);

//! Return the copy of the anchor set
   SPECTRUM_DEVICE_FUNC void GetAnchors(GeoVector* old_anc) const;

//! Return the physical coordinates of the anchors (as a read-only pointer)
   SPECTRUM_DEVICE_FUNC const GeoVector* GetAnchors(void) const;

//! Return the reference coordinates of the anchors (as a read-only pointer)
   SPECTRUM_DEVICE_FUNC const double* GetRefAnchors(void) const;
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\return The number of anchor points
*/
template <int n_anchors>
SPECTRUM_DEVICE_FUNC inline int MappedElement<n_anchors>::GetNAnchors(void)
{
   return n_anchors;
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[in] new_anc New physical coordinates of the anchors
*/
template <int n_anchors>
SPECTRUM_DEVICE_FUNC inline void MappedElement<n_anchors>::SetAnchors(const GeoVector* new_anc)
{
   memcpy(anchors, new_anc, n_anchors * sizeof(GeoVector));
};

/*!
\author Vladimir Florinski
\date 04/01/2024
\param[out] old_anc Physical coordinates of the anchors currently defined
*/
template <int n_anchors>
SPECTRUM_DEVICE_FUNC inline void MappedElement<n_anchors>::GetAnchors(GeoVector* old_anc) const
{
   memcpy(old_anc, anchors, n_anchors * sizeof(GeoVector));
};

/*!
\author Vladimir Florinski
\date 04/03/2024
return Physical coordinates of the anchors currently defined
*/
template <int n_anchors>
SPECTRUM_DEVICE_FUNC inline const GeoVector* MappedElement<n_anchors>::GetAnchors(void) const
{
   return anchors;
};

/*!
\author Vladimir Florinski
\date 04/03/2024
return Reference coordinates of the anchors
*/
template <int n_anchors>
SPECTRUM_DEVICE_FUNC inline const double* MappedElement<n_anchors>::GetRefAnchors(void) const
{
   return ref_anchors;
};

};

#endif
