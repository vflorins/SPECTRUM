/*!
\file background_discontinuity.hh
\brief Declares a simple planar MHD discontinuity field background
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a SHOCK coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_DISCONTINUITY_HH
#define SPECTRUM_BACKGROUND_DISCONTINUITY_HH

#include "common/vectors.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundDiscontinuity class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Planar MHD dicontinuity
\author Juan G Alonso Guzman

Parameters: (BackgroundBase), GeoVector n_discont, double v_discont, GeoVector u1, GeoVector B1
*/
template <typename HConfig_>
class BackgroundDiscontinuity {
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundDiscontinuity";

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;

// secular config:
   static constexpr bool requires_setup = false;
   static constexpr bool stochastic = false;

protected:

   static constexpr double dmax0 = Config::dmax0;

   static constexpr GeoVector r0 = Config::r0;

   //! Discontinuity normal (persistent)
   static constexpr GeoVector n_discont = Config::n_discont.Normalize();

//! Discontinuity velocity (persistent)
   static constexpr double v_discont = Config::v_discont;

//! Downstream flow vector (persistent), "u0" is upstream flow vector
   static constexpr GeoVector u0 = Config::u0;
   static constexpr GeoVector u1 = Config::u1;

//! Downstream magnetic field (persistent), "B0" is upstream magnetic field
   static constexpr GeoVector B0 = Config::B0;
   static constexpr GeoVector B1 = Config::B1;

public:

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackgroundDerivatives(Coordinates&, Fields&);

//! Compute the maximum distance per time step
   template <typename Coordinates>
   static status_t EvaluateDmax(Coordinates&, double*);


};

};

// Something like this is needed for templated classes
#include "background_discontinuity.cc"

#endif
