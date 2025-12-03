/*!
\file background_dipole.hh
\brief Declares a dipole magnetic field background without a flow
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_DIPOLE_HH
#define SPECTRUM_BACKGROUND_DIPOLE_HH

#include "common/vectors.hh"
#include "common/data_container.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundDipole class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Magnetic field of a dipole
\author Vladimir Florinski

Parameters: (BackgroundBase), double r_ref, double dmax_fraction
*/
template <typename HConfig_>
class BackgroundDipole {
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundDipole";

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;

// secular config:
   static constexpr bool requires_setup = true;
   static constexpr bool stochastic = false;

   // TODO: there appears to be no way to proceed based on SimpleArray's design using this "stateless" construction.
   //  To test, set 'requires_setup' to false, uncomment, and comment 'SetupBackground'.

//   static constexpr double dmax0 = Config::dmax0;
//
////! Maximum fraction of the radial distance per step (persistent)
//   static constexpr double dmax_fraction = Config::dmax_fraction;
//
//   static constexpr GeoVector r0 = Config::r0;
//
//   static constexpr double r_ref = Config::r_ref;
//
//   static constexpr GeoVector B0 = Config::B0;
//
////! Dipole moment (persistent)
//   static constexpr GeoVector M = B0*Cube(r_ref);

protected:

   double t0;

   GeoVector r0;

   GeoVector u0;

   double dmax0;

//! Maximum fraction of the radial distance per step (persistent)
   double dmax_fraction;

   double r_ref;

   GeoVector B0;

//! Dipole moment (persistent)
   GeoVector M;

public:

   //! Set up the field evaluator based on "params"
   void SetupBackground(DataContainer& container_in);

//! Compute the maximum distance per time step
   template <typename Coordinates>
   status_t EvaluateDmax(Coordinates&, double*);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   status_t EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   status_t EvaluateBackgroundDerivatives(Coordinates&, Fields&);

};

};

// Something like this is needed for templated classes
#include "background_dipole.cc"

#endif
