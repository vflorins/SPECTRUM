/*!
\file background_shock.hh
\brief Declares a simple planar shock field background
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a SHOCK coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SHOCK_HH
#define SPECTRUM_BACKGROUND_SHOCK_HH

#include "common/vectors.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundShock class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Planar MHD shock
\author Juan G Alonso Guzman

Parameters: (BackgroundBase), GeoVector n_shock, double v_shock, double compression
*/
template <typename HConfig_>
class BackgroundShock {
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundShock";

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;

// secular config:
   static constexpr bool requires_setup = false;
   static constexpr bool stochastic = false;

private:

   static constexpr GeoVector compute_u1() {
      auto u0t = u0 - (u0 * n_shock) * n_shock;
      return u0n / compression + v_shock * n_shock + u0t;
   }

   static constexpr GeoVector compute_B1() {
      auto B0n = (B0 * n_shock) * n_shock;
      auto B0t = B0 - B0n;
      return B0n + B0t * compression;
   }

   static constexpr GeoVector compute_n_shock() {
      auto ns = Config::n_shock;
      return ns.Normalize();
   }

protected:

   static constexpr double dmax0 = Config::dmax0;

   static constexpr GeoVector r0 = Config::r0;

   static constexpr GeoVector B0 = Config::B0;

   static constexpr GeoVector u0 = Config::u0;

   //! Shock normal (persistent)
   static constexpr GeoVector n_shock = compute_n_shock();

   //! Shock velocity (persistent)
   static constexpr double v_shock = Config::v_shock;

   static constexpr GeoVector u0n = (u0 * n_shock - v_shock) * n_shock;

   //! Compression ratio (persistent)
   static constexpr double compression = Config::compression;

//! Downstream velocity (persistent), "u0" is upstream flow vector
   static constexpr GeoVector u1 = compute_u1();

//! Downstream magnetic field (persistent), "B0" is upstream magnetic field
   static constexpr GeoVector B1 = compute_B1();

   static_assert(u0n * n_shock <= 0.0, "Upstream normal velocity is in the wrong direction");

public:

//! Compute the maximum distance per time step
   template <typename Coordinates>
   static status_t EvaluateDmax(Coordinates&, double*);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackgroundDerivatives(Coordinates&, Fields&);

};

};

// Something like this is needed for templated classes
#include "background_shock.cc"

#endif
