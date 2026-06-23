/*!
\file background_waves.hh
\brief Declares a background consisting of a superposition of waves
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_WAVES_HH
#define SPECTRUM_BACKGROUND_WAVES_HH

#include "common/vectors.hh"
#include "common/matrix.hh"
#include "common/turb_prop.hh"
#include "common/random.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundWaves class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Electromagnetic fields of an ensemble of waves
\author Vladimir Florinski

Parameters: (BackgroundBase), [double kmin, double kmax, int n_waves, double variance, double slope] x n_turb_types
*/
template <typename HConfig_>
class BackgroundWaves {
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundWaves";

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;

// secular config:
   static constexpr bool requires_setup = true;
   static constexpr bool stochastic = true;

//   using BackgroundBase = BackgroundBase<HConfig>;
//   using BackgroundBase::_status;
//   using BackgroundBase::container;
//   using BackgroundBase::_ddata;
//   using BackgroundBase::dmax0;
//   using BackgroundBase::r0;
//   using BackgroundBase::u0;
//   using BackgroundBase::B0;
//   // this background uses rng
//   using BackgroundBase::rng;
//   // methods
//   using BackgroundBase::EvaluateDmax;
//   using BackgroundBase::GetDmax;
//   using BackgroundBase::StopServerFront;
//   using BackgroundBase::SetupBackground;

   static constexpr double dmax0 = Config::dmax0;

   static constexpr GeoVector r0 = Config::r0;

   static constexpr GeoVector B0 = Config::B0;

   static constexpr int n_turb_types = TurbProp::n_turb_types;

//! Random number generator object (persistent)
   std::shared_ptr<RNG> rng = nullptr;

protected:

//! Number of waves of each kind (persistent)
   int n_waves[n_turb_types];

//! Wave amplitudes (persistent)
   std::vector<double> Ampl[n_turb_types];

//! Wavenumbers
   std::vector<double> k[n_turb_types];

//! Cosines of the polarization angles (persistent)
   std::vector<double> cosa[n_turb_types];

//! Sines of the polarization angles (persistent)
   std::vector<double> sina[n_turb_types];

//! Phase angles (persistent)
   std::vector<double> phase[n_turb_types];

//! Basis vectors in the frame aligned with the wavevector (persistent)
   std::vector<GeoMatrix> basis[n_turb_types];

//! Shortest wave in the ensemble for time step (persistent)
   double shortest_wave;

////! PSD for component "turb_alfven"
//   static void PSD_Alfven(void);
//
////! PSD for component "turb_transverse"
//   static void PSD_Transverse(void);
//
////! PSD for component "turb_longitudinal"
//   static void PSD_Longitudinal(void);
//
////! PSD for component "turb_isotropic"
//   static void PSD_Isotropic(void);

//! Connect to an existing RNG object
   void ConnectRNG(const std::shared_ptr<RNG> rng_in) {rng = rng_in;};

//! Set up the field evaluator
   void SetupBackground(DataContainer& container_in);

public:

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
#include "background_waves.cc"

#endif
