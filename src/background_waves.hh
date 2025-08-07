/*!
\file background_waves.hh
\brief Declares a background consisting of a superposition of waves
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_WAVES_HH
#define SPECTRUM_BACKGROUND_WAVES_HH

#include "background_base.hh"
#include "common/matrix.hh"
#include "common/turb_prop.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundWaves class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the class
const std::string bg_name_waves = "BackgroundWaves";

/*!
\brief Electromagnetic fields of an ensemble of waves
\author Vladimir Florinski

Parameters: (BackgroundBase), [double kmin, double kmax, int n_waves, double variance, double slope] x n_turb_types
*/
template <typename Fields_>
class BackgroundWaves : public BackgroundBase<Fields_> {
public:

   using Fields = Fields_;
   using BackgroundBase = BackgroundBase<Fields>;
   using BackgroundBase::_status;
   using BackgroundBase::_fields;
   using BackgroundBase::_ddata;
   using BackgroundBase::_pos;
   using BackgroundBase::container;
   using BackgroundBase::r0;
   using BackgroundBase::B0;
   using BackgroundBase::dmax0;
   // methods
   using BackgroundBase::EvaluateBmag;
   using BackgroundBase::EvaluateDmax;
   using BackgroundBase::StopServerFront;
   using BackgroundBase::SetupBackground;
   using BackgroundBase::EvaluateBackground;
   using BackgroundBase::EvaluateBackgroundDerivatives;
   using BackgroundBase::NumericalDerivatives;

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

//! PSD for component "turb_alfven"
   void PSD_Alfven(void);

//! PSD for component "turb_transverse"
   void PSD_Transverse(void);

//! PSD for component "turb_longitudinal"
   void PSD_Longitudinal(void);

//! PSD for component "turb_isotropic"
   void PSD_Isotropic(void);

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct) override;

//! Compute the internal u, B, and E fields
   void EvaluateBackground(void) override;

//! Compute the internal u, B, and E derivatives
   void EvaluateBackgroundDerivatives(void) override;

//! Compute the maximum distance per time step
   void EvaluateDmax(void) override;

public:

//! Default constructor
   BackgroundWaves(void);

//! Copy constructor
   BackgroundWaves(const BackgroundWaves& other);

//! Destructor
   ~BackgroundWaves() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundWaves);
};

};

// Something like this is needed for templated classes
#include "background_waves.cc"

#endif
