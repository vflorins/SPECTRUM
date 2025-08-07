/*!
\file background_dipole.hh
\brief Declares a dipole magnetic field background without a flow
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_DIPOLE_HH
#define SPECTRUM_BACKGROUND_DIPOLE_HH

#include "background_base.hh"

namespace Spectrum {

//! Method for computing derivatives of B (0: analytical, 1: numerical)
#define DIPOLE_DERIVATIVE_METHOD 0

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundDipole class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the class
const std::string bg_name_dipole = "BackgroundDipole";

/*!
\brief Magnetic field of a dipole
\author Vladimir Florinski

Parameters: (BackgroundBase), double r_ref, double dmax_fraction
*/
template <typename Fields_>
class BackgroundDipole : public BackgroundBase<Fields_> {
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

//! Dipole moment (persistent)
   GeoVector M;

//! Maximum fraction of the radial distance per step (persistent)
   double dmax_fraction;

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
   BackgroundDipole(void);

//! Copy constructor
   BackgroundDipole(const BackgroundDipole& other);

//! Destructor
   ~BackgroundDipole() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundDipole);
};

};

// Something like this is needed for templated classes
#include "background_dipole.cc"

#endif
