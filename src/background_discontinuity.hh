/*!
\file background_discontinuity.hh
\brief Declares a simple planar MHD discontinuity field background
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a SHOCK coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_DISCONTINUITY_HH
#define SPECTRUM_BACKGROUND_DISCONTINUITY_HH

#include "background_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundDiscontinuity class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Planar MHD dicontinuity
\author Juan G Alonso Guzman

Parameters: (BackgroundBase), GeoVector n_discont, double v_discont, GeoVector u1, GeoVector B1
*/
template <typename HyperParams_>
class BackgroundDiscontinuity : public BackgroundBase<HyperParams_> {
private:

//! Readable name of the BackgroundDiscontinuity class
   static constexpr std::string_view bg_name = "BackgroundDiscontinuity";

public:

   using HyperParams = HyperParams_;
   using BackgroundBase = BackgroundBase<HyperParams>;
   using BackgroundBase::_status;
   using BackgroundBase::_fields;
   using BackgroundBase::_ddata;
   using BackgroundBase::_pos;
   using BackgroundBase::_t;
   using BackgroundBase::container;
   using BackgroundBase::u0;
   using BackgroundBase::r0;
   using BackgroundBase::B0;
   using BackgroundBase::dmax0;
   // methods
   using BackgroundBase::EvaluateBmag;
   using BackgroundBase::EvaluateDmax;
   using BackgroundBase::GetDmax;
   using BackgroundBase::StopServerFront;
   using BackgroundBase::SetupBackground;
//   using BackgroundBase::EvaluateBackground;
//   using BackgroundBase::EvaluateBackgroundDerivatives;
   using BackgroundBase::NumericalDerivatives;

protected:

//! Discontinuity normal (persistent)
   GeoVector n_discont;

//! Discontinuity velocity (persistent)
   double v_discont;

//! Downstream flow vector (persistent), "u0" is upstream flow vector
   GeoVector u1;

//! Downstream magnetic field (persistent), "B0" is upstream magnetic field
   GeoVector B1;

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct) override;

//! Compute the internal u, B, and E fields
   template <typename Fields>
   void EvaluateBackground(Fields&);

//! Compute the internal derivatives of the fields
   template <typename Fields>
   void EvaluateBackgroundDerivatives(Fields&);

public:

//! Default constructor
   BackgroundDiscontinuity(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundDiscontinuity(const std::string& name_in, uint16_t status_in);

//! Copy constructor
   BackgroundDiscontinuity(const BackgroundDiscontinuity& other);

//! Destructor
   ~BackgroundDiscontinuity() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundDiscontinuity);

};

};

// Something like this is needed for templated classes
#include "background_discontinuity.cc"

#endif
