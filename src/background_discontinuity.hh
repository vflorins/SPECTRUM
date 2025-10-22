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
template <typename HConfig_>
class BackgroundDiscontinuity : public BackgroundBase<HConfig_> {
private:

//! Readable name of the class
   static constexpr std::string_view bg_name = "BackgroundDiscontinuity";

public:

   using HConfig = HConfig_;
   using BackgroundConfig = Cond<std::same_as<typename HConfig::BackgroundConfig, Default>, BackgroundDefault<BackgroundDiscontinuity<HConfig>>, typename HConfig::BackgroundConfig>;
   using BackgroundCoordinates = BackgroundConfig::Coordinates;
   using BackgroundBase = BackgroundBase<HConfig>;
   using BackgroundBase::_status;
   using BackgroundBase::container;
   using BackgroundBase::_ddata;
   using BackgroundBase::dmax0;
   using BackgroundBase::r0;
   using BackgroundBase::u0;
   using BackgroundBase::B0;
   // methods
   using BackgroundBase::EvaluateAbsMag;
   using BackgroundBase::EvaluateDmax;
   using BackgroundBase::GetDmax;
   using BackgroundBase::StopServerFront;
   using BackgroundBase::SetupBackground;

   using BackgroundConfig::derivative_method;

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
   void SetupBackground(bool construct);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   void EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   void EvaluateBackgroundDerivatives(Coordinates&, Fields&);

public:

//! Default constructor
   BackgroundDiscontinuity(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundDiscontinuity(const std::string& name_in, uint16_t status_in);

//! Copy constructor
   BackgroundDiscontinuity(const BackgroundDiscontinuity& other);

//! Destructor
   ~BackgroundDiscontinuity() = default;

//! Clone function
   CloneFunctionBackground(BackgroundDiscontinuity);

};

};

// Something like this is needed for templated classes
#include "background_discontinuity.cc"

#endif
