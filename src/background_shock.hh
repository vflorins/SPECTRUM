/*!
\file background_shock.hh
\brief Declares a simple planar shock field background
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a SHOCK coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SHOCK_HH
#define SPECTRUM_BACKGROUND_SHOCK_HH

#include "background_base.hh"

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
class BackgroundShock : public BackgroundBase<HConfig_> {
private:

//! Readable name of the class
   static constexpr std::string_view bg_name = "BackgroundShock";

public:

   using HConfig = HConfig_;
   using BackgroundConfig = Cond<std::same_as<typename HConfig::BackgroundConfig, Default>, BackgroundDefault<BackgroundShock<HConfig>>, typename HConfig::BackgroundConfig>;
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

//! Shock normal (persistent)
   GeoVector n_shock;

//! Shock velocity (persistent)
   double v_shock;

//! Compression ratio (persistent)
   double compression;

//! Downstream velocity (persistent), "u0" is upstream flow vector
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
   BackgroundShock(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundShock(const std::string& name_in, uint16_t status_in);

//! Copy constructor
   BackgroundShock(const BackgroundShock& other);

//! Destructor
   ~BackgroundShock() = default;

//! Clone function
   CloneFunctionBackground(BackgroundShock);

};

};

// Something like this is needed for templated classes
#include "background_shock.cc"

#endif
