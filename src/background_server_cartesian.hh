/*!
\file background_server_cartesian.hh
\brief Declares a background class using data from uniform Cartesian grid on distributed memory
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SERVER_CARTESIAN_HH
#define SPECTRUM_BACKGROUND_SERVER_CARTESIAN_HH

#include "background_server.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundServerCartesian class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma background interface to a uniform Cartesian grid server
\author Juan G Alonso Guzman

Parameters: (BackgroundServer)
*/
template <typename HyperParams_>
class BackgroundServerCartesian : public BackgroundServer<HyperParams_> {
private:

   //! Readable name of the class
   static constexpr std::string_view bg_name = "BackgroundServerCartesian";

public:

   using HyperParams = HyperParams_;
   using BackgroundBase = BackgroundBase<HyperParams>;
   using BackgroundServer = BackgroundServer<HyperParams>;
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
   using BackgroundBase::GetDmax;
   using BackgroundBase::StopServerFront;
   using BackgroundBase::SetupBackground;
//   using BackgroundBase::EvaluateBackground;
//   using BackgroundBase::EvaluateBackgroundDerivatives;
   using BackgroundBase::NumericalDerivatives;

public:

//! Default constructor
   BackgroundServerCartesian(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundServerCartesian(const std::string& name_in, uint16_t status_in);

//! Copy constructor
   BackgroundServerCartesian(const BackgroundServerCartesian& other);

//! Destructor
   ~BackgroundServerCartesian() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundServerCartesian);
};

};

// Something like this is needed for templated classes
#include "background_server_cartesian.cc"

#endif
