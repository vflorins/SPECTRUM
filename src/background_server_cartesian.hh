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
\author Lucius Schoenbaum

Parameters: (BackgroundServer)
*/
template <typename HConfig_, typename ServerFront_>
class BackgroundServerCartesian : public BackgroundServer<HConfig_, ServerFront_> {
private:

//! Readable name of the class
   static constexpr std::string_view bg_name = "BackgroundServerCartesian";

public:

   using HConfig = HConfig_;
   using ServerFront = ServerFront_;
   using BackgroundConfig = Cond<std::same_as<typename HConfig::BackgroundConfig, Default>, BackgroundDefault<BackgroundServerCartesian<HConfig, ServerFront>>, typename HConfig::BackgroundConfig>;
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
   using BackgroundServer = BackgroundServer<HConfig, ServerFront>;


public:

//! Default constructor
   BackgroundServerCartesian(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundServerCartesian(const std::string& name_in, uint16_t status_in);

//! Copy constructor
   BackgroundServerCartesian(const BackgroundServerCartesian& other);

//! Destructor
   ~BackgroundServerCartesian() = default;

//! Clone function
   CloneFunctionBackground(BackgroundServerCartesian);
};

};

// Something like this is needed for templated classes
#include "background_server_cartesian.cc"

#endif
