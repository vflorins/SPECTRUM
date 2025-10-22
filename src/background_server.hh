/*!
\file background_server.hh
\brief Declares a background class using data from a grid on distributed memory
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SERVER_HH
#define SPECTRUM_BACKGROUND_SERVER_HH

#include "background_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundServer class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma background interface to a grid server
\author Juan G Alonso Guzman

Parameters: (BackgroundBase)
*/
template <typename HConfig_, typename ServerFront_>
class BackgroundServer : public BackgroundBase<HConfig_> {
private:

//! Readable name of the class
   static constexpr std::string_view bg_name = "BackgroundServer";

public:

   using HConfig = HConfig_;
   using ServerFront = ServerFront_;
   using BackgroundConfig = Cond<std::same_as<typename HConfig::BackgroundConfig, Default>, BackgroundDefault<BackgroundServer<HConfig, ServerFront>>, typename HConfig::BackgroundConfig>;
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

protected:

//! Pointer to a server object
   std::unique_ptr<ServerFront> server_front = nullptr;

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct);

   //! Calculate magnetic field magnitude
   void EvaluateAbsMag(void);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   void EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   void EvaluateBackgroundDerivatives(Coordinates&, Fields&);

//! Default constructor (protected, class not designed to be instantiated)
   BackgroundServer(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundServer(const std::string& name_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   BackgroundServer(const BackgroundServer& other);

public:

//! Destructor
   ~BackgroundServer() = default;

//! Signal the backend this client no longer needs its service
   void StopServerFront(void);

//! Return the vector to one of the corners of the block
   GeoVector GetDomainMin(void) const;

//! Return the vector to the corner opposite to that returned in "GetDomainMin()"
   GeoVector GetDomainMax(void) const;

};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
\return Coordinates closest to the origin
*/
template <typename HConfig, typename ServerFront>
inline GeoVector BackgroundServer<HConfig, ServerFront>::GetDomainMin(void) const
{
   if (server_front) return server_front->GetDomainMin();
   else return gv_zeros;
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
\return Coordinates farthest from the origin
*/
template <typename HConfig, typename ServerFront>
inline GeoVector BackgroundServer<HConfig, ServerFront>::GetDomainMax(void) const
{
   if (server_front) return server_front->GetDomainMax();
   else return gv_zeros;
};

};

// Something like this is needed for templated classes
#include "background_server.cc"

#endif
