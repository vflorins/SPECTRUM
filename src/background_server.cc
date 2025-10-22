/*!
\file background_server.cc
\brief Implements a background class using data from a grid on distributed memory
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_server.hh"

namespace Spectrum {

using namespace BackgroundOptions;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundServer methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
*/
template <typename HConfig, typename ServerFront>
BackgroundServer<HConfig, ServerFront>::BackgroundServer(void)
                : BackgroundBase("", STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
*/
template <typename HConfig, typename ServerFront>
BackgroundServer<HConfig, ServerFront>::BackgroundServer(const std::string& name_in, uint16_t status_in)
                : BackgroundBase(name_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename HConfig, typename ServerFront>
BackgroundServer<HConfig, ServerFront>::BackgroundServer(const BackgroundServer& other)
                : BackgroundBase(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig, typename ServerFront>
void BackgroundServer<HConfig, ServerFront>::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BackgroundBase::SetupBackground(false);
   
#ifdef NEED_SERVER

   if (MPI_Config::is_worker) {
      server_front = std::make_unique<ServerFrontType>();
      server_front->ServerStart();
   };
#endif
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
*/
template <typename HConfig, typename ServerFront>
template <typename Coordinates, typename Fields, typename RequestedFields>
void BackgroundServer<HConfig, ServerFront>::EvaluateBackground(Coordinates& coords, Fields& fields)
{
#ifdef NEED_SERVER
   server_front->GetVariables(coords.Time(), coords.Pos(), fields);
#endif
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
*/
template <typename HConfig, typename ServerFront>
template <typename Coordinates, typename Fields, typename RequestedFields>
void BackgroundServer<HConfig, ServerFront>::EvaluateBackgroundDerivatives(Coordinates& coords, Fields& fields)
{
#ifdef NEED_SERVER
   server_front->GetGradients(fields, _ddata);
#endif
   if (_ddata.BACKGROUND_grad_FAIL)
      BackgroundBase::template NumericalDerivatives<Coordinates, Fields, RequestedFields>(coords, fields);
};

/*!
\author Juan G Alonso Guzman
\date 06/17/2024
*/
template <typename HConfig, typename ServerFront>
void BackgroundServer<HConfig, ServerFront>::EvaluateAbsMag(void)
{
   if constexpr (BackgroundConfig::server_interpolation_order > 0) {
      BackgroundBase::EvaluateAbsMag();
   }
   else {
// If variables are interpolated, override this function to be empty, and interpolation will occur within "GetVariables".
      ;
   }
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
*/
template <typename HConfig, typename ServerFront>
void BackgroundServer<HConfig, ServerFront>::StopServerFront(void)
{
#ifdef GEO_DEBUG
   server_front->PrintStencilOutcomes();
#endif
   server_front->ServerFinish();
};

};
