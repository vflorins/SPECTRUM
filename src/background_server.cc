/*!
\file background_server.cc
\brief Implements a background class using data from a grid on distributed memory
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_server.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundServer methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
*/
BackgroundServer::BackgroundServer(void)
                : BackgroundBase("", 0, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
*/
BackgroundServer::BackgroundServer(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                : BackgroundBase(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundServer::BackgroundServer(const BackgroundServer& other)
                : BackgroundBase(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BackgroundServer::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BackgroundBase::SetupBackground(false);
   
#ifdef NEED_SERVER
   if(MPI_Config::is_worker) {
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
void BackgroundServer::EvaluateBackground(void)
{
#ifdef NEED_SERVER
   server_front->GetVariables(_t, _pos, _spdata);
#endif
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
*/
void BackgroundServer::EvaluateBackgroundDerivatives(void)
{
#ifdef NEED_SERVER
   server_front->GetGradients(_spdata);
#endif
   if(BITS_RAISED(_spdata._mask, BACKGROUND_grad_FAIL)) NumericalDerivatives();
};

#if SERVER_INTERP_ORDER > 0
/*!
\author Juan G Alonso Guzman
\date 06/17/2024
*/
void BackgroundServer::EvaluateBmag(void)
{
// If variables are interpolated, override this function to be empty, and interpolation will happen within "GetVariables".
};
#endif

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
*/
void BackgroundServer::StopServerFront(void)
{
#ifdef GEO_DEBUG
   server_front->PrintStencilOutcomes();
#endif
   server_front->ServerFinish();
};

};
