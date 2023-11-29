/*!
\file background_cartesian.cc
\brief Implements a background class using data from a uniform Cartesian grid
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_cartesian.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundCartesian methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
*/
BackgroundCartesian::BackgroundCartesian(void)
                   : BackgroundBase(bg_name_cartesian, 0, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/27/2023
*/
BackgroundCartesian::BackgroundCartesian(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                   : BackgroundBase(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundCartesian::BackgroundCartesian(const BackgroundCartesian& other)
                   : BackgroundBase(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BackgroundCartesian::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BackgroundBase::SetupBackground(false);

   std::shared_ptr<MPI_Config>* mpi_config_ptr;
   container.Read(&mpi_config_ptr);
   
#ifdef NEED_SERVER
   if((*mpi_config_ptr)->is_worker) {
      server_front = std::make_unique<ServerFrontType>();
      server_front->ConnectMPIConfig(*mpi_config_ptr);
      server_front->ServerStart();
   };
#endif
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void BackgroundCartesian::EvaluateBackground(void)
{
#ifdef NEED_SERVER
   server_front->GetVariables(_t, _pos, _spdata);
#endif
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void BackgroundCartesian::EvaluateBackgroundDerivatives(void)
{
#ifdef NEED_SERVER
   server_front->GetGradients(_spdata);
#endif
   if(BITS_RAISED(_spdata._mask, BACKGROUND_grad_FAIL)) NumericalDerivatives();
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void BackgroundCartesian::EvaluateDmax(void)
{
// "_dmax" is actually evaluated in "EvaluateBackground", which is always called after "EvaluateDmax" by the trajectory. It is an empty override here for efficiency, since the base version executes "_dmax = dmax0" every time.
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
*/
void BackgroundCartesian::StopServerFront(void)
{
   server_front->ServerFinish();
};

};
