/*!
\file background_batl.cc
\brief Implements a background class using data from BATL adaptive mesh
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_batl.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundBATL methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 10/05/2022
*/
BackgroundBATL::BackgroundBATL(void)
              : BackgroundCartesian(bg_name_batl, 0, MODEL_STATIC)
{
};

/*!
\author Vladimir Florinski
\date 10/05/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundBATL::BackgroundBATL(const BackgroundBATL& other)
              : BackgroundCartesian(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/27/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BackgroundBATL::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BackgroundBase::SetupBackground(false);

   std::shared_ptr<MPI_Config>* mpi_config_ptr;
   container.Read(&mpi_config_ptr);

   if((*mpi_config_ptr)->is_worker) {
      server_front = std::make_unique<ServerBATLFront>();
      server_front->ConnectMPIConfig(*mpi_config_ptr);
      server_front->ServerStart();
   };
};

};
