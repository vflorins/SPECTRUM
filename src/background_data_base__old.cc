/*!
\file background_data_base.cc
\brief Implements a background class using data from a grid on distributed memory
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_data_base__old.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundServer methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

///*!
//\author Juan G Alonso Guzman
//\date 07/19/2023
//*/
//template <typename HConfig, typename ServerFront>
//BackgroundServer<HConfig, ServerFront>::BackgroundServer(void)
//                : BackgroundBase("", STATE_NONE)
//{
//};
//
///*!
//\author Juan G Alonso Guzman
//\date 07/27/2023
//*/
//template <typename HConfig, typename ServerFront>
//BackgroundServer<HConfig, ServerFront>::BackgroundServer(const std::string_view& name_in, status_t status_in)
//                : BackgroundBase(name_in, status_in)
//{
//};

///*!
//\author Juan G Alonso Guzman
//\date 07/19/2023
//\param[in] other Object to initialize from
//
//A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
//*/
//template <typename HConfig, typename ServerFront>
//BackgroundServer<HConfig, ServerFront>::BackgroundServer(const BackgroundServer& other)
//                : BackgroundBase(other)
//{
//   RAISE_BITS(_status, MODEL_STATIC);
//   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
//};

///*!
//\author Vladimir Florinski
//\author Juan G Alonso Guzman
//\date 01/04/2024
//\param [in] construct Whether called from a copy constructor or separately
//
//This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
//*/
//template <typename HConfig, typename ServerFront>
//void BackgroundServer<HConfig, ServerFront>::SetupBackground(bool construct)
//{
//
//   // todo - this macro-block
//   //  can be executed constexpr-conditionally during init by trajectory,
//   //  that way backgrounds are lightweight and
//   //  much easier to reason about.
// todo UPDATE: this is replaced by the protocol that trajectory (simulation???) should call background.Start() and background.Finish() generically (for any and all backgrounds).
//    ...here, we are doing this in SetupBackground ----> TASK: trace back and see who calls SetupBackground
//#ifdef NEED_SERVER
//
//   if (MPI_Config::is_worker) {
//      server_front = std::make_unique<ServerFrontType>();
//      server_front->ServerStart();
//   };
//#endif
//};

///*!
//\author Vladimir Florinski
//\author Juan G Alonso Guzman
//\date 01/04/2024
//*/
//template <typename HConfig>
//template <typename Coordinates, typename Fields, typename RequestedFields>
//void BackgroundData<HConfig>::EvaluateBackground(Coordinates& coords, Fields& fields)
//{
//   server_front->GetVariables(coords.Time(), coords.Pos(), fields);
//};
//
///*!
//\author Vladimir Florinski
//\author Juan G Alonso Guzman
//\date 01/04/2024
//*/
//template <typename HConfig>
//template <typename Coordinates, typename Fields, typename RequestedFields>
//void BackgroundData<HConfig>::EvaluateBackgroundDerivatives(Coordinates& coords, Fields& fields)
//{
//   server_front->GetGradients(fields, _ddata);
//};





///*!
//\author Vladimir Florinski
//\author Juan G Alonso Guzman
//\author Lucius Schoenbaum
//\date 11/22/2025
//*/
//template <typename HConfig>
//void BackgroundDataBase<HConfig>::Stop(void)
//{
//#ifdef GEO_DEBUG
//   PrintStencilOutcomes();
//#endif
//   ServerFinish();
//};
//
//
///*!
//\author Vladimir Florinski
//\date 10/28/2022
//*/
//template <typename HConfig>
//void BackgroundDataBase<HConfig>::Start(void)
//{
//   ServerInterfaceBase::ServerInterfaceStart();
//   cache_line.Empty();
//};
//
///*!
//\author Vladimir Florinski
//\date 10/28/2022
//*/
//template <typename HConfig>
//void BackgroundDataBase<HConfig>::Finish(void)
//{
//   // todo wrap in interface
//   MPI_Send(nullptr, 0, MPI_BYTE, 0, tag_stopserve, MPI::node_comm);
//   ServerInterfaceBase::ServerInterfaceFinish();
//};
//
///*!
//\author Vladimir Florinski
//\date 10/28/2022
//*/
//template <typename HConfig>
//int BackgroundDataBase<HConfig>::ServerFunctions(void)
//{
//   return 0;
//};
//
///*!
//\author Vladimir Florinski
//\date 06/19/2020
//\return Number of blocks in the cache
//*/
//template <typename HConfig>
//int BackgroundDataBase<HConfig>::GetNCachedBlocks(void)
//{
//   return cache_line.size();
//};
//
///*!
//\author Vladimir Florinski
//\date 01/27/2023
//*/
//template <typename HConfig>
//void BackgroundDataBase<HConfig>::InvalidateCache(void)
//{
//   cache_line.Empty();
//};



};
