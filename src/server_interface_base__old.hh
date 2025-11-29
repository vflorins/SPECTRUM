/*!
\file interface_base.hh
\brief base class for server-worker interface
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_INTERFACE_BASE_HH
#define SPECTRUM_INTERFACE_BASE_HH

#include "common/vectors.hh"
#include "common/mpi_config.hh"

namespace Spectrum {




////----------------------------------------------------------------------------------------------------------------------------------------------------
//// InterfaceBase class declaration
////----------------------------------------------------------------------------------------------------------------------------------------------------
//

//template <typename HConfig_>
//class InterfaceBase {
//public:
//
//   using HConfig = HConfig_;
//
//protected:
//
////! Number of variables
////   int n_variables = 0;
//
////! Smallest position in the domain
////   GeoVector domain_min;
//
////! Largest position in the domain
////   GeoVector domain_max;
//
//////! Current inquiry
////   Inquiry _inquiry;
////
////#ifdef USE_MPI
////

////
////#endif
//
////! Default constructor
//   InterfaceBase(void) = default;
//
//public:
//
////! Destructor
//   virtual ~InterfaceBase(void) = default;
//
////! Common set up prior to main loop
//   virtual void ServerInterfaceStart(void);
//
////! Common clean up tasks after the main loop
//   virtual void ServerInterfaceFinish(void);
//
////! Common tasks during the main loop
//   virtual int ServerFunctions(void);
//
////! Return the vector to one of the corners of the block
//   GeoVector GetDomainMin(void) const;
//
////! Return the vector to the corner opposite to that returned in "GetDomainMin()"
//   GeoVector GetDomainMax(void) const;
//
//};
//
///*!
//\author Vladimir Florinski
//\date 06/19/2020
//\return Coordinates closest to the origin
//*/
//template <typename HConfig>
//inline GeoVector InterfaceBase<HConfig>::GetDomainMin(void) const {
//   return domain_min;
//};
//
///*!
//\author Vladimir Florinski
//\date 06/19/2020
//\return Coordinates farthest from the origin
//*/
//template <typename HConfig>
//inline GeoVector InterfaceBase<HConfig>::GetDomainMax(void) const {
//   return domain_max;
//};

};


#endif
