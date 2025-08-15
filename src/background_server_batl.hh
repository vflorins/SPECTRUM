/*!
\file background_server_batl.hh
\brief Declares a background class using data from BATL adaptive mesh on distributed memory
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SERVER_BATL_HH
#define SPECTRUM_BACKGROUND_SERVER_BATL_HH

#include "background_server_cartesian.hh"

namespace Spectrum {

//! Readable name of the class
const std::string bg_name_server_batl = "BackgroundServerBATL";


//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundServerBATL class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma background interface to a BATL server
\author Vladimir Florinski
\author Juan G Alonso Guzman

Parameters: (BackgroundServerCartesian)
*/
template <typename Fields_>
class BackgroundServerBATL : public BackgroundServerCartesian<Fields_> {
public:

   using Fields = Fields_;
   using BackgroundBase = BackgroundBase<Fields>;
   using BackgroundServerCartesian = BackgroundServerCartesian<Fields>;
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
   using BackgroundBase::EvaluateBackground;
   using BackgroundBase::EvaluateBackgroundDerivatives;
   using BackgroundBase::NumericalDerivatives;

public:

//! Default constructor
   BackgroundServerBATL(void);

//! Copy constructor
   BackgroundServerBATL(const BackgroundServerBATL& other);

//! Destructor
   ~BackgroundServerBATL() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundServerBATL);
};

};

// Something like this is needed for templated classes
#include "background_server_batl.cc"

#endif
