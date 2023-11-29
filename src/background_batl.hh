/*!
\file background_batl.hh
\brief Declares a background class using data from BATL adaptive mesh
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef _BACKGROUND_BATL_HH
#define _BACKGROUND_BATL_HH

#include "background_cartesian.hh"

namespace Spectrum {

//! Readable name of the class
const std::string bg_name_batl = "BackgroundBATL";

/*!
//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundBATL class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma background interface to the AMR server
\author Vladimir Florinski
\author Juan G Alonso Guzman

Parameters: (BackgroundBase)
*/
class BackgroundBATL : public BackgroundCartesian {

public:

//! Default constructor
   BackgroundBATL(void);

//! Copy constructor
   BackgroundBATL(const BackgroundBATL& other);

//! Destructor
   ~BackgroundBATL() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundBATL);
};

};

#endif
