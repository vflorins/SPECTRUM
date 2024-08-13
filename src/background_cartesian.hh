/*!
\file background_cartesian.hh
\brief Declares a background class using data from uniform Cartesian grid on distributed memory
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef _BACKGROUND_CARTESIAN_HH
#define _BACKGROUND_CARTESIAN_HH

#include "background_server.hh"

namespace Spectrum {

//! Readable name of the class
const std::string bg_name_cartesian = "BackgroundCartesian";

/*!
//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundCartesian class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma background interface to a uniform Cartesian grid server
\author Juan G Alonso Guzman

Parameters: (BackgroundServer)
*/
class BackgroundCartesian : public BackgroundServer {

public:

//! Default constructor
   BackgroundCartesian(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundCartesian(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor
   BackgroundCartesian(const BackgroundCartesian& other);

//! Destructor
   ~BackgroundCartesian() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundCartesian);
};

};

#endif
