/*!
\file background_cartesian.hh
\brief Declares a background class using data from uniform Cartesian grid
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef _BACKGROUND_CARTESIAN_HH
#define _BACKGROUND_CARTESIAN_HH

#include "background_base.hh"

namespace Spectrum {

//! Readable name of the class
const std::string bg_name_cartesian = "BackgroundCartesian";

/*!
//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundCartesian class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma background interface to the AMR server
\author Juan G Alonso Guzman

Parameters: (BackgroundBase)
*/
class BackgroundCartesian : public BackgroundBase {

protected:

//! Pointer to a server object
   std::unique_ptr<ServerFrontType> server_front = nullptr;

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct) override;

//! Compute the internal u, B, and E fields
   void EvaluateBackground(void) override;

//! Compute the internal derivatives of the fields
   void EvaluateBackgroundDerivatives(void) override;

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

//! Signal the backend this client no longer needs its service
   void StopServerFront(void) override;

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
inline GeoVector BackgroundCartesian::GetDomainMin(void) const
{
   if(server_front) return server_front->GetDomainMin();
   else return gv_zeros;
};

/*!
\author Juan G Alonso Guzman
\date 07/19/2023
\return Coordinates farthest from the origin
*/
inline GeoVector BackgroundCartesian::GetDomainMax(void) const
{
   if(server_front) return server_front->GetDomainMax();
   else return gv_zeros;
};

};

#endif
