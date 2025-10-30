/*!
\file background_shock.hh
\brief Declares a simple planar MHD shock field background
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a SHOCK coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SHOCK_HH
#define SPECTRUM_BACKGROUND_SHOCK_HH

#include "background_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundShock class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BackgroundShock class
const std::string bg_name_shock = "BackgroundShock";

/*!
\brief Planar MHD shock
\author Juan G Alonso Guzman

Parameters: (BackgroundBase), GeoVector n_shock, double v_shock, double compression
*/
class BackgroundShock : public BackgroundBase {

protected:

//! Shock normal (persistent)
   GeoVector n_shock;

//! Shock velocity (persistent)
   double v_shock;

//! Compression ratio (persistent)
   double compression;

//! Downstream velocity (persistent), "u0" is upstream flow vector
   GeoVector u1;

//! Downstream magnetic field (persistent), "B0" is upstream magnetic field
   GeoVector B1;

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct) override;

//! Compute the internal u, B, and E fields
   void EvaluateBackground(void) override;

//! Compute the internal u, B, and E derivatives
   void EvaluateBackgroundDerivatives(void) override;

public:

//! Default constructor
   BackgroundShock(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundShock(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor
   BackgroundShock(const BackgroundShock& other);

//! Destructor
   ~BackgroundShock() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundShock);
};

};

#endif
