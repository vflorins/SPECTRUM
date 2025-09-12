/*!
\file source_other.hh
\brief Declares several classes to compute source terms
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SOURCE_OTHER_HH
#define SPECTRUM_SOURCE_OTHER_HH

#include "source_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SourceConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the SourceConstant class
const std::string source_name_constant = "SourceConstant";

/*!
\brief Uniform source
\author Juan G Alonso Guzman

Parameters: (SourceBase), double S0
*/
class SourceConstant : public SourceBase {

protected:

//! Source value (persistent)
   double S0;

//! Set up the diffusion model based on "params"
   void SetupSource(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateSource(void) override;

public:

//! Default constructor
   SourceConstant(void);

//! Copy constructor
   SourceConstant(const SourceConstant& other);

//! Destructor
   ~SourceConstant() override = default;

//! Clone function
   CloneFunctionSource(SourceConstant);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SourceMomentumInjection class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the SourceMomentumInjection class
const std::string source_name_momentum_injection = "SourceMomentumInjection";

/*!
\brief Uniform injection at some threshold momentum
\author Juan G Alonso Guzman

Parameters: (SourceBase), double p_inj, double rate
*/
class SourceMomentumInjection : public SourceBase {

protected:

//! Injection momentum (persistent)
   double p_inj;

//! Injection rate as a scaling factor (persistent)
   double rate;

//! Variables to keep track of injection threshold crossing (transient)
   double del_mom, del_mom_old;

//! Set up the diffusion model based on "params"
   void SetupSource(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateSource(void) override;

public:

//! Default constructor
   SourceMomentumInjection(void);

//! Copy constructor
   SourceMomentumInjection(const SourceMomentumInjection& other);

//! Destructor
   ~SourceMomentumInjection() override = default;

//! Clone function
   CloneFunctionSource(SourceMomentumInjection);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SourceSphericalShockInjection class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the SourceSphericalShockInjection class
const std::string source_name_spherical_shock_injection = "SourceSphericalShockInjection";

/*!
\brief Uniform injection at some threshold momentum for a spherical shock
\author Juan G Alonso Guzman

Parameters: (SourceBase), double p_inj, GeoVector r0, double r_sh, double w_sh, double rate
*/
class SourceSphericalShockInjection : public SourceBase {

protected:

//! Injection momentum (persistent)
   double p_inj;

//! Center of spherical shock (persistent)
   GeoVector r0;

//! Radius of spherical shock (persistent)
   double r_sh;

//! Width of spherical shock (persistent)
   double w_sh;

//! Injection rate as a scaling factor (persistent)
   double rate;

//! Variables to keep track of injection threshold crossing (transient)
   double del_mom, del_mom_old;

//! Set up the diffusion model based on "params"
   void SetupSource(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateSource(void) override;

public:

//! Default constructor
   SourceSphericalShockInjection(void);

//! Copy constructor
   SourceSphericalShockInjection(const SourceSphericalShockInjection& other);

//! Destructor
   ~SourceSphericalShockInjection() override = default;

//! Clone function
   CloneFunctionSource(SourceSphericalShockInjection);
};

};

#endif
