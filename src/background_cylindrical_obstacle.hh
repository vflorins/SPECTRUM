/*!
\file background_cylindrical_obstacle.hh
\brief Declares a magnetic field background around a cylindrical object without a flow
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_CYLINDRICAL_OBSTACLE_HH
#define SPECTRUM_BACKGROUND_CYLINDRICAL_OBSTACLE_HH

#include "background_base.hh"

namespace Spectrum {

//! Method for computing derivatives of B (0: analytical, 1: numerical)
#define CYLINDRICAL_OBSTACLE_DERIVATIVE_METHOD 0

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundCylindricalObstacle class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the class
const std::string bg_name_cylindrical_obstacle = "BackgroundCylindricalObstacle";

/*!
\brief Magnetic field around a cylindrical obstacle
\author Juan G Alonso Guzman

Parameters: (BackgroundBase), GeoVector axis, double r_obstacle, double dmax_fraction
*/
template <typename Fields_>
class BackgroundCylindricalObstacle : public BackgroundBase<Fields_>
{
public:

   using Fields = Fields_;
   using BackgroundBase = BackgroundBase<Fields>;
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
   using BackgroundBase::StopServerFront;
   using BackgroundBase::SetupBackground;
   using BackgroundBase::EvaluateBackground;
   using BackgroundBase::EvaluateBackgroundDerivatives;
   using BackgroundBase::NumericalDerivatives;

protected:

//! Axis of the cylinder (persistent)
   GeoVector axis;

//! Radius of cylindrical obstacle (persistent)
   double r_cylinder;

//! Maximum fraction of the radial distance per step (persistent)
   double dmax_fraction;

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct) override;

//! Compute the internal u, B, and E fields
   void EvaluateBackground(void) override;

//! Compute the internal u, B, and E derivatives
   void EvaluateBackgroundDerivatives(void) override;

//! Compute the maximum distance per time step
   void EvaluateDmax(void) override;

public:

//! Default constructor
   BackgroundCylindricalObstacle(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundCylindricalObstacle(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor
   BackgroundCylindricalObstacle(const BackgroundCylindricalObstacle& other);

//! Destructor
   ~BackgroundCylindricalObstacle() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundCylindricalObstacle);

};

};

// Something like this is needed for templated classes
#include "background_cylindrical_obstacle.cc"

#endif
