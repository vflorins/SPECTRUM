/*!
\file background_spherical_obstacle.hh
\brief Declares a magnetic field background around a spherical object without a flow
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SPHERICAL_OBSTACLE_HH
#define SPECTRUM_BACKGROUND_SPHERICAL_OBSTACLE_HH

#include "background_base.hh"

namespace Spectrum {

//! Method for computing derivatives of B (0: analytical, 1: numerical)
#define SPHERICAL_OBSTACLE_DERIVATIVE_METHOD 0

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSphericalObstacle class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the class
const std::string bg_name_spherical_obstacle = "BackgroundSphericalObstacle";

/*!
\brief Magnetic field around a spherical obstacle
\author Juan G Alonso Guzman

Parameters: (BackgroundBase), double r_sphere, double dmax_fraction
*/
class BackgroundSphericalObstacle : public BackgroundBase {

protected:

//! Radius of spherical obstacle (persistent)
   double r_sphere;

//! Dipole moment (persistent)
   GeoVector M;

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
   BackgroundSphericalObstacle(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundSphericalObstacle(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor
   BackgroundSphericalObstacle(const BackgroundSphericalObstacle& other);

//! Destructor
   ~BackgroundSphericalObstacle() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundSphericalObstacle);
};

};

#endif
