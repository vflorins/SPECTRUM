/*!
\file background_magnetized_sphere.cc
\brief Implements a background of an infinite sphere magnetized perpendicular to its axis
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_magnetized_sphere.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundMagnetizedSphere methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 08/20/2024
*/
BackgroundMagnetizedSphere::BackgroundMagnetizedSphere(void)
                          : BackgroundSphericalObstacle(bg_name_magnetized_sphere, 0, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\date 08/20/2024
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundMagnetizedSphere::BackgroundMagnetizedSphere(const BackgroundMagnetizedSphere& other)
                          : BackgroundSphericalObstacle(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Juan G Alonso Guzman
\date 08/20/2024
*/
void BackgroundMagnetizedSphere::EvaluateBackground(void)
{
   BackgroundSphericalObstacle::EvaluateBackground();
   if(BITS_RAISED(_spdata._mask, BACKGROUND_B)) _spdata.Bvec = B0 - _spdata.Bvec;

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\date 08/20/2024
*/
void BackgroundMagnetizedSphere::EvaluateBackgroundDerivatives(void)
{
#if MAGNETIZED_SPHERE_DERIVATIVE_METHOD == 0

   BackgroundSphericalObstacle::EvaluateBackgroundDerivatives();
   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradB)) _spdata.gradBvec *= -1.0;

#else
   NumericalDerivatives();
#endif
};

};
