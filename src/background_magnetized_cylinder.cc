/*!
\file background_magnetized_cylinder.cc
\brief Implements a background of an infinite cylinder magnetized perpendicular to its axis
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_magnetized_cylinder.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundMagnetizedCylinder methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 02/13/2024
*/
BackgroundMagnetizedCylinder::BackgroundMagnetizedCylinder(void)
                            : BackgroundCylindricalObstacle(bg_name_magnetized_cylinder, 0, MODEL_STATIC)
{
};

/*!
\author Vladimir Florinski
\date 02/13/2024
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundMagnetizedCylinder::BackgroundMagnetizedCylinder(const BackgroundMagnetizedCylinder& other)
                            : BackgroundCylindricalObstacle(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\date 02/13/2024
*/
void BackgroundMagnetizedCylinder::EvaluateBackground(void)
{
   BackgroundCylindricalObstacle::EvaluateBackground();
   if(BITS_RAISED(_spdata._mask, BACKGROUND_B)) _spdata.Bvec = B0 - _spdata.Bvec;

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\date 03/11/2024
*/
void BackgroundMagnetizedCylinder::EvaluateBackgroundDerivatives(void)
{
#if MAGNETIZED_CYLINDER_DERIVATIVE_METHOD == 0

   BackgroundCylindricalObstacle::EvaluateBackgroundDerivatives();
   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradB)) _spdata.gradBvec *= -1.0;

#else
   NumericalDerivatives();
#endif
};

};
