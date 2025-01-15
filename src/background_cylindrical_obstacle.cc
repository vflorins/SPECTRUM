/*!
\file background_cylindrical_obstacle.cc
\brief Implements a magnetic field background around a cylindrical object without a flow
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_cylindrical_obstacle.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundCylindricalObstacle methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 08/20/2024
*/
BackgroundCylindricalObstacle::BackgroundCylindricalObstacle(void)
                             : BackgroundBase(bg_name_cylindrical_obstacle, 0, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\date 08/20/2024
*/
BackgroundCylindricalObstacle::BackgroundCylindricalObstacle(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                             : BackgroundBase(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 08/20/2024
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundCylindricalObstacle::BackgroundCylindricalObstacle(const BackgroundCylindricalObstacle& other)
                             : BackgroundBase(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Juan G Alonso Guzman
\date 08/20/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BackgroundCylindricalObstacle::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BackgroundBase::SetupBackground(false);

   container.Read(axis);
   axis.Normalize();

// "B0" must be normal to the axis
   B0.SubtractParallel(axis);

   container.Read(r_cylinder);
   container.Read(dmax_fraction);
};

/*!
\author Juan G Alonso Guzman
\date 08/20/2024
*/
void BackgroundCylindricalObstacle::EvaluateBackground(void)
{
   GeoVector posprime = _pos - r0;
   posprime.SubtractParallel(axis);
   double posprimenorm = posprime.Norm();

   if (BITS_RAISED(_spdata._mask, BACKGROUND_U)) _spdata.Uvec = gv_zeros;
   if (BITS_RAISED(_spdata._mask, BACKGROUND_B)) {
      if (posprimenorm < r_cylinder) _spdata.Bvec = gv_zeros;
      else {
         double s2 = posprime.Norm2();
         _spdata.Bvec = B0 - Sqr(r_cylinder) / s2 * (2.0 * (posprime * B0) / s2 * posprime - B0);
      };
   };
   if (BITS_RAISED(_spdata._mask, BACKGROUND_E)) _spdata.Evec = gv_zeros;
   if (posprimenorm < r_cylinder) _spdata.region = 0.0;
   else _spdata.region = 1.0;

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\date 08/20/2024
*/
void BackgroundCylindricalObstacle::EvaluateBackgroundDerivatives(void)
{
#if SPHERICAL_OBSTACLE_DERIVATIVE_METHOD == 0
   int xyz;
   GeoVector posprime = _pos - r0;
   posprime.SubtractParallel(axis);
   double posprimenorm = posprime.Norm();

   if (BITS_RAISED(_spdata._mask, BACKGROUND_gradU)) _spdata.gradUvec = gm_zeros;
   if (BITS_RAISED(_spdata._mask, BACKGROUND_gradB)) {
      if (posprimenorm < r_cylinder) _spdata.gradBvec = gm_zeros;
      else {
// TODO: complete
      };
   };
   if (BITS_RAISED(_spdata._mask, BACKGROUND_gradE)) _spdata.gradEvec = gm_zeros;
   if (BITS_RAISED(_spdata._mask, BACKGROUND_dUdt)) _spdata.dUvecdt = gv_zeros;
   if (BITS_RAISED(_spdata._mask, BACKGROUND_dBdt)) _spdata.dBvecdt = gv_zeros;
   if (BITS_RAISED(_spdata._mask, BACKGROUND_dEdt)) _spdata.dEvecdt = gv_zeros;

#else
   NumericalDerivatives(); 
#endif
};

/*!
\author Juan G Alonso Guzman
\date 08/20/2024
*/
void BackgroundCylindricalObstacle::EvaluateDmax(void)
{
   GeoVector posprime = _pos - r0;
   posprime.SubtractParallel(axis);
   _spdata.dmax = fmin(dmax_fraction * posprime.Norm(), dmax0);
   LOWER_BITS(_status, STATE_INVALID);
};

};
