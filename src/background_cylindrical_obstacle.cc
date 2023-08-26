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
\date 05/19/2023
*/
BackgroundCylindricalObstacle::BackgroundCylindricalObstacle(void)
                             : BackgroundBase(bg_name_cylindrical_obstacle, 0, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundCylindricalObstacle::BackgroundCylindricalObstacle(const BackgroundCylindricalObstacle& other)
                             : BackgroundBase(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BackgroundCylindricalObstacle::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BackgroundBase::SetupBackground(false);
   container.Read(&r_obstacle);
   r_obstacle2 = Sqr(r_obstacle);
   container.Read(&dmax_fraction);
   B0mag = B0.Norm();
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
*/
void BackgroundCylindricalObstacle::EvaluateBackground(void)
{
   GeoVector posprime = _pos - r0;
   double posprimenorm = posprime.Norm();
   double denom = Sqr(Sqr(posprime[0]) + Sqr(posprime[1]));
   if(BITS_RAISED(_spdata._mask, BACKGROUND_U)) _spdata.Uvec = gv_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_B)) {
      if(posprimenorm < r_obstacle) _spdata.Bvec = gv_zeros;
      else {
         _spdata.Bvec[0] = B0mag * (1.0 + r_obstacle2 * (Sqr(posprime[1]) - Sqr(posprime[0])) / denom);
         _spdata.Bvec[1] = -B0mag * 2.0 * r_obstacle2 * posprime[1] * posprime[0] / denom;
         _spdata.Bvec[2] = 0.0;
      };
   };
   if(BITS_RAISED(_spdata._mask, BACKGROUND_E)) _spdata.Evec = gv_zeros;
   if(posprimenorm < r_obstacle) _spdata.region = 0.0;
   else _spdata.region = 1.0;

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
*/
void BackgroundCylindricalObstacle::EvaluateBackgroundDerivatives(void)
{
#if SPHERICAL_OBSTACLE_DERIVATIVE_METHOD == 0
   GeoVector posprime = _pos - r0;
   double posprimenorm = posprime.Norm();
   double x2 = Sqr(posprime[0]);
   double y2 = Sqr(posprime[1]);
   double denom = Cube(x2 + y2);
   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradU)) _spdata.gradUvec = gm_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradB)) {
      if(posprimenorm < r_obstacle) _spdata.gradBvec = gm_zeros;
      else {
         _spdata.gradBvec[0][0] = B0mag * r_obstacle2 * 2.0 * posprime[0] * (x2 - 3.0 * y2) / denom;
         _spdata.gradBvec[1][0] = B0mag * r_obstacle2 * 2.0 * posprime[1] * (3.0 * x2 - y2) / denom;
         _spdata.gradBvec[0][1] = _spdata.gradBvec[1][0];
         _spdata.gradBvec[1][1] = -_spdata.gradBvec[0][0];
      };
   };
   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradE)) _spdata.gradEvec = gm_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_dUdt)) _spdata.dUvecdt = gv_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_dBdt)) _spdata.dBvecdt = gv_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_dEdt)) _spdata.dEvecdt = gv_zeros;

#else
   NumericalDerivatives(); 
#endif
};

/*!
\author Juan G Alonso Guzman
\date 05/19/2023
*/
void BackgroundCylindricalObstacle::EvaluateDmax(void)
{
// TODO
   _spdata.dmax = fmin(dmax_fraction * (_pos - r0).Norm(), dmax0);
   LOWER_BITS(_status, STATE_INVALID);
};

};
