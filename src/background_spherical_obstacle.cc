/*!
\file background_spherical_obstacle.cc
\brief Implements a magnetic field background around a spherical object without a flow
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_spherical_obstacle.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSphericalObstacle methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 03/25/2022
*/
BackgroundSphericalObstacle::BackgroundSphericalObstacle(void)
                           : BackgroundBase(bg_name_spherical_obstacle, 0, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\date 08/20/2024
*/
BackgroundSphericalObstacle::BackgroundSphericalObstacle(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                           : BackgroundBase(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 03/25/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundSphericalObstacle::BackgroundSphericalObstacle(const BackgroundSphericalObstacle& other)
                           : BackgroundBase(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Juan G Alonso Guzman
\date 03/25/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BackgroundSphericalObstacle::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BackgroundBase::SetupBackground(false);
   container.Read(r_sphere);
   M = 0.5 * B0 * Cube(r_sphere);
   container.Read(dmax_fraction);
};

/*!
\author Juan G Alonso Guzman
\date 03/25/2022
*/
void BackgroundSphericalObstacle::EvaluateBackground(void)
{
   GeoVector posprime = _pos - r0;
   double posprimenorm = posprime.Norm();

   if (BITS_RAISED(_spdata._mask, BACKGROUND_U)) _spdata.Uvec = gv_zeros;
   if (BITS_RAISED(_spdata._mask, BACKGROUND_B)) {
      if (posprimenorm < r_sphere) _spdata.Bvec = gv_zeros;
      else {
         double r2 = Sqr(posprimenorm);
         double r5 = Cube(posprimenorm) * r2;
         _spdata.Bvec = B0 - (3.0 * (posprime * M) * posprime - r2 * M) / r5;
      };
   };
   if (BITS_RAISED(_spdata._mask, BACKGROUND_E)) _spdata.Evec = gv_zeros;
   if (posprimenorm < r_sphere) _spdata.region = 0.0;
   else _spdata.region = 1.0;

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\date 10/14/2022
*/
void BackgroundSphericalObstacle::EvaluateBackgroundDerivatives(void)
{
#if SPHERICAL_OBSTACLE_DERIVATIVE_METHOD == 0
   GeoVector posprime = _pos - r0;
   double posprimenorm = posprime.Norm();

   if (BITS_RAISED(_spdata._mask, BACKGROUND_gradU)) _spdata.gradUvec = gm_zeros;
   if (BITS_RAISED(_spdata._mask, BACKGROUND_gradB)) {
      if (posprimenorm < r_sphere) _spdata.gradBvec = gm_zeros;
      else {
         double r2 = Sqr(posprimenorm);
         double r5 = Cube(posprimenorm) * r2;
         double mdotr = M * posprime;
         GeoMatrix mr, rm, rr;
         mr.Dyadic(M,posprime);
         rm.Dyadic(posprime,M);
         rr.Dyadic(posprime);

         _spdata.gradBvec = -3.0 * (mr + rm + mdotr * (gm_unit - 5.0 * rr / r2)) / r5;
         _spdata.gradBmag = _spdata.gradBvec * _spdata.bhat;
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
\date 03/25/2022
*/
void BackgroundSphericalObstacle::EvaluateDmax(void)
{
   _spdata.dmax = fmin(dmax_fraction * (_pos - r0).Norm(), dmax0);
   LOWER_BITS(_status, STATE_INVALID);
};

};
