/*!
\file background_magcylinder.cc
\brief Implements a background of an infinite cylinder magnetized perpendicular to its axis
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_magnetized_cylinder.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundDipole methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 02/13/2024
*/
BackgroundMagnetizedCylinder::BackgroundMagnetizedCylinder(void)
                            : BackgroundBase(bg_name_magnetized_cylinder, 0, MODEL_STATIC)
{
};

/*!
\author Vladimir Florinski
\date 02/13/2024
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundMagnetizedCylinder::BackgroundMagnetizedCylinder(const BackgroundMagnetizedCylinder& other)
                            : BackgroundBase(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\date 02/13/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BackgroundMagnetizedCylinder::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BackgroundBase::SetupBackground(false);

   container.Read(axis.Data());
   axis.Normalize();

// "B0" must be normal to the axis
   B0.SubtractParallel(axis);

   container.Read(&r_cylinder);
   container.Read(&dmax_fraction);
};

/*!
\author Vladimir Florinski
\date 02/13/2024
*/
void BackgroundMagnetizedCylinder::EvaluateBackground(void)
{
   GeoVector posprime = _pos - r0;
   posprime.SubtractParallel(axis);
   double posprimenorm = posprime.Norm();

   if(BITS_RAISED(_spdata._mask, BACKGROUND_U)) _spdata.Uvec = gv_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_B)) {
      if(posprimenorm < r_cylinder) _spdata.Bvec = B0;
      else {
         double s2 = posprime.Norm2();
         _spdata.Bvec = Sqr(r_cylinder) / s2 * (2.0 / s2 * (posprime * B0) * posprime - B0);
      };
   };
   if(BITS_RAISED(_spdata._mask, BACKGROUND_E)) _spdata.Evec = gv_zeros;
   if(posprimenorm < r_cylinder) _spdata.region = 0.0;
   else _spdata.region = 1.0;

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\date 02/13/2024
*/
void BackgroundMagnetizedCylinder::EvaluateDmax(void)
{
   GeoVector posprime = _pos - r0;
   posprime.SubtractParallel(axis);
   _spdata.dmax = fmin(dmax_fraction * posprime.Norm(), dmax0);
   LOWER_BITS(_status, STATE_INVALID);
};

};
