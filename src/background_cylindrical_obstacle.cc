/*!
\file background_cylindrical_obstacle.cc
\brief Implements a magnetic field background around a cylindrical object without a flow
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_cylindrical_obstacle.hh"
#include <common/matrix.hh>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundCylindricalObstacle methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/06/2025
*/
template <typename Fields>
BackgroundCylindricalObstacle<Fields>::BackgroundCylindricalObstacle(void)
                             : BackgroundBase(bg_name_cylindrical_obstacle, 0, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/06/2025
*/
template <typename Fields>
BackgroundCylindricalObstacle<Fields>::BackgroundCylindricalObstacle(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                             : BackgroundBase(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/06/2025
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename Fields>
BackgroundCylindricalObstacle<Fields>::BackgroundCylindricalObstacle(const BackgroundCylindricalObstacle& other)
                             : BackgroundBase(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/06/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename Fields>
void BackgroundCylindricalObstacle<Fields>::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) SetupBackground(false);

   container.Read(axis);
   axis.Normalize();

// "B0" must be normal to the axis
   B0.SubtractParallel(axis);

   container.Read(r_cylinder);
   container.Read(dmax_fraction);
};

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/06/2025
Compute the internal u, B, and E fields
*/
template <typename Fields>
void BackgroundCylindricalObstacle<Fields>::EvaluateBackground(void)
{

   GeoVector posprime = _pos - r0;
   posprime.SubtractParallel(axis);
   double posprimenorm = posprime.Norm();

   if constexpr (Fields::Vel_found())
      _fields.Vel() = gv_zeros;
   if constexpr (Fields::Mag_found()) {
      if (posprimenorm < r_cylinder)
         _fields.Mag() = gv_zeros;
      else {
         double s2 = posprime.Norm2();
         _fields.Mag() = B0 - Sqr(r_cylinder) / s2 * (2.0 * (posprime * B0) / s2 * posprime - B0);
      };
   };
   if constexpr (Fields::Elc_found())
      _fields.Elc() = gv_zeros;
   if constexpr (Fields::Iv1_found()) {
      if (posprimenorm < r_cylinder) _fields.Iv1() = 0.0;
      else _fields.Iv1() = 1.0;
   }

   LOWER_BITS(this->_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/06/2025
*/
template <typename Fields>
void BackgroundCylindricalObstacle<Fields>::EvaluateBackgroundDerivatives(void)
{

#if SPHERICAL_OBSTACLE_DERIVATIVE_METHOD == 0
   int xyz;
   GeoVector posprime = _pos - r0;
   posprime.SubtractParallel(axis);
   double posprimenorm = posprime.Norm();

   if constexpr (Fields::DelVel_found())
      _fields.DelVel() = gm_zeros;
   if constexpr (Fields::DelMag_found()) {
      if (posprimenorm < r_cylinder) _fields.DelMag() = gm_zeros;
      else {
// TODO: complete
      };
   };
   if constexpr (Fields::DelElc_found()) _fields.DelElc() = gm_zeros;
   if constexpr (Fields::DdtVel_found()) _fields.DdtVel() = gv_zeros;
   if constexpr (Fields::DdtMag_found()) _fields.DdtMag() = gv_zeros;
   if constexpr (Fields::DdtElc_found()) _fields.DdtElc() = gv_zeros;

#else
   NumericalDerivatives(); 
#endif
};

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/06/2025
*/
template <typename Fields>
void BackgroundCylindricalObstacle<Fields>::EvaluateDmax(void)
{
   GeoVector posprime = _pos - r0;
   posprime.SubtractParallel(axis);
   _ddata.dmax = fmin(dmax_fraction * posprime.Norm(), dmax0);
   LOWER_BITS(this->_status, STATE_INVALID);
};


};

