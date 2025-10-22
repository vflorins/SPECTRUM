/*!
\file background_cylindrical_obstacle.cc
\brief Implements a magnetic field background around a cylindrical object without a flow
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_cylindrical_obstacle.hh"
#include <common/matrix.hh>

namespace Spectrum {

using namespace BackgroundOptions;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundCylindricalObstacle methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/06/2025
*/
template <typename HConfig>
BackgroundCylindricalObstacle<HConfig>::BackgroundCylindricalObstacle(void)
                             : BackgroundBase(bg_name, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/06/2025
*/
template <typename HConfig>
BackgroundCylindricalObstacle<HConfig>::BackgroundCylindricalObstacle(const std::string& name_in, uint16_t status_in)
                             : BackgroundBase(name_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/06/2025
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename HConfig>
BackgroundCylindricalObstacle<HConfig>::BackgroundCylindricalObstacle(const BackgroundCylindricalObstacle& other)
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
template <typename HConfig>
void BackgroundCylindricalObstacle<HConfig>::SetupBackground(bool construct)
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
\author Lucius Schoenbaum
\date 08/06/2025
Compute the internal u, B, and E fields
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
void BackgroundCylindricalObstacle<HConfig>::EvaluateBackground(Coordinates& coords, Fields& fields)
{

   GeoVector posprime = coords.Pos() - r0;
   posprime.SubtractParallel(axis);
   double posprimenorm = posprime.Norm();

   if constexpr (RequestedFields::Fluv_found())
      fields.Fluv() = gv_zeros;
   if constexpr (RequestedFields::Mag_found()) {
      if (posprimenorm < r_cylinder)
         fields.Mag() = gv_zeros;
      else {
         double s2 = posprime.Norm2();
         fields.Mag() = B0 - Sqr(r_cylinder) / s2 * (2.0 * (posprime * B0) / s2 * posprime - B0);
      };
   };
   if constexpr (RequestedFields::Elc_found())
      fields.Elc() = gv_zeros;
   if constexpr (RequestedFields::Iv1_found()) {
      if (posprimenorm < r_cylinder) fields.Iv1() = 0.0;
      else fields.Iv1() = 1.0;
   }

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 09/08/2025
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
void BackgroundCylindricalObstacle<HConfig>::EvaluateBackgroundDerivatives(Coordinates& coords, Fields& fields)
{
   if constexpr (derivative_method == DerivativeMethod::analytic) {
      GeoVector posprime = coords.Pos() - r0;
      posprime.SubtractParallel(axis);
      double posprimenorm = posprime.Norm();

      if constexpr (RequestedFields::DelFluv_found())
         fields.DelFluv() = gm_zeros;
      if constexpr (RequestedFields::DelMag_found()) {
         if (posprimenorm < r_cylinder)
            fields.DelMag() = gm_zeros;
         else {
// TODO: complete
         };
      };
      if constexpr (RequestedFields::DelElc_found()) fields.DelElc() = gm_zeros;
      if constexpr (RequestedFields::DotFluv_found()) fields.DotFluv() = gv_zeros;
      if constexpr (RequestedFields::DotMag_found()) fields.DotMag() = gv_zeros;
      if constexpr (RequestedFields::DotElc_found()) fields.DotElc() = gv_zeros;
   }
   else {
      BackgroundBase::template NumericalDerivatives<Coordinates, Fields, RequestedFields>(coords, fields);
   }
};

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/06/2025
*/
template <typename HConfig>
template <typename Coordinates>
void BackgroundCylindricalObstacle<HConfig>::EvaluateDmax(Coordinates& coords)
{
   GeoVector posprime = coords.Pos() - r0;
   posprime.SubtractParallel(axis);
   _ddata.dmax = fmin(dmax_fraction * posprime.Norm(), dmax0);
   LOWER_BITS(_status, STATE_INVALID);
};


};

