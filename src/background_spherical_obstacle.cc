/*!
\file background_spherical_obstacle.cc
\brief Implements a magnetic field background around a spherical object without a flow
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_spherical_obstacle.hh"
#include "common/matrix.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSphericalObstacle methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 03/25/2022
*/
template <typename HConfig>
BackgroundSphericalObstacle<HConfig>::BackgroundSphericalObstacle(void)
                           : BackgroundBase(bg_name, MODEL_STATIC)
{
};

/*!
\author Juan G Alonso Guzman
\date 08/20/2024
*/
template <typename HConfig>
BackgroundSphericalObstacle<HConfig>::BackgroundSphericalObstacle(const std::string& name_in, uint16_t status_in)
                           : BackgroundBase(name_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 03/25/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename HConfig>
BackgroundSphericalObstacle<HConfig>::BackgroundSphericalObstacle(const BackgroundSphericalObstacle& other)
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
template <typename HConfig>
void BackgroundSphericalObstacle<HConfig>::SetupBackground(bool construct)
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
template <typename HConfig>
template <typename Fields>
void BackgroundSphericalObstacle<HConfig>::EvaluateBackground(Coordinates& coords, Fields& fields)
{
   GeoVector posprime = coords.Pos() - r0;
   double posprimenorm = posprime.Norm();

   if constexpr (Fields::Vel_found()) fields.Vel() = gv_zeros;
   if constexpr (Fields::Mag_found()) {
      if (posprimenorm < r_sphere) fields.Mag() = gv_zeros;
      else {
         double r2 = Sqr(posprimenorm);
         double r5 = Cube(posprimenorm) * r2;
         fields.Mag() = B0 - (3.0 * (posprime * M) * posprime - r2 * M) / r5;
      };
   };
   if constexpr (Fields::Elc_found()) fields.Elc() = gv_zeros;
   if constexpr (Fields::Iv0_found()) {
      if (posprimenorm < r_sphere) fields.Iv0() = 0.0;
      else fields.Iv0() = 1.0;
   }

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\date 10/14/2022
*/
template <typename HConfig>
template <typename Fields>
void BackgroundSphericalObstacle<HConfig>::EvaluateBackgroundDerivatives(Coordinates& coords, Specie& specie, Fields& fields)
{
   if constexpr (HConfig::derivative_method == DerivativeMethod::analytic) {
      GeoVector posprime = coords.Pos() - r0;
      double posprimenorm = posprime.Norm();

      if constexpr (Fields::DelVel_found()) fields.DelVel() = gm_zeros;
      if constexpr (Fields::DelMag_found() || Fields::DelAbsMag_found()) {
         if (posprimenorm < r_sphere) {
            if constexpr (Fields::DelMag_found()) fields.DelMag() = gm_zeros;
            if constexpr (Fields::DelAbsMag_found()) fields.DelAbsMag() = 0.0;
         }
         else {
            double r2 = Sqr(posprimenorm);
            double r5 = Cube(posprimenorm) * r2;
            double mdotr = M * posprime;
            GeoMatrix mr, rm, rr;
            mr.Dyadic(M,posprime);
            rm.Dyadic(posprime,M);
            rr.Dyadic(posprime);

            auto DelMag = -3.0 * (mr + rm + mdotr * (gm_unit - 5.0 * rr / r2)) / r5;
            if constexpr (Fields::DelMag_found())
               fields.DelMag() = DelMag;
            if constexpr (Fields::DelAbsMag_found()) {
               GeoVector bhat;
               if constexpr (Fields::HatMag_found) bhat = fields.HatMag();
               else bhat = fields.Mag() / fields.Mag().Norm();
               fields.DelAbsMag() = DelMag * bhat;
            }
         };
      };
      if constexpr (Fields::DelElc_found()) fields.DelElc() = gm_zeros;
      if constexpr (Fields::DotVel_found()) fields.DotVel() = gv_zeros;
      if constexpr (Fields::DotMag_found()) fields.DotMag() = gv_zeros;
      if constexpr (Fields::DotElc_found()) fields.DotElc() = gv_zeros;
   }
   else {
      NumericalDerivatives();
   };
};

/*!
\author Juan G Alonso Guzman
\date 03/25/2022
*/
template <typename HConfig>
void BackgroundSphericalObstacle<HConfig>::EvaluateDmax(Coordinates& coords)
{
   _ddata.dmax = fmin(dmax_fraction * (coords.Pos() - r0).Norm(), dmax0);
   LOWER_BITS(_status, STATE_INVALID);
};

};
