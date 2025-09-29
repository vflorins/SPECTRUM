/*!
\file background_vlism_bochum.cc
\brief Implements a plasma background class for the Very Local Interstellar Medium with the Roken/Kleimann analytic field model
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_vlism_bochum.hh"
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_errno.h>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundVLISMBochum methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 09/28/2021
*/
template <typename HConfig>
BackgroundVLISMBochum<HConfig>::BackgroundVLISMBochum(void)
                     : BackgroundBase(bg_name, MODEL_STATIC)
{
// https://www.gnu.org/software/gsl/doc/html/err.html#c.gsl_set_error_handler
// Turn gsl error handler off
   gsl_error_handler_t* gsl_default_error_handler = gsl_set_error_handler_off();
};

/*!
\author Jens Kleimann
\date 10/30/2019
\param[in] z Normalized z-component of position
\return Normalized transverse field
*/
template <typename HConfig>
double BackgroundVLISMBochum<HConfig>::RelBtrans(double z) const
{
   return 1.0 / sqrt(1.0 - 1.0 / Sqr(z));
};

/*!
\author Vladimir Florinski
\date 09/28/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename HConfig>
BackgroundVLISMBochum<HConfig>::BackgroundVLISMBochum(const BackgroundVLISMBochum& other)
                     : BackgroundBase(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/26/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BackgroundVLISMBochum<HConfig>::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BackgroundBase::SetupBackground(false);
   container.Read(z_nose);

// Build the new coordinate system with z axis along -u0 and convert the magnetic field to the primed frame.
   eprime[2] = -UnitVec(u0);
   eprime[0] = GetSecondUnitVec(eprime[2]);
   eprime[1] = eprime[2] ^ eprime[0];

   u0.ChangeToBasis(eprime);
   B0.ChangeToBasis(eprime);

   sin_theta_B0 = fmax(sp_tiny, sqrt(1.0 - Sqr(UnitVec(B0) * eprime[2])));

   if constexpr (HConfig::VLISMBochum_mod_type == 3) {
      double ht = Sqr(scB / sin_theta_B0) * (Sqr(ztr) - 1.0);
      fzoom = sqrt((ht - 1.0) / (ht - Sqr(ztr)));
   }

};

/*!
\author Jens Kleimann
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/26/2022
\return Field amplification factor
*/
template <typename HConfig>
double BackgroundVLISMBochum<HConfig>::GetAmpFactor(double zeta) const
{
   if constexpr (HConfig::VLISMBochum_mod_type == 0) {
      return 1.0;
   }

// Where to evaluate the custom amplification function.
   double zev;

// If "mod_rpos" is zero, the scaling is done with respect to the intersection point of the isochrone with the z-axis (s=0). To find "zev" the nonlinear equaiton zev+(1/2)ln[(zev-1)/(zev+1)]=zeta is solved via Newton iterations.
   if constexpr (HConfig::VLISMBochum_mod_rpos == 0) {
      // Current error of iterator loop
      double f0, f1, delta = 1.0;

// Newton iterations
      zev = 2.0;
      while (fabs(delta) > sp_tiny) {
         f0 = (zev + 0.5 * log((zev - 1.0) / (zev + 1.0)) - zeta);
         f1 = (zev * zev) / (zev * zev - 1.0);
         delta = f0 / f1;

// Catch overshoot before doing log[<0] by just going halfway to 1.0
         if (zev - delta <= 1.0) delta = 0.5 * (zev - 1.0);
         zev -= delta;
      };
   }
   else {
// If "mod_rpos" is one, the scaling is done with respect to the isochrone asymptotic position at s=inf, so zev=zeta.
      zev = zeta;
   }

   if constexpr (HConfig::VLISMBochum_mod_type == 1) {
      // No modification for zev>ztr, zero for zev<ztr
      return (zev > ztr ? 1.0 : 0.0);
   }
   else if constexpr (HConfig::VLISMBochum_mod_type == 2) {
      // No modification for zev>ztr, constant for zev<ztr
      if (fabs(sin_theta_B0) < sp_tiny) return 1.0;
      else return (zev > ztr ? 1.0 : RelBtrans(ztr) / RelBtrans(zev));
   }
   else if constexpr (HConfig::VLISMBochum_mod_type == 3) {
      // No modification for zev>ztr, scaled to reach scB at the nose for zev<ztr
      if ((fabs(sin_theta_B0) < sp_tiny) || (zev > ztr)) return 1.0;
      else return RelBtrans(fzoom * zev) / RelBtrans(fzoom * ztr) / RelBtrans(zev) * RelBtrans(ztr);
   }
   else {
      // Unmodified field
      return 1.0;
   }

};

/*!
\author Jens Kleimann
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/26/2022
*/
template <typename HConfig>
template <typename Fields>
void BackgroundVLISMBochum<HConfig>::EvaluateBackground(BackgroundCoordinates& coords, Fields& fields)
{
// Convert position into flow aligned coordinates and scale to "z_nose"
   GeoVector posprime = coords.Pos() - r0;
   posprime.ChangeToBasis(eprime);
   posprime /= z_nose;

// Position in cylindrical coordinates
   double s, z, s2, r, r3, a, a2, a2s2, a1s1;
   z = posprime[2];
   s2 = Sqr(posprime[0]) + Sqr(posprime[1]);
   s = sqrt(s2);
   r = sqrt(s2 + Sqr(z));
   r3 = Cube(r);

// This is the "impact parameter". If it is less than zero then the point is inside the HP, which is an error.
   a2 = s2 - 2.0 * (1.0 - z / r);

// Test for inside of the HP
   if ((z < 1.0) && ((s < sp_tiny) || (a2 < 0.0))) {
      fields.Mag() = gv_zeros;
// No need to check computation flags since it's just assigning zeros
      fields.Elc() = gv_zeros;
      fields.Iv0() = -1.0;

      RAISE_BITS(_status, STATE_INVALID);
      return;
   }
   else fields.Iv0() = 1.0;

// Compute the velocity (valid in every region). The velocity is transformed to the global frame and un-normalized.
   fields.Vel() = (posprime / r3 - gv_nz) * u0.Norm();
   fields.Vel().ChangeFromBasis(eprime);

//----------------------------------------------------------------------------------------------------------------------------------------------------

   double zeta, sinphi, cosphi, lambda, kappa, arg1, arg2, calT, Bs_inf, Bp_inf, Bs, Bp, Bz, wus;
   bool small_s_eval = false;

   a = sqrt(a2);
   a2s2 = 1.0 - 2.0 / (r * (r + z));
   a1s1 = sqrt(a2s2);

// Separate treatment for z axis
   if (s < sp_tiny) {
      if (z > 1.0 + sp_tiny) {
         wus = sqrt(1.0 - 1.0 / Sqr(z));
         fields.Mag()[0] = B0[0] / wus;
         fields.Mag()[1] = B0[1] / wus;
         zeta = z + 0.5 * log((z - 1.0) / (z + 1.0));
      }
      else {

         if constexpr (HConfig::VLISMBochum_mod_type == 3 && HConfig::VLISMBochum_mod_rpos == 0) {
            // GetAmpFactor() could diverge, so we cannot use it right at the nose. Instead the amplification is computed by hand and "zeta" is artificially set to be larger than "ztr" so that GetAmpFactor() returns 1
            wus = 0.0;
            fields.Mag()[0] = B0[0] * scB / sin_theta_B0;
            fields.Mag()[1] = B0[1] * scB / sin_theta_B0;
            zeta = ztr + sp_small;
         }
         else {
// FIXME Derive asymptotics for other taper functions
         }

      };
      small_s_eval = true;
   }

// Compute the general solution
   else {

// Arguments to F and E
      lambda = sqrt(1.0 - a2s2);
      kappa = sqrt(1.0 + a2 / 4.0);
      arg1 = asin(lambda * kappa);
      arg2 = 1.0 / kappa;

// New implementation for circumventing divergence of F and E
      gsl_sf_result gsl_result_F, gsl_result_E;
      int statusF = gsl_sf_ellint_F_e(arg1, arg2, GSL_PREC_DOUBLE, &gsl_result_F);
      int statusE = gsl_sf_ellint_E_e(arg1, arg2, GSL_PREC_DOUBLE, &gsl_result_E);
      double ellintF, ellintE;

// A nonzero status indicates an error in elliptic integral evaluation. A value of 20 for F was obtained experimentally (the function is increasing logarithmically)
      if (statusF || statusE) {
         ellintF = 20.0;
         ellintE = 1.0;
      }
      else {
         ellintF = gsl_result_F.val;
         ellintE = gsl_result_E.val;
      };

// Compute caligraphy T
      calT = (2.0 * kappa - 1.0 / kappa) * ellintE - 2.0 * (kappa - 1.0 / kappa) * ellintF;

// Cylindrical components of the field at infinity for the same x and y as the point
      sinphi = posprime[1] / fmax(s, sp_tiny);
      cosphi = posprime[0] / fmax(s, sp_tiny);

      Bs_inf = B0[0] * cosphi + B0[1] * sinphi;
      Bp_inf = B0[1] * cosphi - B0[0] * sinphi;

// Transverse part of B.
      Bs = Bs_inf * (calT / a1s1 / r3 + a1s1 * (1.0 + z / r3));
      Bp = Bp_inf / a1s1;
      Bz = Bs_inf * (calT * (z / r3 - 1.0) + a2s2 * Sqr(z) / r3) / a;

// Build isochrone label. The points (s,z) and (inf,zeta) are on the same isochrone.
      zeta = z - 2.0 * lambda * sqrt((1.0 - Sqr(lambda * kappa)) / (1.0 - Sqr(lambda)))
             - (2.0 * kappa - 1.0 / kappa) * ellintF + 2.0 * kappa * ellintE;

// Transform components from cylindrical back to Cartesian.
      fields.Mag()[0] = Bs * cosphi - Bp * sinphi;
      fields.Mag()[1] = Bs * sinphi + Bp * cosphi;
      fields.Mag()[2] = Bz;
   };

// Apply isochrone-based scaling.
   fields.Mag() *= GetAmpFactor(zeta);

// Add unmodified longitudinal (flow-parallel) part.
   if (!small_s_eval) {
      fields.Mag()[0] -= B0[2] * posprime[0] / r3;
      fields.Mag()[1] -= B0[2] * posprime[1] / r3;
      fields.Mag()[2] -= B0[2] * (z / r3 - 1.0);
   }
   else fields.Mag()[2] = B0[2] * Sqr(wus);

// Convert back to global frame and compute the electric field.
   fields.Mag().ChangeFromBasis(eprime);


// --------------------------------------------------------------------------------
// TESTING: CAP THE MAGNETIC FIELD MAGNITUDE
// --------------------------------------------------------------------------------
   // double Bmaglim = 8.0e-6 / unit_magnetic_fluid;
   // _spdata.Bmag = _spdata.Bvec.Norm();
   // if (_spdata.Bmag > Bmaglim) _spdata.Bvec *= Bmaglim / _spdata.Bmag; 
// --------------------------------------------------------------------------------


// Note that the flags to compute U and B should be enabled in order to compute E
   fields.Elc() = -(fields.Vel() ^ fields.Mag()) / c_code;

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/17/2022
*/
template <typename HConfig>
template <typename Fields>
void BackgroundVLISMBochum<HConfig>::EvaluateBackgroundDerivatives(BackgroundCoordinates& coords, Fields& fields)
{
   NumericalDerivatives(coords, fields);
};

};
