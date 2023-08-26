/*!
\file background_vlism_bochum.cc
\brief Implements a plasma background class for the Very Local Interstellar Medium with the Roken/Kleimann analytic field model
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_vlism_bochum.hh"
#include <gsl/gsl_sf_ellint.h>
#include <iostream>
#include <iomanip>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundVLISMBochum methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 09/28/2021
*/
BackgroundVLISMBochum::BackgroundVLISMBochum(void)
                     : BackgroundBase(bg_name_bochum, 0, MODEL_STATIC)
{
};

/*!
\author Jens Kleimann
\date 10/30/2019
\param[in] z Normalized z-component of position
\return Normalized transverse field
*/
double BackgroundVLISMBochum::RelBtrans(double z) const
{
   return 1.0 / sqrt(1.0 - 1.0 / Sqr(z));
};

/*!
\author Vladimir Florinski
\date 09/28/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundVLISMBochum::BackgroundVLISMBochum(const BackgroundVLISMBochum& other)
                     : BackgroundBase(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/26/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BackgroundVLISMBochum::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) BackgroundBase::SetupBackground(false);
   container.Read(&z_nose);

// Build the new coordinate system with z axis along -u0 and convert the magnetic field to the primed frame.
   eprime[2] = -UnitVec(u0);
   eprime[0] = GetSecondUnitVec(eprime[2]);
   eprime[1] = eprime[2] ^ eprime[0];

   u0.ChangeToBasis(eprime);
   B0.ChangeToBasis(eprime);

   sin_theta_B0 = fmax(tiny, sqrt(1.0 - Sqr(UnitVec(B0) * eprime[2])));

#if MOD_TYPE == 3
   double ht = Sqr(scB / sin_theta_B0) * (Sqr(ztr) - 1.0);
   fzoom = sqrt((ht - 1.0) / (ht - Sqr(ztr)));
#endif
};

/*!
\author Jens Kleimann
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/26/2022
\return Field amplification factor
*/
double BackgroundVLISMBochum::GetAmpFactor(double zeta) const
{
#if MOD_TYPE == 0
   return 1.0;
#endif

// Where to evaluate the custom amplification function.
   double zev;

// If "MOD_RPOS" is zero, the scaling is done with respect to the intersection point of the isochrone with the z-axis (s=0). To find "zev" the nonlinear equaiton zev+(1/2)ln[(zev-1)/(zev+1)]=zeta is solved via Newton iterations.
#if MOD_RPOS == 0

// Current error of iterator loop
   double f0, f1, delta = 1.0;

// Newton iterations
   zev = 2.0;
   while(fabs(delta) > tiny) {
      f0 = (zev + 0.5 * log((zev - 1.0) / (zev + 1.0)) - zeta);
      f1 = (zev * zev) / (zev * zev - 1.0); 
      delta = f0 / f1;

// Catch overshoot before doing log[<0] by just going halfway to 1.0
      if(zev - delta <= 1.0) delta = 0.5 * (zev - 1.0);
      zev -= delta;
   };

// If "MOD_RPOS" is one, the scaling is done with respect to the isochrone asymptotic position at s=inf, so zev=zeta.
#else
   zev = zeta;
#endif

// No modification for zev>ztr, zero for zev<ztr
#if MOD_TYPE == 1
   return (zev > ztr ? 1.0 : 0.0);

// No modification for zev>ztr, constant for zev<ztr
#elif MOD_TYPE == 2

   if(fabs(sin_theta_B0) < tiny) return 1.0;
   else return (zev > ztr ? 1.0 : RelBtrans(ztr) / RelBtrans(zev));

// No modification for zev>ztr, scaled to reach scB at the nose for zev<ztr
#elif MOD_TYPE == 3

   if((fabs(sin_theta_B0) < tiny) || (zev > ztr)) return 1.0;
   else return RelBtrans(fzoom * zev) / RelBtrans(fzoom * ztr) / RelBtrans(zev) * RelBtrans(ztr);

// Unmodified field
#else
   return 1.0;
#endif

};

/*!
\author Jens Kleimann
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/26/2022
*/
void BackgroundVLISMBochum::EvaluateBackground(void)
{
// Convert position into flow aligned coordinates and scale to "z_nose"
   GeoVector posprime = _pos - r0;
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
   if((z < 1.0) && ((s < tiny) || (a2 < 0.0))) {
      _spdata.Bvec = gv_zeros;
// No need to check computation flags since it's just assigning zeros
      _spdata.Evec = gv_zeros;
      _spdata.region = -1.0;

      RAISE_BITS(_status, STATE_INVALID);
      return;
   }
   else _spdata.region = 1.0;

// Compute the velocity (valid in every region). The velocity is transformed to the global frame and un-normalized.
   if(BITS_RAISED(_spdata._mask, BACKGROUND_U)) {
      _spdata.Uvec = (posprime / r3 - gv_nz) * u0.Norm();
      _spdata.Uvec.ChangeFromBasis(eprime);
   };

// If the flag for computing B is down, exit function because E depends on B, so it cannot be computed without it
   if(BITS_LOWERED(_spdata._mask, BACKGROUND_B)) return;

//----------------------------------------------------------------------------------------------------------------------------------------------------

   double zeta, sinphi, cosphi, lambda, kappa, arg1, arg2, ellintF, ellintE, calT, Bs_inf, Bp_inf, Bs, Bp, Bz, wus;
   bool small_s_eval = false;

   a = sqrt(a2);
   a2s2 = 1.0 - 2.0 / (r * (r + z));
   a1s1 = sqrt(a2s2);

// Separate treatment for z axis
   if(s < tiny) {
      if(z > 1.0 + tiny) {
         wus = sqrt(1.0 - 1.0 / Sqr(z));
         _spdata.Bvec[0] = B0[0] / wus;
         _spdata.Bvec[1] = B0[1] / wus;
         zeta = z + 0.5 * log((z - 1.0) / (z + 1.0));
      }
      else {

// FIXME Derive asymptotics for other taper functions
#if (MOD_RPOS == 0) && (MOD_TYPE == 3)

// GetAmpFactor() could diverge, so we cannot use it right at the nose. Instead the amplification is computed by hand and "zeta" is artificially set to be larger than "ztr" so that GetAmpFactor() returns 1
         wus = 0.0;
         _spdata.Bvec[0] = B0[0] * scB / sin_theta_B0;
         _spdata.Bvec[1] = B0[1] * scB / sin_theta_B0;
         zeta = ztr + small;

#endif

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
      if(statusF || statusE) {
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
      sinphi = posprime[1] / fmax(s, tiny);
      cosphi = posprime[0] / fmax(s, tiny);

      Bs_inf = B0[0] * cosphi + B0[1] * sinphi;
      Bp_inf = B0[1] * cosphi - B0[0] * sinphi;

// Transverse part of B.
      Bs = Bs_inf * (calT / a1s1 / r3 + a1s1 * (1.0 + z / r3));
      Bp = Bp_inf / a1s1;
      Bz = Bs_inf * (calT * (z / r3 - 1.0) + a2s2 * Sqr(z) / r3) / a;

// Build isochone label. The points (s,z) and (inf,zeta) are on the same isochrone.
      zeta = z - 2.0 * lambda * sqrt((1.0 - Sqr(lambda * kappa)) / (1.0 - Sqr(lambda)))
             - (2.0 * kappa - 1.0 / kappa) * ellintF + 2.0 * kappa * ellintE;

// Transform components from cylindrical back to Cartesian.
      _spdata.Bvec[0] = Bs * cosphi - Bp * sinphi;
      _spdata.Bvec[1] = Bs * sinphi + Bp * cosphi;
      _spdata.Bvec[2] = Bz;
   };

// Apply isochrone-based scaling.
   _spdata.Bvec *= GetAmpFactor(zeta);

// Add unmodified longitudinal (flow-parallel) part.
   if(!small_s_eval) {
      _spdata.Bvec[0] -= B0[2] * posprime[0] / r3;
      _spdata.Bvec[1] -= B0[2] * posprime[1] / r3;
      _spdata.Bvec[2] -= B0[2] * (z / r3 - 1.0);
   }
   else _spdata.Bvec[2] = B0[2] * Sqr(wus);

// Convert back to global frame and compute the electric field.
   _spdata.Bvec.ChangeFromBasis(eprime);


// --------------------------------------------------------------------------------
// TESTING: CAP THE MAGNETIC FIELD MAGNITUDE
// --------------------------------------------------------------------------------
   // double Bmaglim = 8.0e-6 / unit_magnetic_fluid;
   // _spdata.Bmag = _spdata.Bvec.Norm();
   // if (_spdata.Bmag > Bmaglim) _spdata.Bvec *= Bmaglim / _spdata.Bmag; 
// --------------------------------------------------------------------------------


// Note that the flags to compute U and B should be enabled in order to compute E
   if(BITS_RAISED(_spdata._mask, BACKGROUND_E)) _spdata.Evec = -(_spdata.Uvec ^ _spdata.Bvec) / c_code;

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/17/2022
*/
void BackgroundVLISMBochum::EvaluateBackgroundDerivatives(void)
{
   NumericalDerivatives();
};

};
