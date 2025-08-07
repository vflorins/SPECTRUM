/*!
\file background_solarwind.hh
\brief Implements a plasma background class for the constant speed supersonic wind of a rotating star
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_solarwind.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSolarWind methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/28/2022
*/
template <typename Fields>
BackgroundSolarWind<Fields>::BackgroundSolarWind(void)
                   : BackgroundBase<Fields>(bg_name_solarwind, 0, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 02/22/2024
*/
template <typename Fields>
BackgroundSolarWind<Fields>::BackgroundSolarWind(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                   : BackgroundBase<Fields>(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 01/26/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename Fields>
BackgroundSolarWind<Fields>::BackgroundSolarWind(const BackgroundSolarWind& other)
                   : BackgroundBase<Fields>(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\date 01/26/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename Fields>
void BackgroundSolarWind<Fields>::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BackgroundBase<Fields>::SetupBackground(false);
   container.Read(Omega);
   container.Read(r_ref);
   container.Read(dmax_fraction);

// Build the new coordinate system. The z axis is along "Omega" unless w0 = 0.0, in which case the system is non-rotating and the global z axis is used.
   w0 = Omega.Norm(); 
   if (w0 < sp_tiny) eprime[2] = gv_nz;
   else eprime[2] = UnitVec(Omega);
   eprime[0] = GetSecondUnitVec(eprime[2]);
   eprime[1] = eprime[2] ^ eprime[0];

// Only the first components of velocity is used.
   ur0 = fabs(u0[0]);

// Compute auxiliary quantities for fast-slow wind calculation
#if SOLARWIND_SPEED_LATITUDE_PROFILE == 1
   fsl_pls = fast_slow_lat_sw + 2.0 / fast_slow_dlat_sw;
   fsl_mns = fast_slow_lat_sw - 2.0 / fast_slow_dlat_sw;
#elif SOLARWIND_SPEED_LATITUDE_PROFILE == 2
   fsr_pls = 0.5 * (fast_slow_ratio_sw + 1.0);
   fsr_mns = 0.5 * (fast_slow_ratio_sw - 1.0);
#endif
};

/*!
\author Juan G Alonso Guzman
\date 03/14/2024
\param[in]  r      radial distance
\param[out] ur_mod modified radial flow
*/
template <typename Fields>
void BackgroundSolarWind<Fields>::ModifyUr(const double r, double &ur_mod)
{
};

/*!
\author Juan G Alonso Guzman
\date 06/21/2024
\param[in]  r radial distance
\param[out] time lag of propagation from solar surface to current position
*/
template <typename Fields>
double BackgroundSolarWind<Fields>::TimeLag(const double r)
{
   return r / ur0;
};

/*!
\author Vladimir Florinski
\date 06/21/2024
*/
template <typename Fields>
void BackgroundSolarWind<Fields>::EvaluateBackground(void)
{
   double r, s, costheta, sintheta, sinphi, cosphi;
   double r_mns, phase0, phase, sinphase, cosphase;
   double tilt_amp, t_lag, ur, Br, Bt, Bp, arg;

// Convert position into solar rotation frame
   posprime = _pos - r0;
   posprime.ChangeToBasis(eprime);
   r = posprime.Norm();
   r_mns = r - r_ref;

// Compute time lag due to solar wind propagation.
   t_lag = (_t - t0) - TimeLag(r);

// Find the magnitude of the magnetic field at the radial source surface accounting for the time lag. The value for B could be negative (for negative cycles).
#if SOLARWIND_CURRENT_SHEET == 3
   arg = 2.0 * W0_sw * t_lag;
   Br0 = B0[0] + B0[1] * cos(arg);
#else
   Br0 = B0[0];
#endif
  
// Compute latitude and enforce equatorial symmetry.
   costheta = posprime[2] / r;
   double fs_theta_sym = acos(costheta);
   if (fs_theta_sym > M_PI_2) fs_theta_sym = M_PI - fs_theta_sym;

// indicator variables: region[0] = heliosphere(+1)/LISM(-1); region[1] = sectored(+1)/unipolar field(-1); region[2] = time-lagged solar cycle phase
   _spdata.region[0] = (r < hp_rad_sw ? 1.0 : -1.0);
   _spdata.region[1] = -1.0;
   _spdata.region[2] = arg;
  
// Assign magnetic mixing region
#if SOLARWIND_CURRENT_SHEET >= 2
// Constant tilt
   tilt_amp = tilt_ang_sw;
#if SOLARWIND_CURRENT_SHEET == 3
// Variable tilt
   tilt_amp += dtilt_ang_sw * cos(CubicStretch(arg - M_2PI * floor(arg / M_2PI)));
#endif
#if SOLARWIND_SECTORED_REGION == 1
   if (M_PI_2 - fs_theta_sym < tilt_amp) _spdata.region[1] = 1.0;
#endif
#endif

// Calculate speed (fast/slow) based on latitude
   ur = ur0;
#if SOLARWIND_SPEED_LATITUDE_PROFILE == 1
      if (fs_theta_sym < fsl_mns) ur *= fast_slow_ratio_sw;
      else if (fs_theta_sym < fsl_pls) ur *= fast_slow_ratio_sw - 0.25 * fast_slow_dlat_sw * (fs_theta_sym - fsl_mns);
#elif SOLARWIND_SPEED_LATITUDE_PROFILE == 2
      ur *= fsr_pls - fsr_mns * tanh(fast_slow_dlat_sw * (fs_theta_sym - fast_slow_lat_sw));
#endif
// Modify "ur" with radial distance
   ModifyUr(r, ur);

// Compute the (radial) velocity and convert back to global frame
   if (BITS_RAISED(_spdata._mask, BACKGROUND_U)) {
      _spdata.Uvec = ur * UnitVec(posprime);
      _spdata.Uvec.ChangeFromBasis(eprime);
   };
   
// Compute (Parker spiral) magnetic field and convert back to global frame
   if (BITS_RAISED(_spdata._mask, BACKGROUND_B)) {
// Coordinates for conversion
      sintheta = sqrt(fmax(1.0 - Sqr(costheta), 0.0));
      s = sqrt(Sqr(posprime[0]) + Sqr(posprime[1]));
      if (s < r * sp_tiny) {
         cosphi = 0.0;
         sinphi = 0.0;
      }
      else {
         cosphi = posprime[0] / s;
         sinphi = posprime[1] / s;
      };
// Field components
      Br = Br0 * Sqr(r_ref / r);
      Bp = -Br * sintheta * r_mns * w0 / ur;
#if SOLARWIND_POLAR_CORRECTION == 1
      Bp -= Br * delta_omega_sw * r * w0 / ur;
#elif SOLARWIND_POLAR_CORRECTION == 2
      phase = r_mns * w0 / ur;
      phase0 = r_mns * w0 / ur0;
      sinphase = sin(phase0);
      cosphase = cos(phase0);
      Bt = Br * phase * dwt_sw * (sinphi * cosphase + cosphi * sinphase);
      Bp += Br * phase * (dwp_sw * sintheta + dwt_sw * costheta * (cosphi * cosphase - sinphi * sinphase));
#elif SOLARWIND_POLAR_CORRECTION == 3
// TODO
#endif

      _spdata.Bvec[0] = Br * sintheta * cosphi - Bp * sinphi;
      _spdata.Bvec[1] = Br * sintheta * sinphi + Bp * cosphi;
      _spdata.Bvec[2] = Br * costheta;
#if SOLARWIND_POLAR_CORRECTION > 1
// Add theta component from polar correction
      _spdata.Bvec[0] += Bt * costheta * cosphi;
      _spdata.Bvec[1] += Bt * costheta * sinphi;
      _spdata.Bvec[2] -= Bt * sintheta;
#endif

// Correct polarity based on current sheet
#if SOLARWIND_CURRENT_SHEET == 1
// Flat current sheet at the equator
      if (acos(costheta) > M_PI_2) _spdata.Bvec *= -1.0;
#elif SOLARWIND_CURRENT_SHEET >= 2
// Wavy current sheet
      phase0 = w0 * t_lag;
      sinphase = sin(phase0);
      cosphase = cos(phase0);
      if (acos(costheta) > M_PI_2 + atan(tan(tilt_amp) * (sinphi * cosphase + cosphi * sinphase))) _spdata.Bvec *= -1.0;
#if SOLARWIND_CURRENT_SHEET == 3
// Solar cycle polarity changes
      if (sin(W0_sw * t_lag) < 0.0) _spdata.Bvec *= -1.0;
#endif
#endif
      _spdata.Bvec.ChangeFromBasis(eprime);
   };

// Compute electric field, already in global frame. Note that the flags to compute U and B should be enabled in order to compute E.
   if (BITS_RAISED(_spdata._mask, BACKGROUND_E)) _spdata.Evec = -(_spdata.Uvec ^ _spdata.Bvec) / c_code;

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/14/2022
*/
template <typename Fields>
void BackgroundSolarWind<Fields>::EvaluateBackgroundDerivatives(void)
{
#if SOLARWIND_DERIVATIVE_METHOD == 0
   GeoVector posprime;
   GeoMatrix rr;

   if (BITS_RAISED(_spdata._mask, BACKGROUND_gradU)) {
// Expression valid only for radial flow
      posprime = _pos - r0;
      rr.Dyadic(posprime);
      _spdata.gradUvec = (_spdata.Uvec.Norm() / posprime.Norm()) * (gm_unit - rr);
   };
   if (BITS_RAISED(_spdata._mask, BACKGROUND_gradB)) {
//TODO: complete
   };
   if (BITS_RAISED(_spdata._mask, BACKGROUND_gradE)) {
      _spdata.gradEvec = -((_spdata.gradUvec ^ _spdata.Bvec) + (_spdata.Uvec ^ _spdata.gradBvec)) / c_code;
   };
   if (BITS_RAISED(_spdata._mask, BACKGROUND_dUdt)) _spdata.dUvecdt = gv_zeros;
   if (BITS_RAISED(_spdata._mask, BACKGROUND_dBdt)) _spdata.dBvecdt = gv_zeros;
   if (BITS_RAISED(_spdata._mask, BACKGROUND_dEdt)) _spdata.dEvecdt = gv_zeros;

#else
   NumericalDerivatives();
#endif
};

/*!
\author Vladimir Florinski
\date 03/10/2022
*/
template <typename Fields>
void BackgroundSolarWind<Fields>::EvaluateDmax(void)
{
   _spdata.dmax = fmin(dmax_fraction * (_pos - r0).Norm(), dmax0);
   LOWER_BITS(_status, STATE_INVALID);
};

};
