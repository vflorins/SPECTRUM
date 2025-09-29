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
template <typename HConfig>
BackgroundSolarWind<HConfig>::BackgroundSolarWind(void)
                           : BackgroundBase(bg_name, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 02/22/2024
*/
template <typename HConfig>
BackgroundSolarWind<HConfig>::BackgroundSolarWind(const std::string& name_in, uint16_t status_in)
                           : BackgroundBase(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 01/26/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename HConfig>
BackgroundSolarWind<HConfig>::BackgroundSolarWind(const BackgroundSolarWind& other)
                           : BackgroundBase(other)
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
template <typename HConfig>
void BackgroundSolarWind<HConfig>::SetupBackground(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) BackgroundBase::SetupBackground(false);
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
   if constexpr (HConfig::solarwind_speed_latitude_profile == SolarWindOptions::SpeedLatitudeProfile::constant) {
      // do nothing.
   }
   if constexpr (HConfig::solarwind_speed_latitude_profile == SolarWindOptions::SpeedLatitudeProfile::linear_step) {
      fsl_pls = fast_slow_lat_sw + 2.0 / fast_slow_dlat_sw;
      fsl_mns = fast_slow_lat_sw - 2.0 / fast_slow_dlat_sw;
   }
   else if constexpr (HConfig::solarwind_speed_latitude_profile == SolarWindOptions::SpeedLatitudeProfile::smooth_step) {
      fsl_pls = 0.5 * (fast_slow_ratio_sw + 1.0);
      fsl_mns = 0.5 * (fast_slow_ratio_sw - 1.0);
   }
   else {
      // (undefined)
   }
};

/*!
\author Juan G Alonso Guzman
\date 03/14/2024
\param[in]  r      radial distance
\param[out] ur_mod modified radial flow
*/
template <typename HConfig>
void BackgroundSolarWind<HConfig>::ModifyUr(const double r, double &ur_mod)
{
};

/*!
\author Juan G Alonso Guzman
\date 06/21/2024
\param[in]  r radial distance
\param[out] time lag of propagation from solar surface to current position
*/
template <typename HConfig>
double BackgroundSolarWind<HConfig>::TimeLag(const double r)
{
   return r / ur0;
};

/*!
\author Vladimir Florinski
\date 06/21/2024
*/
template <typename HConfig>
template <typename Fields>
void BackgroundSolarWind<HConfig>::EvaluateBackground(BackgroundCoordinates& coords, Fields& fields)
{
   double r, s, costheta, sintheta, sinphi, cosphi;
   double r_mns, phase0, phase, sinphase, cosphase;
   double tilt_amp, t_lag, ur, Br, Bt, Bp, arg;

// Convert position into solar rotation frame
   posprime = coords.Pos() - r0;
   posprime.ChangeToBasis(eprime);
   r = posprime.Norm();
   r_mns = r - r_ref;

// Compute time lag due to solar wind propagation.
   t_lag = (coords.Time() - t0) - TimeLag(r);

// Find the magnitude of the magnetic field at the radial source surface accounting for the time lag. The value for B could be negative (for negative cycles).
   if constexpr (HConfig::SolarWindCurrentSheet == SolarWindOptions::CurrentSheet::wavy_time_dependent) {
      arg = 2.0 * W0_sw * t_lag;
      Br0 = B0[0] + B0[1] * cos(arg);
   }
   else {
      Br0 = B0[0];
   }
  
// Compute latitude and enforce equatorial symmetry.
   costheta = posprime[2] / r;
   double fs_theta_sym = acos(costheta);
   if (fs_theta_sym > M_PI_2) fs_theta_sym = M_PI - fs_theta_sym;

// indicator variables: region[0] = heliosphere(+1)/LISM(-1); region[1] = sectored(+1)/unipolar field(-1); region[2] = time-lagged solar cycle phase
   if constexpr (Fields::Iv0_found())
      fields.Iv0() = (r < hp_rad_sw ? 1.0 : -1.0);
   if constexpr (Fields::Iv1_found())
      fields.Iv1() = -1.0;
   if constexpr (Fields::Iv2_found())
      fields.Iv2() = arg;

// Assign magnetic mixing region
   if constexpr (HConfig::SolarWindCurrentSheet == SolarWindOptions::CurrentSheet::wavy_static
      || HConfig::SolarWindCurrentSheet == SolarWindOptions::CurrentSheet::wavy_time_dependent
   ) {
// Set tilt amplitude
      tilt_amp = tilt_ang_sw;
      if constexpr (HConfig::SolarWindCurrentSheet == SolarWindOptions::CurrentSheet::wavy_time_dependent) {
// Variable tilt
         tilt_amp += dtilt_ang_sw * cos(CubicStretch(arg - M_2PI * floor(arg / M_2PI)));
      }
      if constexpr (HConfig::SolarWindSectoredRegion == SolarWindOptions::SectoredRegion::HCS) {
         if (M_PI_2 - fs_theta_sym < tilt_amp) {
            // todo placeholder - review after spdata/fields update
//            _spdata.region[1] = 1.0;
            fields.Iv1() = 1.0;
         }
      }
   }

// Calculate speed (fast/slow) based on latitude
   ur = ur0;

   if constexpr (HConfig::solarwind_speed_latitude_profile == SolarWindOptions::SpeedLatitudeProfile::constant) {
      // do nothing.
   }
   if constexpr (HConfig::solarwind_speed_latitude_profile == SolarWindOptions::SpeedLatitudeProfile::linear_step) {
      if (fs_theta_sym < fsl_mns) ur *= fast_slow_ratio_sw;
      else if (fs_theta_sym < fsl_pls) ur *= fast_slow_ratio_sw - 0.25 * fast_slow_dlat_sw * (fs_theta_sym - fsl_mns);
   }
   else if constexpr (HConfig::solarwind_speed_latitude_profile == SolarWindOptions::SpeedLatitudeProfile::smooth_step) {
      ur *= fsl_pls - fsl_mns * tanh(fast_slow_dlat_sw * (fs_theta_sym - fast_slow_lat_sw));
   }
   else {
      // (undefined)
   }

// Modify "ur" with radial distance
   ModifyUr(r, ur);

// Compute the (radial) velocity and convert back to global frame
   if constexpr (Fields::Vel_found()) {
      fields.Vel() = ur * UnitVec(posprime);
      fields.Vel().ChangeFromBasis(eprime);
   };
   
// Compute (Parker spiral) magnetic field and convert back to global frame
   if constexpr (Fields::Mag_found()) {
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

      if constexpr (HConfig::SolarWindPolarCorrection == SolarWindOptions::PolarCorrection::none) {
         // do nothing.
      }
      else if constexpr (HConfig::SolarWindPolarCorrection == SolarWindOptions::PolarCorrection::Smith_Bieber) {
         Bp -= Br * delta_omega_sw * r * w0 / ur;
      }
      else if constexpr (HConfig::SolarWindPolarCorrection == SolarWindOptions::PolarCorrection::Zurbuchen_etal) {
         phase = r_mns * w0 / ur;
         phase0 = r_mns * w0 / ur0;
         sinphase = sin(phase0);
         cosphase = cos(phase0);
         Bt = Br * phase * dwt_sw * (sinphi * cosphase + cosphi * sinphase);
         Bp += Br * phase * (dwp_sw * sintheta + dwt_sw * costheta * (cosphi * cosphase - sinphi * sinphase));
      }
      else if (HConfig::SolarWindPolarCorrection == SolarWindOptions::PolarCorrection::Schwadron_McComas) {
         // TODO
      }

      fields.Mag()[0] = Br * sintheta * cosphi - Bp * sinphi;
      fields.Mag()[1] = Br * sintheta * sinphi + Bp * cosphi;
      fields.Mag()[2] = Br * costheta;


      if constexpr (HConfig::SolarWindPolarCorrection != SolarWindOptions::PolarCorrection::none
               && HConfig::SolarWindPolarCorrection != SolarWindOptions::PolarCorrection::Smith_Bieber
      ) {
// Add theta component from polar correction
         fields.Mag()[0] += Bt * costheta * cosphi;
         fields.Mag()[1] += Bt * costheta * sinphi;
         fields.Mag()[2] -= Bt * sintheta;
      }

// Correct polarity based on current sheet
      if constexpr (HConfig::SolarWindCurrentSheet == SolarWindOptions::CurrentSheet::disabled) {
         // do nothing.
      }
      else if constexpr (HConfig::SolarWindCurrentSheet == SolarWindOptions::CurrentSheet::flat) {
// Flat current sheet at the equator
         if (acos(costheta) > M_PI_2) fields.Mag() *= -1.0;
      }
      else {
// Wavy current sheet
         phase0 = w0 * t_lag;
         sinphase = sin(phase0);
         cosphase = cos(phase0);
         if (acos(costheta) > M_PI_2 + atan(tan(tilt_amp) * (sinphi * cosphase + cosphi * sinphase)))
            fields.Mag() *= -1.0;
         if constexpr (HConfig::SolarWindCurrentSheet == SolarWindOptions::CurrentSheet::wavy_time_dependent) {
// Solar cycle polarity changes
            if (sin(W0_sw * t_lag) < 0.0) fields.Mag() *= -1.0;
         }
      }
      fields.Mag().ChangeFromBasis(eprime);
   };

// Compute electric field, already in global frame. Note that the flags to compute U and B should be enabled in order to compute E.
   if constexpr (Fields::Elc_found()) fields.Elc() = -(fields.Vel() ^ fields.Vel()) / c_code;

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/14/2022
*/
template <typename HConfig>
template <typename Fields>
void BackgroundSolarWind<HConfig>::EvaluateBackgroundDerivatives(BackgroundCoordinates& coords, Fields& fields)
{
   if constexpr (HConfig::derivative_method == DerivativeMethod::analytic) {
      GeoMatrix rr;

      if constexpr (Fields::DelVel_found()) {
// Expression valid only for radial flow
         posprime = coords.Pos() - r0;
         rr.Dyadic(posprime);
         fields.DelVel() = (fields.Vel().Norm() / posprime.Norm()) * (gm_unit - rr);
      };
      if constexpr (Fields::DelMag_found()) {
//TODO: complete
      };
      if constexpr (Fields::DelElc_found()) {
         fields.DelElc() = -((fields.DelVel() ^ fields.Mag()) + (fields.Vel() ^ fields.DelMag())) / c_code;
      };
      if constexpr (Fields::DotVel_found()) fields.DotVel() = gv_zeros;
      if constexpr (Fields::DotMag_found()) fields.DotMag() = gv_zeros;
      if constexpr (Fields::DotElc_found()) fields.DotElc() = gv_zeros;
   }
   else {
      NumericalDerivatives(coords, fields);
   };
};

/*!
\author Vladimir Florinski
\date 03/10/2022
*/
template <typename HConfig>
void BackgroundSolarWind<HConfig>::EvaluateDmax(BackgroundCoordinates& coords)
{
   _ddata.dmax = fmin(dmax_fraction * (coords.Pos() - r0).Norm(), dmax0);
   LOWER_BITS(_status, STATE_INVALID);
};

};
