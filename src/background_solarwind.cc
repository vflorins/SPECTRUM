/*!
\file background_solarwind.hh
\brief Implements a plasma background class for the constant speed supersonic wind of a rotating star
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_solarwind.hh"

namespace Spectrum {

using namespace BackgroundOptions;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSolarWind methods
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\author Vladimir Florinski
\date 01/26/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BackgroundSolarWind<HConfig>::SetupBackground(DataContainer& container)
{
// This is a backwards-compatible read-in of the container.
   container.Reset();
   container.Read(t0);
   container.Read(r0);
   container.Read(u0);
   container.Read(B0);
   container.Read(dmax0);
// end of BackgroundBase container read.
   container.Read(Omega);
   container.Read(r_ref);
   container.Read(dmax_fraction);

// Only the first components of velocity is used.
   ur0 = fabs(u0[0]);

// Build the new coordinate system. The z axis is along "Omega" unless w0 = 0.0, in which case the system is non-rotating and the global z axis is used.
   w0 = Omega.Norm(); 
   if (w0 < sp_tiny) eprime[2] = gv_nz;
   else eprime[2] = UnitVec(Omega);
   eprime[0] = GetSecondUnitVec(eprime[2]);
   eprime[1] = eprime[2] ^ eprime[0];

   if constexpr (with_termination_shock) {
      container.Read(TS.r);
      container.Read(TS.w);
      container.Read(TS.s);

      TS.s_inv = 1.0 / TS.s;
      TS.dmax = dmax_fraction * TS.w;
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
   if constexpr (with_termination_shock) {
      if (r > TS.r) {
         if (r > TS.r + TS.w) {
            if constexpr (termshock_speed_exponent == TermShockSpeedExponent::zero) {
               ur_mod *= TS.s_inv;
            }
            else if constexpr (termshock_speed_exponent == TermShockSpeedExponent::one) {
               ur_mod *= TS.s_inv * (TS.r + TS.w) / r;
            }
            else if constexpr (termshock_speed_exponent == TermShockSpeedExponent::square) {
               ur_mod *= TS.s_inv * Sqr((TS.r + TS.w) / r);
            }
            else {
// not implemented
               ur_mod *= 1e30;
            }
         }
         else {
            ur_mod *= 1.0 + (TS.s_inv - 1.0) * (r - TS.r) / TS.w;
         }
      };
   }
};

/*!
\author Juan G Alonso Guzman
\date 05/14/2025
\param[in]  r      radial distance
*/
template <typename HConfig>
double BackgroundSolarWind<HConfig>::dUrdr(const double r, const double v_norm)
{
   if constexpr (with_termination_shock) {
      if (r > TS.r) {
         if (r > TS.r + TS.w) {
            if constexpr (termshock_speed_exponent == TermShockSpeedExponent::zero) {
               return 0.0;
            }
            else if constexpr (termshock_speed_exponent == TermShockSpeedExponent::one) {
               return -v_norm / r;
            }
            else if constexpr (termshock_speed_exponent == TermShockSpeedExponent::square) {
               return -2.0 * v_norm / r;
            }
            else {
// not implemented
               return 1e30;
            }
         }
         else {
            return (TS.s_inv - 1.0) * (TS.r / TS.w) * ur0;
         }
      }
   }
   return 0.0;
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
   if constexpr (with_termination_shock) {
      if (r < TS.r) return r / ur0;
      else {
         if constexpr (termshock_speed_exponent == TermShockSpeedExponent::zero) {
            return (TS.r + TS.s * (r - TS.r)) / ur0;
         }
         else if constexpr (termshock_speed_exponent == TermShockSpeedExponent::one) {
            return (TS.r + TS.s * (Sqr(r) - Sqr(TS.r)) / (2.0 * TS.r)) / ur0;
         }
         else if constexpr (termshock_speed_exponent == TermShockSpeedExponent::square) {
            return (TS.r + TS.s * (Cube(r) - Cube(TS.r)) / (3.0 * Sqr(TS.r))) / ur0;
         }
         else {
// not implemented
            return 1e30;
         }
      }
   }
   else {
      return r / ur0;
   }
};

/*!
\author Vladimir Florinski
\date 06/21/2024
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundSolarWind<HConfig>::EvaluateBackground(Coordinates& coords, Fields& fields)
{
   double r, s, costheta, sintheta, sinphi, cosphi;
   double r_mns, phase0, phase, sinphase, cosphase;
   double tilt_amp, t_lag, ur, Br, Bt, Bp, arg;

// Convert position into solar rotation frame
   auto posprime = coords.Pos() - r0;
   posprime.ChangeToBasis(eprime);
   r = posprime.Norm();
   r_mns = r - r_ref;

// Compute time lag due to solar wind propagation.
   t_lag = (coords.Time() - t0) - TimeLag(r);

// Find the magnitude of the magnetic field at the radial source surface accounting for the time lag. The value for B could be negative (for negative cycles).
   double Br0;
   if constexpr (solarwind_current_sheet == CurrentSheet::wavy_time_dependent) {
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
   if constexpr (RequestedFields::Iv0_found())
      fields.Iv0('w') = (r < hp_rad_sw ? 1.0 : -1.0);
   if constexpr (RequestedFields::Iv1_found())
      fields.Iv1('w') = -1.0;
   if constexpr (RequestedFields::Iv2_found())
      fields.Iv2('w') = arg;

// Assign magnetic mixing region
   if constexpr (solarwind_current_sheet == CurrentSheet::wavy_static
      || solarwind_current_sheet == CurrentSheet::wavy_time_dependent
   ) {
// Set tilt amplitude
      tilt_amp = tilt_ang_sw;
      if constexpr (solarwind_current_sheet == CurrentSheet::wavy_time_dependent) {
// Variable tilt
         tilt_amp += dtilt_ang_sw * cos(CubicStretch(arg - M_2PI * floor(arg / M_2PI)));
      }
      if constexpr (solarwind_sectored_region == SectoredRegion::HCS) {
         if (M_PI_2 - fs_theta_sym < tilt_amp) {
            fields.Iv1('w') = 1.0;
         }
      }
   }

// Calculate speed (fast/slow) based on latitude
   ur = ur0;

   if constexpr (solarwind_speed_latitude_profile == SpeedLatitudeProfile::constant) {
      // do nothing.
   }
   if constexpr (solarwind_speed_latitude_profile == SpeedLatitudeProfile::linear_step) {
      if (fs_theta_sym < fsl_mns) ur *= fast_slow_ratio_sw;
      else if (fs_theta_sym < fsl_pls) ur *= fast_slow_ratio_sw - 0.25 * fast_slow_dlat_sw * (fs_theta_sym - fsl_mns);
   }
   else if constexpr (solarwind_speed_latitude_profile == SpeedLatitudeProfile::smooth_step) {
      ur *= fsl_pls - fsl_mns * tanh(fast_slow_dlat_sw * (fs_theta_sym - fast_slow_lat_sw));
   }
   else {
      // (undefined)
   }

// Modify "ur" with radial distance
   ModifyUr(r, ur);

// Compute the (radial) velocity and convert back to global frame
   if constexpr (RequestedFields::Fluv_found()) {
      fields.Fluv('w') = ur * UnitVec(posprime);
      fields.Fluv('w').ChangeFromBasis(eprime);
   };
   
// Compute (Parker spiral) magnetic field and convert back to global frame
   if constexpr (RequestedFields::Mag_found()) {
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

      if constexpr (solarwind_polar_correction == PolarCorrection::none) {
         // do nothing.
      }
      else if constexpr (solarwind_polar_correction == PolarCorrection::Smith_Bieber) {
         Bp -= Br * delta_omega_sw * r * w0 / ur;
      }
      else if constexpr (solarwind_polar_correction == PolarCorrection::Zurbuchen_etal) {
         phase = r_mns * w0 / ur;
         phase0 = r_mns * w0 / ur0;
         sinphase = sin(phase0);
         cosphase = cos(phase0);
         Bt = Br * phase * dwt_sw * (sinphi * cosphase + cosphi * sinphase);
         Bp += Br * phase * (dwp_sw * sintheta + dwt_sw * costheta * (cosphi * cosphase - sinphi * sinphase));
      }
      else if (solarwind_polar_correction == PolarCorrection::Schwadron_McComas) {
         // TODO
      }

      fields.Mag('w')[0] = Br * sintheta * cosphi - Bp * sinphi;
      fields.Mag('w')[1] = Br * sintheta * sinphi + Bp * cosphi;
      fields.Mag('w')[2] = Br * costheta;


      if constexpr (solarwind_polar_correction != PolarCorrection::none
               && solarwind_polar_correction != PolarCorrection::Smith_Bieber
      ) {
// Add theta component from polar correction
         fields.Mag('w')[0] += Bt * costheta * cosphi;
         fields.Mag('w')[1] += Bt * costheta * sinphi;
         fields.Mag('w')[2] -= Bt * sintheta;
      }

// Correct polarity based on current sheet
      if constexpr (solarwind_current_sheet == CurrentSheet::disabled) {
         // do nothing.
      }
      else if constexpr (solarwind_current_sheet == CurrentSheet::flat) {
// Flat current sheet at the equator
         if (acos(costheta) > M_PI_2) fields.Mag('w') *= -1.0;
      }
      else {
// Wavy current sheet
         phase0 = w0 * t_lag;
         sinphase = sin(phase0);
         cosphase = cos(phase0);
         if (acos(costheta) > M_PI_2 + atan(tan(tilt_amp) * (sinphi * cosphase + cosphi * sinphase)))
            fields.Mag('w') *= -1.0;
         if constexpr (solarwind_current_sheet == CurrentSheet::wavy_time_dependent) {
// Solar cycle polarity changes
            if (sin(W0_sw * t_lag) < 0.0) fields.Mag('w') *= -1.0;
         }
      }
      fields.Mag('w').ChangeFromBasis(eprime);
   };

// Compute electric field, already in global frame. Note that the flags to compute U and B should be enabled in order to compute E.
   if constexpr (RequestedFields::Ele_found()) fields.Ele('w') = -(fields.Fluv() ^ fields.Fluv()) / c_code;

   return 0;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/14/2022
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundSolarWind<HConfig>::EvaluateBackgroundDerivatives(Coordinates& coords, Fields& fields)
{
   if constexpr (RequestedFields::DelFluv_found()) {
// Expression valid only for radial flow
      GeoVector posprime = coords.Pos() - r0;
      if constexpr (with_termination_shock) {
         double r = posprime.Norm();
         double v_norm = fields.AbsFluv();
         GeoMatrix rr = Dyadic(posprime / r);
         fields.DelFluv('w') = dUrdr(r, v_norm) * rr + (v_norm / r) * (gm_unit - rr);
      }
      else {
         GeoMatrix rr = Dyadic(posprime);
         fields.DelFluv('w') = (fields.AbsFluv() / posprime.Norm()) * (gm_unit - rr);
      }
   };
   if constexpr (RequestedFields::DelMag_found()) {
//TODO: complete
   };
   if constexpr (RequestedFields::DelEle_found()) {
      fields.DelEle('w') = -((fields.DelFluv() ^ fields.Mag()) + (fields.Fluv() ^ fields.DelMag())) / c_code;
   };
   if constexpr (RequestedFields::DotFluv_found()) fields.DotFluv('w') = gv_zeros;
   if constexpr (RequestedFields::DotMag_found()) fields.DotMag('w') = gv_zeros;
   if constexpr (RequestedFields::DotEle_found()) fields.DotEle('w') = gv_zeros;
   return 0;
};

/*!
\author Vladimir Florinski
\date 03/10/2022
*/
template <typename HConfig>
template <typename Coordinates>
status_t BackgroundSolarWind<HConfig>::EvaluateDmax(Coordinates& coords, double* dmax)
{
   *dmax = fmin(dmax_fraction * (coords.Pos() - r0).Norm(), dmax0);
   if constexpr (with_termination_shock) {
// Reduce "dmax" around the shock. This implemenation assumes that "dmax" = "dmax0" near "TS.r" by default.
      double r = (coords.Pos() - r0).Norm();
      if (TS.r - dmax0 < r && r < TS.r + TS.w + dmax0) {
         if (r < TS.r) *dmax += (TS.dmax - dmax0) * (r - TS.r + dmax0) / dmax0;
         else if (r > TS.r + TS.w) *dmax -= (TS.dmax - dmax0) * (r - TS.r - TS.w - dmax0) / dmax0;
         else *dmax = TS.dmax;
      };
   }
   return 0;
};

};

