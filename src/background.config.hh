/*!
\file background.config.hh
\brief (Hyper)parameters and config(uration) options for a SPECTRUM MHD background class
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_CONFIG_HH
#define SPECTRUM_BACKGROUND_CONFIG_HH

#include "common/compiletime_lists.hh"
#include "common/fields.hh"

namespace Spectrum {

/*!
\brief (Hyper)parameters and config(uration) options for a SPECTRUM MHD background class
\author Lucius Schoenbaum
\date 09/29/2025
*/
template<
      SpecieId specieid,
      BackgroundOptions::DerivativeMethod derivative_method_,
      int num_numeric_grad_evals_,
      Ratio incr_dmax_ratio_,
// todo impl: SERVER_INTERP_ORDER
      int server_interpolation_order_,
      int smooth_discontinuity_order_,
      int server_num_ghost_cells_,
      BackgroundOptions::CurrentSheet solarwind_current_sheet_,
      BackgroundOptions::SectoredRegion solarwind_sectored_region_,
      BackgroundOptions::PolarCorrection solarwind_polar_correction_,
      BackgroundOptions::SpeedLatitudeProfile solarwind_speed_latitude_profile_,
      BackgroundOptions::TermShockSpeedExponent solarwind_termshock_speed_exponent_,
      BackgroundOptions::ModType mod_type_,
      BackgroundOptions::ModRPos mod_rpos_
>
struct BackgroundConfig {
/*
 * The coordinates for ALL backgrounds are: Pos_t, Time_t, with cartesian coordinates for position.
 * This information is stored as config info, but this is only done in order to make an update easier if this assumption changes.
 */
   using Coordinates = Fields<FConfig<specieid, CoordinateSystem::cartesian>, Time_t, Pos_t>;
/*
 * The following is boilerplate compile-time storage of the template arguments.
 */
   static constexpr auto derivative_method = derivative_method_;
   static constexpr auto num_numeric_grad_evals = num_numeric_grad_evals_;
   static constexpr auto incr_dmax_ratio = get_fp<incr_dmax_ratio_>();
   static constexpr auto server_interpolation_order = server_interpolation_order_;
   static constexpr auto smooth_discontinuity_order = smooth_discontinuity_order_;
   static constexpr auto server_num_ghost_cells = server_num_ghost_cells_;
   static constexpr auto solarwind_current_sheet = solarwind_current_sheet_;
   static constexpr auto solarwind_sectored_region = solarwind_sectored_region_;
   static constexpr auto solarwind_polar_correction = solarwind_polar_correction_;
   static constexpr auto solarwind_speed_latitude_profile = solarwind_speed_latitude_profile_;
   static constexpr auto solarwind_termshock_speed_exponent = solarwind_termshock_speed_exponent_;
   static constexpr auto mod_type = mod_type_;
   static constexpr auto mod_rpos = mod_rpos_;
};


}


#endif