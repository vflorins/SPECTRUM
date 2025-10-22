/*!
\file diffusion.config.hh
\brief (Hyper)parameters and config(uration) options for a SPECTRUM diffusion class
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DIFFUSION_CONFIG_HH
#define SPECTRUM_DIFFUSION_CONFIG_HH


#include "common/compiletime_lists.hh"
#include "common/fields.hh"

namespace Spectrum {


/*!
\brief (Hyper)parameters and config(uration) options for a SPECTRUM diffusion class
\author Lucius Schoenbaum
\date 09/29/2025
*/
template <
      SpecieId specieid,
//! Flag to use QLT pitch angle scattering with WLNT perpendicular diffusion
      bool use_qlt_scatt_,
//! diffusion coefficient (if constant)
//! Scattering frequency (persistent)
      Ratio D0_,
      Ratio Dperp_,
      Ratio Dpara_,
      Ratio A2A_,
      Ratio l_max_,
      Ratio k_min_,
      Ratio ps_index_,
      Ratio ps_minus_,
      Ratio A2T_,
      Ratio A2L_,
      Ratio ps_plus_,
      Ratio l_max_HP_,
      Ratio dl_max_,
      Ratio nose_z_nose_,
      Ratio nose_z_sheath_,
      Ratio nose_dz_,
      Ratio kappa0_,
      Ratio kappa_ratio_,
      Ratio U0_,
      Ratio p0_,
      Ratio T0_,
      Ratio r0_,
      Ratio lam0_,
      Ratio R0_,
      Ratio B0_,
      Ratio pow_law_U_,
      Ratio pow_law_p_,
      Ratio pow_law_T_,
      Ratio pow_law_r_,
      Ratio pow_law_R_,
      Ratio pow_law_B_,
      int LISM_idx_,
      Ratio LISM_ind_,
      Ratio kappa_ratio_inner_,
      Ratio kappa_ratio_outer_,
      Ratio lam_inner_,
      Ratio lam_outer_,
      Ratio kappa_inner_,
      Ratio kappa_outer_,
      Ratio lam_para_,
      Ratio lam_perp_,
      Ratio kappa_ratio_red_,
      Ratio radial_limit_perp_red_,
      Ratio solar_cycle_idx_,
      Ratio solar_cycle_effect_,
      int Bmix_idx_,
      Ratio Bmix_ind_
>
struct DiffusionConfig {
   // todo reduced pitchangle ? (do not compute phi)
   using Coordinates = Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>;
   // todo Flum or AbsFlum in DiffusionFlowMomentumPowerLaw
   // todo some have indicator fields
   // todo implement: create default
   using Fields = Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>;

   static constexpr auto use_qlt_scatt = get_fp<use_qlt_scatt_>();
   static constexpr auto D0 = get_fp<D0_>();
   static constexpr auto Dperp = get_fp<Dperp_>();
   static constexpr auto Dpara = get_fp<Dpara_>();
   static constexpr auto A2A = get_fp<A2A_>();
   static constexpr auto l_max = get_fp<l_max_>();
   static constexpr auto k_min = get_fp<k_min_>();
   static constexpr auto ps_index = get_fp<ps_index_>();
   static constexpr auto ps_minus = get_fp<ps_minus_>();
   static constexpr auto A2T = get_fp<A2T_>();
   static constexpr auto A2L = get_fp<A2L_>();
   static constexpr auto ps_plus = get_fp<ps_plus_>();
   static constexpr auto l_max_HP = get_fp<l_max_HP_>();
   static constexpr auto dl_max = get_fp<dl_max_>();
   static constexpr auto nose_z_nose = get_fp<nose_z_nose_>();
   static constexpr auto nose_z_sheath = get_fp<nose_z_sheath_>();
   static constexpr auto nose_dz = get_fp<nose_dz_>();
   static constexpr auto kappa0 = get_fp<kappa0_>();
   static constexpr auto kappa_ratio = get_fp<kappa_ratio_>();
   static constexpr auto U0 = get_fp<U0_>();
   static constexpr auto p0 = get_fp<p0_>();
   static constexpr auto T0 = get_fp<T0_>();
   static constexpr auto r0 = get_fp<r0_>();
   static constexpr auto lam0 = get_fp<lam0_>();
   static constexpr auto R0 = get_fp<R0_>();
   static constexpr auto B0 = get_fp<B0_>();
   static constexpr auto pow_law_U = get_fp<pow_law_U_>();
   static constexpr auto pow_law_p = get_fp<pow_law_p_>();
   static constexpr auto pow_law_T = get_fp<pow_law_T_>();
   static constexpr auto pow_law_r = get_fp<pow_law_r_>();
   static constexpr auto poaw_law_R = get_fp<pow_law_R_>();
   static constexpr auto pow_law_B = get_fp<pow_law_B_>();
   static constexpr auto LISM_idx = LISM_idx_;
   static constexpr auto LISM_ind = get_fp<LISM_ind_>();
   static constexpr auto kappa_ratio_inner = get_fp<kappa_ratio_inner_>();
   static constexpr auto kappa_ratio_outer = get_fp<kappa_ratio_outer_>();
   static constexpr auto lam_inner = get_fp<lam_inner_>();
   static constexpr auto lam_outer = get_fp<lam_outer_>();
   static constexpr auto kappa_inner = get_fp<kappa_inner_>();
   static constexpr auto kappa_outer = get_fp<kappa_outer_>();
   static constexpr auto lam_para = get_fp<lam_para_>();
   static constexpr auto lam_perp = get_fp<lam_perp_>();
   static constexpr auto kappa_ratio_red = get_fp<kappa_ratio_red_>();
   static constexpr auto radial_limit_perp_red = get_fp<radial_limit_perp_red_>();
   static constexpr auto solar_cycle_idx = get_fp<solar_cycle_idx_>();
   static constexpr auto solar_cycle_effect = get_fp<solar_cycle_effect_>();
   static constexpr auto Bmix_idx = Bmix_idx_;
   static constexpr auto Bmix_ind = get_fp<Bmix_ind_>();
};



}


#endif