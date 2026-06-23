/*!
\file background_vlism_bochum.hh
\brief Declares a plasma background class for the Very Local Interstellar Medium with the Roken/Kleimann analytic field model
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_VLISM_BOCHUM_HH
#define SPECTRUM_BACKGROUND_VLISM_BOCHUM_HH

#include "common/vectors.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundVLISMBochum class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma background calculator for the Roken-Kleimann model of field draping around the heliopause
\author Vladimir Florinski

This class calculates the velocity and magnetic fields around the heliopause, represented by a Rankine half-body (a product of an interaction between two potential flows, a point source and a uniform flow). The reference is: Roken, C., Kleimann, J., and Fichtner, H., An exact analytical solution for the interstellar magnetic field in the vicinity of the heliosphere, Astrophys. J., v. 805, p. 173 (2015).

Parameters: (BackgroundBase), double z_nose
*/
template <typename HConfig_>
class BackgroundVLISMBochum {
public:

   //! Readable name of the class
   const std::string name = "BackgroundVLISMBochum";

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;

// secular config:
   static constexpr bool requires_setup = true;
   static constexpr bool stochastic = false;

   static constexpr double dmax0 = Config::dmax0;

   static constexpr auto mod_rpos = Config::mod_rpos;
   static constexpr auto mod_type = Config::mod_type;

private:

   static_assert(!(mod_rpos != ModRPos::scale_rel_zero && mod_rpos != ModRPos::scale_rel_inf), "Invalid mod_rpos");
   static_assert(!(mod_type == ModType::zero && mod_rpos != ModRPos::scale_rel_inf), "Invalid combination of mod_rpos and mod_type");
   static_assert(!((mod_type == ModType::constant || mod_type == ModType::scaled) && mod_rpos != ModRPos::scale_rel_zero), "Invalid combination of mod_rpos and mod_type");

   static constexpr double Getztr() {
      if constexpr (mod_type == 1)
         return -5.0;
      else if constexpr (mod_type == 2)
         return 1.3;
      else if constexpr (mod_type == 3)
         return 1.3;
      else
         return 1.0;
   }

/*!
\author Jens Kleimann
\date 10/30/2019
\param[in] z Normalized z-component of position
\return Normalized transverse field
*/
//! Strength of unmodified transversal B field at s=0 normalized to B0 for use in MOD_TYPE={2,3}.
   static constexpr double RelBtrans(double z)
   {
      return 1.0 / csqrt(1.0 - 1.0 / Sqr(z));
   };

   static constexpr double construct_sin_theta_B0() {
      constexpr double tmp = csqrt(1.0 - Sqr(UnitVec(B0) * eprime[2]));
      if constexpr (sp_tiny > tmp)
         return tmp;
      else
         return sp_tiny;
   }

   static constexpr double construct_fzoom() {
      if constexpr (mod_type == 3) {
         double ht = Sqr(scB / sin_theta_B0) * (Sqr(ztr) - 1.0);
         return csqrt((ht - 1.0) / (ht - Sqr(ztr)));
      }
      else
         return 0;
   }

   static constexpr std::array<GeoVector, 3> construct_flowaligned_cs() {
      std::array<GeoVector, 3> out;
      out[2] = -UnitVec(Config::u0);
      out[0] = GetSecondUnitVec(eprime[2]);
      out[1] = out[2] ^ out[0];
      return out;
   }

   static constexpr GeoVector construct_u0() {
      GeoVector out = Config::u0;
      out.ChangeToBasis(eprime.data());
      return out;
   }

//static constexpr double scB = 8.958 / 3.0; // 60 deg gives 8/3 ratio
//static constexpr double scB = 11.6 / 3.0; // 40 deg gives 8/3 ratio
   static constexpr double scB = 8.0 / 3.0; // 90 deg gives 8/3 ratio <- use this!
//static constexpr double scB = 30.0 / 3.0;

public:

   static constexpr GeoVector r0 = Config::r0;

   static constexpr GeoVector B0 = Config::B0;

   //! Distance to the heliopause in the nose direction
   static constexpr double z_nose = Config::z_nose;

//! todo docstring
   static constexpr double ztr = Getztr();

//! Amplification factor at the nose
   static constexpr double fzoom = construct_fzoom();

//! Local flow-aligned coordinate system (persistent)
   static constexpr std::array<GeoVector, 3> eprime = construct_flowaligned_cs();

   static constexpr GeoVector u0 = construct_u0();

   //! Sine of the angle between u0 and B0 (persistent)
   static constexpr double sin_theta_B0 = construct_sin_theta_B0();

//! Returns the amplification factor for current isochrone
   static double GetAmpFactor(double zeta);

//! Set up the field evaluator based on "params"
   void SetupBackground(DataContainer& container_in);

//! Compute the maximum distance per time step
   template <typename Coordinates>
   static status_t EvaluateDmax(Coordinates& coords, double*);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackgroundDerivatives(Coordinates&, Fields&);

};

};

// Something like this is needed for templated classes
#include "background_vlism_bochum.cc"

#endif
