/*!
\file background_vlism_bochum.hh
\brief Declares a plasma background class for the Very Local Interstellar Medium with the Roken/Kleimann analytic field model
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_VLISM_BOCHUM_HH
#define SPECTRUM_BACKGROUND_VLISM_BOCHUM_HH

#include "background_base.hh"

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
class BackgroundVLISMBochum : public BackgroundBase<HConfig_> {
private:

   //! Readable name of the class
   const std::string bg_name = "BackgroundVLISMBochum";

public:

   using HConfig = HConfig_;
   using BackgroundConfig = Cond<std::same_as<typename HConfig::BackgroundConfig, Default>, BackgroundDefault<BackgroundVLISMBochum<HConfig>>, typename HConfig::BackgroundConfig>;
   using BackgroundCoordinates = BackgroundConfig::Coordinates;
   using BackgroundBase = BackgroundBase<HConfig>;
   using BackgroundBase::_status;
   using BackgroundBase::container;
   using BackgroundBase::_ddata;
   using BackgroundBase::dmax0;
   using BackgroundBase::r0;
   using BackgroundBase::u0;
   using BackgroundBase::B0;
   // methods
   using BackgroundBase::EvaluateAbsMag;
   using BackgroundBase::EvaluateDmax;
   using BackgroundBase::GetDmax;
   using BackgroundBase::StopServerFront;
   using BackgroundBase::SetupBackground;

   using BackgroundConfig::derivative_method;
   using BackgroundConfig::mod_rpos;
   using BackgroundConfig::mod_type;

private:

   static_assert(!(mod_rpos != 0 && mod_rpos != 1), "Invalid mod_rpos");
   static_assert(!(mod_type == 1 && mod_rpos != 1), "Invalid combination of mod_rpos and mod_type");
   static_assert(!((mod_type == 2 || mod_type == 3) && mod_rpos != 0), "Invalid combination of mod_rpos and mod_type");

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

   static constexpr double ztr = Getztr();

//static constexpr double scB = 8.958 / 3.0; // 60 deg gives 8/3 ratio
//static constexpr double scB = 11.6 / 3.0; // 40 deg gives 8/3 ratio
   static constexpr double scB = 8.0 / 3.0; // 90 deg gives 8/3 ratio <- use this!
//static constexpr double scB = 30.0 / 3.0;


protected:

//! Distance to the heliopause in the nose direction (persistent)
   double z_nose;

//! Sine of the angle between u0 and B0 (persistent)
   double sin_theta_B0;

//! Amplification factor at the nose (persistent)
   double fzoom;

//! Local flow-aligned coordinate system (persistent)
   GeoVector eprime[3];

//! Strength of unmodified transversal B field at s=0 normalized to B0 for use in MOD_TYPE={2,3}.
   double RelBtrans(double z) const;

//! Returns the amplification factor for current isochrone
   double GetAmpFactor(double zeta) const;

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   void EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   void EvaluateBackgroundDerivatives(Coordinates&, Fields&);

public:

//! Default constructor
   BackgroundVLISMBochum(void);

//! Copy constructor
   BackgroundVLISMBochum(const BackgroundVLISMBochum& other);

//! Destructor
   ~BackgroundVLISMBochum() = default;

//! Clone function
   CloneFunctionBackground(BackgroundVLISMBochum);

};

};

// Something like this is needed for templated classes
#include "background_vlism_bochum.cc"

#endif
