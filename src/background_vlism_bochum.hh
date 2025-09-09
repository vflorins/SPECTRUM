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
template <typename HyperParams_>
class BackgroundVLISMBochum : public BackgroundBase<HyperParams_> {
private:

   //! Readable name of the class
   const std::string bg_name = "BackgroundVLISMBochum";

#if (MOD_RPOS != 0) && (MOD_RPOS != 1)
#error Invalid MOD_RPOS
#endif

#if (MOD_TYPE == 1) && (MOD_RPOS != 1)
#error Invalid combination of MOD_RPOS and MOD_TYPE
#endif

#if ((MOD_TYPE == 2) || (MOD_TYPE == 3)) && (MOD_RPOS != 0)
#error Invalid combination of MOD_RPOS and MOD_TYPE
#endif

#if MOD_TYPE == 1
   const double ztr = -5.0;
#elif MOD_TYPE == 2
   const double ztr = 1.3;
#elif MOD_TYPE == 3
   const double ztr = 1.3;
#endif
   static constexpr double ztr = finish_me;

//static constexpr double scB = 8.958 / 3.0; // 60 deg gives 8/3 ratio
//static constexpr double scB = 11.6 / 3.0; // 40 deg gives 8/3 ratio
   static constexpr double scB = 8.0 / 3.0; // 90 deg gives 8/3 ratio <- use this!
//static constexpr double scB = 30.0 / 3.0;

public:

   using HyperParams = HyperParams_;
   using BackgroundBase = BackgroundBase<HyperParams>;
   using BackgroundBase::_status;
   using BackgroundBase::_fields;
   using BackgroundBase::_ddata;
   using BackgroundBase::_pos;
   using BackgroundBase::container;
   using BackgroundBase::r0;
   using BackgroundBase::u0;
   using BackgroundBase::B0;
   using BackgroundBase::dmax0;
   // methods
   using BackgroundBase::EvaluateBmag;
   using BackgroundBase::EvaluateDmax;
   using BackgroundBase::GetDmax;
   using BackgroundBase::StopServerFront;
   using BackgroundBase::SetupBackground;
//   using BackgroundBase::EvaluateBackground;
//   using BackgroundBase::EvaluateBackgroundDerivatives;
   using BackgroundBase::NumericalDerivatives;

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
   void SetupBackground(bool construct) override;

//! Compute the internal u, B, and E fields
   template <typename Fields>
   void EvaluateBackground(Fields&);

//! Compute the internal derivatives of the fields
   template <typename Fields>
   void EvaluateBackgroundDerivatives(Fields&);

public:

//! Default constructor
   BackgroundVLISMBochum(void);

//! Copy constructor
   BackgroundVLISMBochum(const BackgroundVLISMBochum& other);

//! Destructor
   ~BackgroundVLISMBochum() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundVLISMBochum);

};

};

// Something like this is needed for templated classes
#include "background_vlism_bochum.cc"

#endif
