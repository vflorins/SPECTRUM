/*!
\file background_vlism_bochum.hh
\brief Declares a plasma background class for the Very Local Interstellar Medium with the Roken/Kleimann analytic field model
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef _BACKGROUND_VLISM_BOCHUM_HH
#define _BACKGROUND_VLISM_BOCHUM_HH

#include "background_base.hh"
#include <gsl/gsl_errno.h>

namespace Spectrum {

//! Readable name of the class
const std::string bg_name_bochum = "BackgroundVLISMBochum";

//! What function to use within 'get_ampfactor' (0 = none)
#define MOD_TYPE 3

//! Whether to scale relative to s=0 (0) or s=+inf (1)
#define MOD_RPOS 0

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

//const double scB = 8.958 / 3.0; // 60 deg gives 8/3 ratio
//const double scB = 11.6 / 3.0; // 40 deg gives 8/3 ratio
const double scB = 8.0 / 3.0; // 90 deg gives 8/3 ratio <- use this!

//const double scB = 30.0 / 3.0;

// Turn gsl error handler off
static const gsl_error_handler_t* gsl_default_error_handler = gsl_set_error_handler_off();

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundVLISMBochum class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma background calculator for the Roken-Kleimann model of field draping around the heliopause
\author Vladimir Florinski

This class calculates the velocity and magnetic fields around the heliopause, represented by a Rankine half-body (a product of an interaction between two potential flows, a point source and a uniform flow). The reference is: Roken, C., Kleimann, J., and Fichtner, H., An exact analytical solution for the interstellar magnetic field in the vicinity of the heliosphere, Astrophys. J., v. 805, p. 173 (2015).

Parameters: (BackgroundBase), double z_nose
*/
class BackgroundVLISMBochum : public BackgroundBase {

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
   void EvaluateBackground(void) override;

//! Compute the internal derivatives of the fields
   void EvaluateBackgroundDerivatives(void) override;

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

#endif
