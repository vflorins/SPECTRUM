/*!
\file background_solarwind.hh
\brief Declares a plasma background class for the constant speed supersonic wind of a rotating star
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SOLARWIND_HH
#define SPECTRUM_BACKGROUND_SOLARWIND_HH

#include "background_base.hh"

namespace Spectrum {

//! Method for computing derivatives (0: analytical, 1: numerical)
#define SOLARWIND_DERIVATIVE_METHOD 1

//! Heliospheric current sheet (0: disabled, 1: flat, 2: wavy (Jokipii-Thomas 1981) and static, 3: wavy and time-dependent).
#define SOLARWIND_CURRENT_SHEET 0

//! Magnetic topology region (0: nowhere, 1: same as HCS)
#define SOLARWIND_SECTORED_REGION 0

//! Correction to Parker Spiral, mainly for polar regions (0: none, 1: Smith-Bieber 1991, 2: Zurbuchen et al. 1997, 3: Schwadron-McComas 2003)
#define SOLARWIND_POLAR_CORRECTION 0

//! Latitudinal profile for bulk speed (0: constant, 1: linear step, 2: smooth step)
#define SOLARWIND_SPEED_LATITUDE_PROFILE 0

//! Heliopause radius
#ifdef USE_GSL
const double hp_rad_sw = 117.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
#else
const double hp_rad_sw = -1.0;
#endif

//! Magnetic axis tilt angle relative to the solar rotation axis
const double tilt_ang_sw = 40.0 * M_PI / 180.0;

#if SOLARWIND_CURRENT_SHEET == 3
//! Amplitude of variation to magnetic axis tilt angle
const double dtilt_ang_sw = 35.0 * M_PI / 180.0;

//! Solar cycle frequency
const double W0_sw = M_2PI / (60.0 * 60.0 * 24.0 * 365.0 * 22.0) / unit_frequency_fluid;

//! Factor to thin peaks and widen troughs (0.0: largest modification, 1.0: no modification)
const double stilt_ang_sw = 1.0;
#endif

#if SOLARWIND_POLAR_CORRECTION > 0
//! Differential rotation factor
const double delta_omega_sw = 0.05;
#endif

#if SOLARWIND_POLAR_CORRECTION == 2
//! Polar correction angle
const double polar_offset_sw = 30.0 * M_PI / 180.0;

//! Ratio of polar differential rotation to angular frequency of rotation
const double dwt_sw = delta_omega_sw * sin(polar_offset_sw);

//! Ratio of azimuthal differential rotation to angular frequency of rotation
const double dwp_sw = delta_omega_sw * cos(polar_offset_sw);
#endif

#if SOLARWIND_SPEED_LATITUDE_PROFILE > 0
//! Ratio of fast to slow wind speed
const double fast_slow_ratio_sw = 2.0;

//! Angle to transition between fast and slow
const double fast_slow_lat_sw = M_PI_2 - 30.0 * M_PI / 180.0 - tilt_ang_sw;

//! Transition speed coefficient
const double fast_slow_dlat_sw = 20.0;
#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSolarWind class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the class
const std::string bg_name_solarwind = "BackgroundSolarWind";

/*!
\brief Plasma background calculator for a radially expanding solar wind
\author Vladimir Florinski

Parameters: (BackgroundBase), GeoVector Omega, double r_ref, double dmax_fraction
*/
template <typename Fields_>
class BackgroundSolarWind : public BackgroundBase<Fields_> {
public:

   using Fields = Fields_;
   using BackgroundBase = BackgroundBase<Fields>;
   using BackgroundBase::_status;
   using BackgroundBase::_fields;
   using BackgroundBase::_ddata;
   using BackgroundBase::_pos;
   using BackgroundBase::_t;
   using BackgroundBase::container;
   using BackgroundBase::t0;
   using BackgroundBase::r0;
   using BackgroundBase::u0;
   using BackgroundBase::B0;
   using BackgroundBase::dmax0;
   // methods
   using BackgroundBase::EvaluateBmag;
   using BackgroundBase::EvaluateDmax;
   using BackgroundBase::StopServerFront;
   using BackgroundBase::SetupBackground;
   using BackgroundBase::EvaluateBackground;
   using BackgroundBase::EvaluateBackgroundDerivatives;
   using BackgroundBase::NumericalDerivatives;

protected:

//! Angular velocity vector of a rotating star (persistent)
   GeoVector Omega;

//! Reference radius (persistent)
   double r_ref;

//! Maximum fraction of the radial distance per step (persistent)
   double dmax_fraction;

//! Local coordinate system tied to the rotation axis (persistent)
   GeoVector eprime[3];

//! Velocity magnitude for slow wind (persistent)
   double ur0;

//! Radial magnetic field at "r_ref" (persistent)
   double Br0;

//! Angular frequency magnitude (persistent)
   double w0;

//! Position relative to origin (transient)
   GeoVector posprime;

#if SOLARWIND_SPEED_LATITUDE_PROFILE == 1
//! Latitude separating transition region from slow wind (persistent)
   double fsl_pls;

//! Latitude separating transition region from fast wind (persistent)
   double fsl_mns;
#elif SOLARWIND_SPEED_LATITUDE_PROFILE == 2
//! Half of fast-slow ratio plus 1 (persistent)
   double fsr_pls;

//! Half of fast-slow ratio minus 1 (persistent)
   double fsr_mns;
#endif

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct) override;

//! Modify radial flow (if necessary)
   virtual void ModifyUr(const double r, double &ur_mod);

//! Get time lag for time dependent current sheet (if necessary)
   virtual double TimeLag(const double r);

//! Compute the internal u, B, and E fields
   void EvaluateBackground(void) override;

//! Compute the internal u, B, and E derivatives
   void EvaluateBackgroundDerivatives(void) override;

//! Compute the maximum distance per time step
   void EvaluateDmax(void) override;

#if SOLARWIND_CURRENT_SHEET == 3
//! Function to compress peaks and stretch troughs
   double CubicStretch(double t) const;
#endif

public:

//! Default constructor
   BackgroundSolarWind(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundSolarWind(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor
   BackgroundSolarWind(const BackgroundSolarWind& other);

//! Destructor
   ~BackgroundSolarWind() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundSolarWind);
};

#if SOLARWIND_CURRENT_SHEET == 3
/*!
\author Juan G Alonso Guzman
\date 07/16/2024
\param[in] t periodic time to stretch between 0 and M_2PI
\return stretched time
*/
inline double BackgroundSolarWind::CubicStretch(double t) const
{
   double t_pi = t / M_PI;
   return t * ((1.0 - stilt_ang_sw) * t_pi * (t_pi - 3.0) + 3.0 - 2.0 * stilt_ang_sw);
};
#endif

};

// Something like this is needed for templated classes
#include "background_solarwind.cc"

#endif
