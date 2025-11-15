/*!
\file physical_constants.hh
\brief Defines fundamental physical constants such as speed of light and electron charge
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_PHYSICAL_CONSTANTS_HH
#define SPECTRUM_PHYSICAL_CONSTANTS_HH

#include "config.h"

#ifdef USE_GSL
#include <gsl/gsl_const_cgsm.h>
#endif

namespace Spectrum {

//! Elementary charge (Fr)
#define SPC_CONST_CGSM_ELECTRON_CHARGE 4.8032044E-10

//! Speed of light (cm/s)
#ifdef USE_GSL
#define SPC_CONST_CGSM_SPEED_OF_LIGHT GSL_CONST_CGSM_SPEED_OF_LIGHT
#else
#define SPC_CONST_CGSM_SPEED_OF_LIGHT 2.99792458e+10
#endif

//! Electron mass (g)
#ifdef USE_GSL
#define SPC_CONST_CGSM_MASS_ELECTRON GSL_CONST_CGSM_MASS_ELECTRON
#else
#define SPC_CONST_CGSM_MASS_ELECTRON 9.1093837139E-28
#endif

//! Muon mass (g)
#ifdef USE_GSL
#define SPC_CONST_CGSM_MASS_MUON GSL_CONST_CGSM_MASS_MUON
#else
#define SPC_CONST_CGSM_MASS_MUON 1.883531627E-25
#endif

//! Proton mass (g)
#ifdef USE_GSL
#define SPC_CONST_CGSM_MASS_PROTON GSL_CONST_CGSM_MASS_PROTON
#else
#define SPC_CONST_CGSM_MASS_PROTON 1.67262192596E-24
#endif

//! Neutron mass (g)
#ifdef USE_GSL
#define SPC_CONST_CGSM_MASS_NEUTRON GSL_CONST_CGSM_MASS_NEUTRON
#else
#define SPC_CONST_CGSM_MASS_NEUTRON 1.67492750057E-24
#endif

//! Electronvolt (erg)
#ifdef USE_GSL
#define SPC_CONST_CGSM_ELECTRON_VOLT GSL_CONST_CGSM_ELECTRON_VOLT
#else
#define SPC_CONST_CGSM_ELECTRON_VOLT 1.602176634E-12
#endif

#define SPC_CONST_CGSM_KILO_ELECTRON_VOLT (1.0E3 * SPC_CONST_CGSM_ELECTRON_VOLT)
#define SPC_CONST_CGSM_MEGA_ELECTRON_VOLT (1.0E6 * SPC_CONST_CGSM_ELECTRON_VOLT)
#define SPC_CONST_CGSM_GIGA_ELECTRON_VOLT (1.0E9 * SPC_CONST_CGSM_ELECTRON_VOLT)

//! Boltzmann constant
#ifdef USE_GSL
#define SPC_CONST_CGSM_BOLTZMANN GSL_CONST_CGSM_BOLTZMANN
#else
#define SPC_CONST_CGSM_BOLTZMANN 1.380649E-16
#endif

//! Astronomical unit
#ifdef USE_GSL
#define SPC_CONST_CGSM_ASTRONOMICAL_UNIT GSL_CONST_CGSM_ASTRONOMICAL_UNIT
#else
#define SPC_CONST_CGSM_ASTRONOMICAL_UNIT 1.495978707E+13
#endif

//! Conversion of seconds into hours
#define HOURS(sec) ((sec) / 3600.0)

//! Conversion of seconds into days
#define DAYS(sec) ((sec) / 86400.0)

//! Conversion of seconds into years
#define YEARS(sec) ((sec) / 31536000.0)

};

#endif
