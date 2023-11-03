/*!
\file turb_prop.hh
\brief Declares a structure storing turbulence properties
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TURB_PROP_HH
#define SPECTRUM_TURB_PROP_HH

namespace Spectrum {

//! Number of turbulence types
const int n_turb_types = 4;

//! Turbulence types that can be generated
enum turb_type {turb_alfven, turb_transverse, turb_longitudinal, turb_isotropic};

/*!
\brief Simple structure storing turbulence properties
\author Vladimir Florinski
\author Juan G Alonso Guzman
*/
struct TurbProp {

//! Smallest wavenumber
   double kmin;

//! Largest wavenumber
   double kmax;

//! Characteristic length
   double l0;

//! Number of modes
   double n_waves;

//! Variance
   double variance;

//! Power law slope
   double slope;
};

};

#endif
