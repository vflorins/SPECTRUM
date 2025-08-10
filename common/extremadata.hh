/*!
\file spatial_data.hh
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_EXTREMA_DATA_HH
#define SPECTRUM_EXTREMA_DATA_HH


namespace Spectrum {


struct ExtremaData {

//! Maximum magnetic field magnitude along a trajectory
   double Bmag_max;

//! Minimum magnetic field magnitude along a trajectory
   double Bmag_min;

};


};


#endif
