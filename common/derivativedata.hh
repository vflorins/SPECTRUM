/*!
\file spatial_data.hh
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DERIVATIVE_DATA_HH
#define SPECTRUM_DERIVATIVE_DATA_HH

#include <common/vectors.hh>


namespace Spectrum {


//! Flag to indicate spatial derivatives were not computed
const uint16_t BACKGROUND_grad_FAIL = 0x1000;

//! Flag to indicate time derivatives were not computed
const uint16_t BACKGROUND_ddt_FAIL = 0x2000;



class DerivativeData {

public:

 //! "Safe" box for computing directional derivatives
   GeoVector _dr;

//! Spatial maximum distance per time step, grid dependent
   double dmax;

//! "Safe" time increment for computing time derivatives
   double _dt;

//! Status information for derivative computations
// todo: NOTE: this information is formerly stored in spdata._mask - need to review initialization
   uint16_t _status;

//! Flag for forward increment when computing directional derivatives
   bool _dr_forw_fail[3];

//! Flag for backward increment when computing directional derivatives
   bool _dr_back_fail[3];

//! Flag for forward increment when computing time derivatives
   bool _dt_forw_fail;

//! Flag for backward increment when computing time derivatives
   bool _dt_back_fail;

    DerivativeData() = default;

    DerivativeData& operator=(const DerivativeData& other) = default;

};


};


#endif
