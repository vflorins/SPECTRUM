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


template <typename ... Ts>
class DerivativeData {

public:

    //! "Safe" box for computing directional derivatives
   GeoVector _dr;

//! Flag for forward increment when computing directional derivatives
   bool _dr_forw_fail[3];

//! Flag for backward increment when computing directional derivatives
   bool _dr_back_fail[3];

//! "Safe" time increment for computing time derivatives
   double _dt;

//! Flag for forward increment when computing time derivatives
   bool _dt_forw_fail;

//! Flag for backward increment when computing time derivatives
   bool _dt_back_fail;

//! Spatial maximum distance per time step, grid dependent
   double dmax;

    DerivativeData() = default;

    DerivativeData& operator=(const DerivativeData& other) {
        dmax = other.dmax;
        if (/* variables: Der, d/dt */ true) {
           // Needed if and only if time derivatives are computed numerically
            _dt = other._dt;
            _dt_forw_fail = other._dt_forw_fail;
            _dt_back_fail = other._dt_back_fail;
        }
        if (/* variables: Del, d/dx */ true) {
           // Needed if and only if spatial derivatives (gradients) are computed numerically
           // cf. trajectory_parker.cc
            _dr = other._dr;
            for (int xyz = 0; xyz < 3; ++xyz) {
                _dr_forw_fail[xyz] = other._dr_forw_fail[xyz];
                _dr_back_fail[xyz] = other._dr_back_fail[xyz];
            }
        }
        return *this;
    }


};


};


#endif
