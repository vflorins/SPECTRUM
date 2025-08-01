/*!
\file spatial_data.hh
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SPATIAL_DATA_HH
#define SPECTRUM_SPATIAL_DATA_HH

#include "mhdtuple/mhdtuple.hh"
#include <common/matrix.hh>

namespace Spectrum {


template <typename ... Ts>
class SpatialData: public MHDtuple<Ts...> {

   using MHDtuple = MHDtuple<Ts...>;

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

   //! Region
// TODO: provide for more regions (perhaps use a SimpleArray)
   GeoVector region;

//! Spatial maximum distance per time step, grid dependent
   double dmax;

    explicit SpatialData():
            MHDtuple()
    {}

   explicit SpatialData(Ts... in):
           MHDtuple(in...)
   {
      // nothing at this time.
   }

    SpatialData& operator=(const SpatialData& other) {
        MHDtuple::operator=(other);
        region = other.region;
        dmax = other.dmax;
        if (/* variables: D, d/dt */ true) {
           // Needed if and only if time derivatives are computed numerically
            _dt = other._dt;
            _dt_forw_fail = other._dt_forw_fail;
            _dt_back_fail = other._dt_back_fail;
        }
        if (/* variables: G, d/dx */ true) {
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

/*!
\author Juan G Alonso Guzman
\date 10/18/2022
\return Divergence of U
*/
    inline double divU(void)
    {
        return MHDtuple::Gvel().Trace();
    };

/*!
\author Juan G Alonso Guzman
\date 10/18/2022
\return Divergence of B
*/
    inline double divB(void)
    {
        return MHDtuple::GBdata().Gmag().Trace();
    };

/*!
\author Juan G Alonso Guzman
\date 10/18/2022
\return Divergence of E
*/
    inline double divE(void)
    {
        return MHDtuple::Gele().Trace();
    };

/*!
\author Juan G Alonso Guzman
\date 10/18/2022
\return Curl of U
*/
    inline GeoVector curlU(void)
    {
        GeoVector vec_tmp;
        GeoMatrix G = MHDtuple::Gvel();
        vec_tmp[0] = G[1][2] - G[2][1];
        vec_tmp[1] = G[2][0] - G[0][2];
        vec_tmp[2] = G[0][1] - G[1][0];
        return vec_tmp;
    };

/*!
\author Juan G Alonso Guzman
\date 10/18/2022
\return Curl of B
*/
    inline GeoVector curlB(void)
    {
        GeoVector vec_tmp;
        GeoMatrix G = MHDtuple::GBdata().Gmag();
        vec_tmp[0] = G[1][2] - G[2][1];
        vec_tmp[1] = G[2][0] - G[0][2];
        vec_tmp[2] = G[0][1] - G[1][0];
        return vec_tmp;
    };

/*!
\author Juan G Alonso Guzman
\date 10/18/2022
\return Curl of E
*/
    inline GeoVector curlE(void)
    {
        GeoVector vec_tmp;
        GeoMatrix G = MHDtuple::Gele();
        vec_tmp[0] = G[1][2] - G[2][1];
        vec_tmp[1] = G[2][0] - G[0][2];
        vec_tmp[2] = G[0][1] - G[1][0];
        return vec_tmp;
    };



/*!
\author Juan G Alonso Guzman
\date 07/02/2024
\return Divergence of bhat
\note The formula comes from applying vector identity (7) in the NRL Plasma formulary
*/
   double divbhat()
   {
      auto bhat = MHDtuple::Bdata().Dmag();
      auto Bmag = MHDtuple::Bdata().Mmag();
      auto gradBmag = MHDtuple::GBdata().GMmag();
      auto Bdiv = divB();
      double x1 = gradBmag * bhat;
      auto x2 = Bdiv - x1;
      auto x3 = x2/Bmag;
      return x3;
   };

/*!
\author Juan G Alonso Guzman
\date 07/02/2024
\return Curl of bhat
\note The formula comes from applying vector identity (8) in the NRL Plasma formulary
*/
   GeoVector curlbhat()
   {
      auto bhat = MHDtuple::Bdata().Dmag();
      // todo fix when using auto or MmagT to set type
      double Bmag = MHDtuple::Bdata().Mmag();
      auto gradBmag = MHDtuple::GBdata().GMmag();
      return (curlB() - (gradBmag ^ bhat)) / Bmag;
   };

/*!
\author Juan G Alonso Guzman
\date 07/02/2024
\return Gradient of bhat
\note The formula comes from expanding \partial_i bhat_j = d/dx^i (B_j / B)
*/
   GeoMatrix gradbhat()
   {
      auto bhat = MHDtuple::Bdata().Dmag();
      double Bmag = MHDtuple::Bdata().Mmag();
      auto gradB = MHDtuple::GBdata().Gmag();
      auto gradBmag = MHDtuple::GBdata().GMmag();
      // todo Dyadic can be made static
      GeoMatrix tmp;
      tmp.Dyadic(gradBmag, bhat);
      return (gradB - tmp) / Bmag;
   };

/*!
\author Juan G Alonso Guzman
\date 07/02/2024
\return Time derivative of bhat
*/
   GeoVector dbhatdt()
   {
      auto dBvecdt = MHDtuple::TBdata().Tmag();
      auto dBmagdt = MHDtuple::TBdata().TMmag();
      auto bhat = MHDtuple::Bdata().Dmag();
      double Bmag = MHDtuple::Bdata().Mmag();
      return (dBvecdt - (dBmagdt * bhat)) / Bmag;
   };

};


using SpatialDataU = SpatialData<vel_t>;
using SpatialDataUE = SpatialData<vel_t, ele_t>;
using SpatialDataUEB = SpatialData<vel_t, ele_t, Bdata_t>;
using SpatialDataALL = SpatialData<den_t,vel_t, ele_t, Bdata_t,Gvel_t, Gele_t, GBdata_t,Dvel_t, Dele_t, DBdata_t>;
//using SpatialDataGlobalBackground =

};


#endif
