/*!
\file spatial_data.hh
\brief Declares a structure storing spatial data
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SPATIAL_DATA_HH
#define SPECTRUM_SPATIAL_DATA_HH

#include <cstdint>
#include "vectors.hh"
#include "matrix.hh"

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Flags for background computation
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Basic quantities
const uint16_t BACKGROUND_U = 0x0001;
const uint16_t BACKGROUND_B = 0x0002;
const uint16_t BACKGROUND_E = 0x0004;
const uint16_t BACKGROUND_ALL = BACKGROUND_U | BACKGROUND_B | BACKGROUND_E;

//! Spatial derivatives
const uint16_t BACKGROUND_gradU = 0x0010;
const uint16_t BACKGROUND_gradB = 0x0020;
const uint16_t BACKGROUND_gradE = 0x0040;
const uint16_t BACKGROUND_gradALL = BACKGROUND_gradU | BACKGROUND_gradB | BACKGROUND_gradE;

//! Time derivatives
const uint16_t BACKGROUND_dUdt = 0x0100;
const uint16_t BACKGROUND_dBdt = 0x0200;
const uint16_t BACKGROUND_dEdt = 0x0400;
const uint16_t BACKGROUND_dALLdt = BACKGROUND_dUdt | BACKGROUND_dBdt | BACKGROUND_dEdt;

//! Flag to indicate spatial derivatives were not computed
const uint16_t BACKGROUND_grad_FAIL = 0x0800;

//! Flag to indicate time derivatives were not computed
const uint16_t BACKGROUND_ddt_FAIL = 0x1000;

//! Shift between mask blocks
const int mask_offset = 4;

namespace Spectrum {

/*!
\brief Simple structure storing physical data defined at a spatial location
\author Vladimir Florinski
\author Juan G Alonso Guzman
*/
struct SpatialData {

//! Mask for field copy
   uint16_t _mask;

//! Plasma velocity vector
   GeoVector Uvec;

//! Magnetic field vector
   GeoVector Bvec;

//! Electric field vector
   GeoVector Evec;

//! Plasma velocity gradient
   GeoMatrix gradUvec;

//! Magnetic field gradient
   GeoMatrix gradBvec;

//! Electric field gradient
   GeoMatrix gradEvec;

//! Plasma velocity gradient
   GeoVector dUvecdt;

//! Magnetic field gradient
   GeoVector dBvecdt;

//! Electric field gradient
   GeoVector dEvecdt;

//! Magnetic field magnitude
   double Bmag;

//! Minimum magnetic field magnitude along a trajectory
   double Bmag_min;

//! Maximum magnetic field magnitude along a trajectory
   double Bmag_max;

//! Magnetic field unit vector
   GeoVector bhat;

//! "Safe" box for computing directional derivatives
   GeoVector _dr;

//! "Safe" time increment for computing time derivatives
   double _dt;

//! Region
   double region;

//! Spatial maximum distance per time step, grid dependent
   double dmax;

//! Default constructor
   SpatialData(void) = default;

//! Assignment operator
   SpatialData& operator =(const SpatialData& other);

//! Divergence of Uvec
   double divU(void);

//! Divergence of Bvec
   double divB(void);

//! Divergence of Evec
   double divE(void);

//! Curl of Uvec
   GeoVector curlU(void);

//! Curl of Bvec
   GeoVector curlB(void);

//! Curl of Evec
   GeoVector curlE(void);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 02/10/2023
\param[in] other Structure to copy from
\return Reference to the object
*/
inline SpatialData& SpatialData::operator =(const SpatialData& other)
{
// Field copy
   if(BITS_RAISED(_mask, BACKGROUND_U)) Uvec = other.Uvec;
   if(BITS_RAISED(_mask, BACKGROUND_B)) {
      Bvec = other.Bvec;
      Bmag = other.Bmag;
      bhat = other.bhat;
   };
   if(BITS_RAISED(_mask, BACKGROUND_E)) Evec = other.Evec;

// Gradient copy
   if(BITS_RAISED(_mask, BACKGROUND_gradU)) gradUvec = other.gradUvec;
   if(BITS_RAISED(_mask, BACKGROUND_gradB)) gradBvec = other.gradBvec;
   if(BITS_RAISED(_mask, BACKGROUND_gradE)) gradEvec = other.gradEvec;
   if(BITS_RAISED(_mask, BACKGROUND_gradALL)) _dr = other._dr;

// Time derivative copy
   if(BITS_RAISED(_mask, BACKGROUND_dUdt)) dUvecdt = other.dUvecdt;
   if(BITS_RAISED(_mask, BACKGROUND_dBdt)) dBvecdt = other.dBvecdt;
   if(BITS_RAISED(_mask, BACKGROUND_dEdt)) dEvecdt = other.dEvecdt;
   if(BITS_RAISED(_mask, BACKGROUND_dALLdt)) _dt = other._dt;

   region = other.region;
   dmax = other.dmax;
   return *this;
};

/*!
\author Juan G Alonso Guzman
\date 10/18/2022
\return Divergence of U
*/
inline double SpatialData::divU(void)
{
   return gradUvec.Trace();
};

/*!
\author Juan G Alonso Guzman
\date 10/18/2022
\return Divergence of B
*/
inline double SpatialData::divB(void)
{
   return gradBvec.Trace();
};

/*!
\author Juan G Alonso Guzman
\date 10/18/2022
\return Divergence of E
*/
inline double SpatialData::divE(void)
{
   return gradEvec.Trace();
};

/*!
\author Juan G Alonso Guzman
\date 10/18/2022
\return Curl of U
*/
inline GeoVector SpatialData::curlU(void)
{
   GeoVector vec_tmp;
   vec_tmp[0] = gradUvec[2][1] - gradUvec[1][2]; 
   vec_tmp[1] = gradUvec[0][2] - gradUvec[2][0];
   vec_tmp[2] = gradUvec[1][0] - gradUvec[0][1];
   return vec_tmp;
};

/*!
\author Juan G Alonso Guzman
\date 10/18/2022
\return Curl of B
*/
inline GeoVector SpatialData::curlB(void)
{
   GeoVector vec_tmp;
   vec_tmp[0] = gradBvec[2][1] - gradBvec[1][2]; 
   vec_tmp[1] = gradBvec[0][2] - gradBvec[2][0];
   vec_tmp[2] = gradBvec[1][0] - gradBvec[0][1];
   return vec_tmp;
};

/*!
\author Juan G Alonso Guzman
\date 10/18/2022
\return Curl of E
*/
inline GeoVector SpatialData::curlE(void)
{
   GeoVector vec_tmp;
   vec_tmp[0] = gradEvec[2][1] - gradEvec[1][2]; 
   vec_tmp[1] = gradEvec[0][2] - gradEvec[2][0];
   vec_tmp[2] = gradEvec[1][0] - gradEvec[0][1];
   return vec_tmp;
};

};

#endif
