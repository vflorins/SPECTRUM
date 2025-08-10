/*!
\file server_server.cc
\brief Implements a class of a data server for a uniform Cartesian grid
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "server_cartesian.hh"
#include "block_cartesian.hh"
#include "reader_cartesian.hh"
#include "common/print_warn.hh"
#include <iostream>
#include <iomanip>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ServerCartesianFront templated methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 12/01/2023
\param[out] spdata Fields, dmax, etc.
*/
template <typename Fields>
void ServerCartesianFront::GetVariablesFromReader(Fields& fields)
{
   int xyz;
   double rho, vars[n_variables] = {0.0};

// Get variables
   MPI_Send(&_inquiry, 1, MPIInquiryType, 0, tag_needvars, MPI_Config::node_comm);
   MPI_Recv(vars, n_variables, MPI_DOUBLE, 0, tag_sendvars, MPI_Config::node_comm, MPI_STATUS_IGNORE);
   stencil_outcomes[2]++;

// Mass density, if provided
#ifdef SERVER_VAR_INDEX_RHO
   rho = vars[SERVER_VAR_INDEX_RHO];
#endif

// Number density, if provided
#ifdef SERVER_VAR_INDEX_DEN
   spdata.n_dens = vars[SERVER_VAR_INDEX_DEN];
#endif

// Thermal pressure, if provided
#ifdef SERVER_VAR_INDEX_PTH
   spdata.p_ther = vars[SERVER_VAR_INDEX_PTH];
#endif

   for (xyz = 0; xyz < 3; xyz++) {

// Bulk flow from mass density and momentum, if provided
#if defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO)
      spdata.Uvec[xyz] = vars[SERVER_VAR_INDEX_MOM + xyz] / rho;
// Bulk flow, if provided
#elif defined(SERVER_VAR_INDEX_FLO)
      spdata.Uvec[xyz] = vars[SERVER_VAR_INDEX_FLO + xyz];
#else
      spdata.Uvec[xyz] = 0.0;
#endif

// The magnetic field must be always provided
      spdata.Bvec[xyz] = vars[SERVER_VAR_INDEX_MAG + xyz];

// Electric field, if provided
#if defined(SERVER_VAR_INDEX_ELE)
      spdata.Evec[xyz] = vars[SERVER_VAR_INDEX_ELE + xyz];
#elif defined(SERVER_VAR_INDEX_FLO) && !(defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO)) 
      spdata.Evec[xyz] = 0.0;
#endif

   };

// Electric field, if B and U provided
#ifndef SERVER_VAR_INDEX_ELE
#if defined(SERVER_VAR_INDEX_FLO) || (defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))   
   spdata.Evec = -(spdata.Uvec ^ spdata.Bvec) / c_code;
#endif
#endif

// Region(s) indicator variable(s), if provided
#ifdef SERVER_VAR_INDEX_REG
   for (xyz = 0; xyz < SERVER_NUM_INDEX_REG; xyz++) spdata.region[xyz] = vars[SERVER_VAR_INDEX_REG + xyz];
#else
   spdata.region = gv_zeros;
#endif
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 12/01/2023
\param[in]  pos    Position
\param[in]  block  Block containing pos
\param[out] spdata Fields, dmax, etc.
*/
tempate <typename Fields>
void ServerCartesianFront::GetVariablesInterp0(const GeoVector& pos, Fields& fields)
{
   int xyz;
   double rho;

// Take the nearest cell value
   MultiIndex zone = block_pri->GetZone(pos);

// Mass density, if provided
#ifdef SERVER_VAR_INDEX_RHO
   rho = block_pri->GetValue(zone, SERVER_VAR_INDEX_RHO);
#endif

// Number density, if provided
#ifdef SERVER_VAR_INDEX_DEN
   spdata.n_dens = block_pri->GetValue(zone, SERVER_VAR_INDEX_DEN);
#endif

// Thermal pressure, if provided
#ifdef SERVER_VAR_INDEX_PTH
   spdata.p_ther = block_pri->GetValue(zone, SERVER_VAR_INDEX_PTH);
#endif

// Convert the variables to SPECTRUM format
   for (xyz = 0; xyz < 3; xyz++) {

// Bulk flow from mass density and momentum, if provided
#if defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO)
      spdata.Uvec[xyz] = block_pri->GetValue(zone, SERVER_VAR_INDEX_MOM + xyz) / rho;
// Bulk flow, if provided
#elif defined(SERVER_VAR_INDEX_FLO)
      spdata.Uvec[xyz] = block_pri->GetValue(zone, SERVER_VAR_INDEX_FLO + xyz);
#else
      spdata.Uvec[xyz] = 0.0;
#endif

// The magnetic field must be always provided
      spdata.Bvec[xyz] = block_pri->GetValue(zone, SERVER_VAR_INDEX_MAG + xyz);

// Electric field, if provided
#if defined(SERVER_VAR_INDEX_ELE)
      spdata.Evec[xyz] = block_pri->GetValue(zone, SERVER_VAR_INDEX_ELE + xyz);
#elif !defined(SERVER_VAR_INDEX_FLO) && !(defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
      spdata.Evec[xyz] = 0.0;
#endif

   };

// Electric field, if B and U provided
#ifndef SERVER_VAR_INDEX_ELE
#if defined(SERVER_VAR_INDEX_FLO) || (defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))   
   spdata.Evec = -(spdata.Uvec ^ spdata.Bvec) / c_code;
#endif
#endif

// Region(s) indicator variable(s), if provided
#ifdef SERVER_VAR_INDEX_REG
   for (xyz = 0; xyz < SERVER_NUM_INDEX_REG; xyz++) spdata.region[xyz] = block_pri->GetValue(zone, SERVER_VAR_INDEX_REG + xyz);
#else
   spdata.region = gv_zeros;
#endif
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 12/01/2023
\param[in]  pos    Position
\param[out] spdata Fields, dmax, etc.
*/
template <typename Fields>
void ServerCartesianFront::GetVariablesInterp1(const GeoVector& pos, Fields& fields)
{
   int xyz, iz, vidx, pri_idx, sec_idx;
   double rho, var, vars[n_variables] = {0.0}, _Bmag2;

// Build the stencil. This is a time consuming operation.
   stencil_status = BuildInterpolationStencil(pos);
   spdata.Bmag = 0.0;

// Internal interpolation
   if (stencil_status == 0) {
      stencil_outcomes[0]++;
      for (iz = 0; iz < stencil.n_elements; iz++) {
         for (vidx = 0; vidx < n_variables; vidx++) {
            var = block_pri->GetValue(stencil.zones[iz], vidx);
            vars[vidx] += stencil.weights[iz] * var;
         };
// B magnitude
         _Bmag2 = Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                + Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                + Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         spdata.Bmag += stencil.weights[iz] * sqrt(_Bmag2);
      };
   }

// Edge or corner interpolation
   else if (stencil_status == 4) {
      stencil_outcomes[2]++;
      pri_idx = block_pri->GetNode();
      sec_idx = block_sec->GetNode();
      for (iz = 0; iz < stencil.n_elements; iz++) {
         if (stencil.blocks[iz] == pri_idx) block_stn = block_pri;
         else if (stencil.blocks[iz] == sec_idx) block_stn = block_sec;
         else block_stn = cache_line[stencil.blocks[iz]];
         for (vidx = 0; vidx < n_variables; vidx++) {
            var = block_stn->GetValue(stencil.zones[iz], vidx);
            vars[vidx] += stencil.weights[iz] * var;
         };
// B magnitude
         _Bmag2 = Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         spdata.Bmag += stencil.weights[iz] * sqrt(_Bmag2);
      };
   }

// Plane interpolation
   else {
      stencil_outcomes[1]++;
      pri_idx = block_pri->GetNode();
      for (iz = 0; iz < stencil.n_elements; iz++) {
         if (stencil.blocks[iz] == pri_idx) block_stn = block_pri;
         else block_stn = block_sec;
         for (vidx = 0; vidx < n_variables; vidx++) {
            var = block_stn->GetValue(stencil.zones[iz], vidx);
            vars[vidx] += stencil.weights[iz] * var;
         };
// B magnitude
         _Bmag2 = Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         spdata.Bmag += stencil.weights[iz] * sqrt(_Bmag2);
      };
   };

// Mass density, if provided
#ifdef SERVER_VAR_INDEX_RHO
   rho = vars[SERVER_VAR_INDEX_RHO];
#endif

// Number density, if provided
#ifdef SERVER_VAR_INDEX_DEN
   spdata.n_dens = vars[SERVER_VAR_INDEX_DEN];
#endif

// Thermal pressure, if provided
#ifdef SERVER_VAR_INDEX_PTH
   spdata.p_ther = vars[SERVER_VAR_INDEX_PTH];
#endif

// Convert the variables to SPECTRUM format
   for (xyz = 0; xyz < 3; xyz++) {

// Bulk flow from mass density and momentum, if provided
#if defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO)
      spdata.Uvec[xyz] = vars[SERVER_VAR_INDEX_MOM + xyz] / rho;
// Bulk flow, if provided
#elif defined(SERVER_VAR_INDEX_FLO)
      spdata.Uvec[xyz] = vars[SERVER_VAR_INDEX_FLO + xyz];
#else
      spdata.Uvec[xyz] = 0.0;
#endif

// The magnetic field must be always provided
      spdata.Bvec[xyz] = vars[SERVER_VAR_INDEX_MAG + xyz];

// Electric field, if provided
#if defined(SERVER_VAR_INDEX_ELE)
      spdata.Evec[xyz] = vars[SERVER_VAR_INDEX_ELE + xyz];
#elif !defined(SERVER_VAR_INDEX_FLO) && !(defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
      spdata.Evec[xyz] = 0.0;
#endif

   };

// Electric field, if B and U provided
#ifndef SERVER_VAR_INDEX_ELE
#if defined(SERVER_VAR_INDEX_FLO) || (defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))   
   spdata.Evec = -(spdata.Uvec ^ spdata.Bvec) / c_code;
#endif
#endif

// Region(s) indicator variable(s), if provided
#ifdef SERVER_VAR_INDEX_REG
   for (xyz = 0; xyz < SERVER_NUM_INDEX_REG; xyz++) spdata.region[xyz] = vars[SERVER_VAR_INDEX_REG + xyz];
#else
   spdata.region = gv_zeros;
#endif
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
\param[in]  t      Time
\param[in]  pos    Position
\param[out] spdata Fields, dmax, etc.
*/
template <typename Fields>
void ServerCartesianFront::GetVariables(double t, const GeoVector& pos, Fields& fields)
{
   int bidx;

// Request the block and the zone size
   _inquiry.type = 1;
   _inquiry.pos = pos;
   bidx = RequestBlock();

// If "block_pri" or "block_sec" is the position owner (based on the call to RequestBlock), we don't need to acccess the cache
   if (block_pri->GetNode() != bidx) {
      if (block_sec->GetNode() == bidx) block_pri = block_sec;
      else block_pri = cache_line[bidx];
   };
   spdata.dmax = fmin(spdata.dmax, block_pri->GetZoneLength().Smallest());

#if SERVER_INTERP_ORDER == -1
// Get variables directly from reader program
   GetVariablesFromReader(spdata);
#elif SERVER_INTERP_ORDER == 0
// Get variables using 0th order interpolation
   GetVariablesInterp0(pos, spdata);
#elif SERVER_INTERP_ORDER == 1
// Get variables using 1st order interpolation
   GetVariablesInterp1(pos, spdata);
#else
#error Unsupported interpolation order!
#endif

// Perform unit conversion for fields and region
#ifdef SERVER_VAR_INDEX_DEN
   spdata.region /= spdata.n_dens;
#endif
   spdata.n_dens *= unit_number_density_server / unit_number_density_fluid;
   spdata.Uvec *= unit_velocity_server / unit_velocity_fluid;
   spdata.Bvec *= unit_magnetic_server / unit_magnetic_fluid;
   spdata.Evec *= unit_electric_server / unit_electric_fluid;
   spdata.p_ther *= unit_pressure_server / unit_pressure_fluid;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 03/11/2024
\param[out] spdata Fields, dmax, etc.
*/
template <typename Fields>
void ServerCartesianFront::GetGradientsInterp1(Fields& fields)
{
   double var, _Bmag2, rho = 0.0;
   int vidx, xyz, uvw, iz, pri_idx, sec_idx;

   double grads[n_variables][3] = {0.0};
   spdata.gradBmag = gv_zeros;

   LOWER_BITS(spdata._mask, BACKGROUND_grad_FAIL);

// Internal interpolation
   if (stencil_status == 0) {
      for (iz = 0; iz < stencil.n_elements; iz++) {
         for (vidx = 0; vidx < n_variables; vidx++) {
            var = block_pri->GetValue(stencil.zones[iz], vidx);
            for (uvw = 0; uvw < 3; uvw++) grads[vidx][uvw] += stencil.derivatives[3 * iz + uvw] * var;
         };
// Mass density, if provided
#ifdef SERVER_VAR_INDEX_RHO
         rho += stencil.weights[iz] * block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_RHO);
#endif
// gradient of B magnitude
         _Bmag2 = Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                + Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                + Sqr(block_pri->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         for (uvw = 0; uvw < 3; uvw++) spdata.gradBmag[uvw] += stencil.derivatives[3 * iz + uvw] * sqrt(_Bmag2);
      };
   }

// Edge or corner interpolation
   else if (stencil_status == 4) {
      pri_idx = block_pri->GetNode();
      sec_idx = block_sec->GetNode();
      for (iz = 0; iz < stencil.n_elements; iz++) {
         if (stencil.blocks[iz] == pri_idx) block_stn = block_pri;
         else if (stencil.blocks[iz] == sec_idx) block_stn = block_sec;
         else block_stn = cache_line[stencil.blocks[iz]];
         for (vidx = 0; vidx < n_variables; vidx++) {
            var = block_stn->GetValue(stencil.zones[iz], vidx);
            for (uvw = 0; uvw < 3; uvw++) grads[vidx][uvw] += stencil.derivatives[3 * iz + uvw] * var;
         };
// Mass density, if provided
#ifdef SERVER_VAR_INDEX_RHO
         rho += stencil.weights[iz] * block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_RHO);
#endif
// gradient of B magnitude
         _Bmag2 = Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         for (uvw = 0; uvw < 3; uvw++) spdata.gradBmag[uvw] += stencil.derivatives[3 * iz + uvw] * sqrt(_Bmag2);
      };
   }

// Plane interpolation
   else {
      pri_idx = block_pri->GetNode();
      for (iz = 0; iz < stencil.n_elements; iz++) {
         if (stencil.blocks[iz] == pri_idx) block_stn = block_pri;
         else block_stn = block_sec;
         for (vidx = 0; vidx < n_variables; vidx++) {
            var = block_stn->GetValue(stencil.zones[iz], vidx);
            for (uvw = 0; uvw < 3; uvw++) grads[vidx][uvw] += stencil.derivatives[3 * iz + uvw] * var;
         };
// Mass density, if provided
#ifdef SERVER_VAR_INDEX_RHO
         rho += stencil.weights[iz] * block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_RHO);
#endif
// gradient of B magnitude
         _Bmag2 = Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG))
                + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 1))
                + Sqr(block_stn->GetValue(stencil.zones[iz], SERVER_VAR_INDEX_MAG + 2));
         for (uvw = 0; uvw < 3; uvw++) spdata.gradBmag[uvw] += stencil.derivatives[3 * iz + uvw] * sqrt(_Bmag2);
      };
   };

// Convert the gradients to SPECTRUM format
   for (uvw = 0; uvw < 3; uvw++) {
      for (xyz = 0; xyz < 3; xyz++) {

// Bulk flow from mass density and momentum, if provided
#if defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO)
         spdata.gradUvec[uvw][xyz] = (grads[SERVER_VAR_INDEX_MOM + xyz][uvw] - spdata.Uvec[xyz] * grads[SERVER_VAR_INDEX_RHO][uvw]) / rho;
// Bulk flow, if provided
#elif defined(SERVER_VAR_INDEX_FLO)
         spdata.gradUvec[uvw][xyz] = grads[SERVER_VAR_INDEX_FLO + xyz][uvw];
#else
         spdata.gradUvec[uvw][xyz] = 0.0;
#endif

// The magnetic field must be always provided
         spdata.gradBvec[uvw][xyz] = grads[SERVER_VAR_INDEX_MAG + xyz][uvw];

// Electric field, if provided
#ifdef SERVER_VAR_INDEX_ELE
         spdata.gradEvec[uvw][xyz] = grads[SERVER_VAR_INDEX_ELE + xyz][uvw];
#elif !defined(SERVER_VAR_INDEX_FLO) && !(defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
         spdata.gradEvec[uvw][xyz] = 0.0;
#endif

      };
   };

// Electric field, if B and U provided
#ifndef SERVER_VAR_INDEX_ELE
#if defined(SERVER_VAR_INDEX_FLO) || (defined(SERVER_VAR_INDEX_MOM) && defined(SERVER_VAR_INDEX_RHO))
   spdata.gradEvec = -((spdata.gradUvec ^ spdata.Bvec) + (spdata.Uvec ^ spdata.gradBvec)) / c_code;
#endif
#endif

};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 07/19/2023
\param[out] spdata Field gradients
*/
template <typename Fields>
void ServerCartesianFront::GetGradients(Fields& fields)
{
#if SERVER_INTERP_ORDER == -1
// Gradients must be computed numerically
   RAISE_BITS(spdata._mask, BACKGROUND_grad_FAIL);
   return;
#elif SERVER_INTERP_ORDER == 0
// All gradients are explicitly set to zero, and the background must not attempt to compute them using "NumericalDerivatives()"
   spdata.gradUvec = gm_zeros;
   spdata.gradBvec = gm_zeros;
   spdata.gradBmag = gv_zeros;
   spdata.gradEvec = gm_zeros;
#elif SERVER_INTERP_ORDER == 1
// Gradients can be obtained from the stencil construction
   GetGradientsInterp1(spdata);
#else
#error Unsupported interpolation order!
#endif

// Perform unit conversion for gradients
   spdata.gradUvec *= unit_velocity_server / unit_velocity_fluid;
   spdata.gradBvec *= unit_magnetic_server / unit_magnetic_fluid;
   spdata.gradEvec *= unit_electric_server / unit_electric_fluid;
};


};
