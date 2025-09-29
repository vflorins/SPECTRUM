/*!
\file riemann_solver.hh
\brief Declares some common one-dimensional Riemann solvers
\author Vladimir Florinski
\author XiaoCheng Guo
\author Dinshaw Balsara
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_RIEMANN_SOLVER_HH
#define SPECTRUM_RIEMANN_SOLVER_HH

#include "fluid_fields.hh"
#include "riemann_solver_config.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RiemannSolverBase class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Riemann solver base class template
\author Vladimir Florinski
\author Lucius Schoenbaum
*/
template <typename RiemannSolverConfig_>
class RiemannSolverBase {
public:

   using RiemannSolverConfig = RiemannSolverConfig_;
   static constexpr Model model = RiemannSolverConfig::model;
   static constexpr Specie specie = RiemannSolverConfig::specie;
   static constexpr Passivity passivity = RiemannSolverConfig::passivity;

   using PrimitiveState = FluidSpecie<model, Form::primitive, specie.id, passivity>;
   using ConservedState = FluidSpecie<model, Form::conserved, specie.id, passivity>;
   using FluxFunction = FluidSpecie<model, Form::flux, specie.id, passivity>;

protected:

//! Left primitive state
   PrimitiveState left_prim;

//! Right primitive state
   PrimitiveState rght_prim;

//! Resolved primitive state
   PrimitiveState resv_prim;

//! Left conserved state
   ConservedState left_cons;

//! Right conserved state
   ConservedState rght_cons;

//! Resolved conserved state
   ConservedState resv_cons;

//! Left flux function
   FluxFunction left_flux;

//! Right flux function
   FluxFunction rght_flux;

//! Resolved flux function
   FluxFunction resv_flux;

//! Internal solver for the Riemann problem
   SPECTRUM_DEVICE_FUNC void SolveInternal(void) {};

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC RiemannSolverBase(void) = default;

//! Solve the Riemann problem using primitive reconstruction
   SPECTRUM_DEVICE_FUNC void Solve(const PrimitiveState& left_prim_in, const PrimitiveState& rght_prim_in);

//! Solve the Riemann problem using conserved reconstruction
   SPECTRUM_DEVICE_FUNC void Solve(const ConservedState& left_cons_in, const ConservedState& rght_cons_in);

//! Return the resolved state
   template <Form form>
   SPECTRUM_DEVICE_FUNC decltype(auto) GetResolved(void) const {
      if constexpr (form == Form::primitive) {
         return resv_prim;
      }
      else if constexpr (form == Form::conserved) {
         return resv_cons;
      }
      else {
         return resv_flux;
      }
   }

};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RiemannSolverRusanov class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Rusanov Riemann solver class template
\author Vladimir Florinski
*/
template <typename RiemannSolverConfig_>
class RiemannSolverRusanov : public RiemannSolverBase<RiemannSolverConfig_>
{
   using RiemannSolverConfig = RiemannSolverConfig_;
   using RiemannSolverBase = RiemannSolverBase<RiemannSolverConfig>;
   using RiemannSolverBase::model;
   using RiemannSolverBase::specie;
   using RiemannSolverBase::passivity;
   using RiemannSolverBase::left_prim;
   using RiemannSolverBase::rght_prim;
   using RiemannSolverBase::resv_prim;
   using RiemannSolverBase::left_cons;
   using RiemannSolverBase::rght_cons;
   using RiemannSolverBase::resv_cons;
   using RiemannSolverBase::left_flux;
   using RiemannSolverBase::rght_flux;
   using RiemannSolverBase::resv_flux;

protected:

//! Speed of the extremal waves
   double S;

//! Internal solver for the Riemann problem
   SPECTRUM_DEVICE_FUNC void SolveInternal(void);

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC RiemannSolverRusanov(void) = default;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RiemannSolverHLLE class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief HLLE Riemann solver class template
\author Vladimir Florinski
*/
template <typename RiemannSolverConfig_>
class RiemannSolverHLLE : public RiemannSolverBase<RiemannSolverConfig_>
{
   using RiemannSolverConfig = RiemannSolverConfig_;
   using RiemannSolverBase = RiemannSolverBase<RiemannSolverConfig>;
   using RiemannSolverBase::model;
   using RiemannSolverBase::specie;
   using RiemannSolverBase::passivity;
   using RiemannSolverBase::left_prim;
   using RiemannSolverBase::rght_prim;
   using RiemannSolverBase::resv_prim;
   using RiemannSolverBase::left_cons;
   using RiemannSolverBase::rght_cons;
   using RiemannSolverBase::resv_cons;
   using RiemannSolverBase::left_flux;
   using RiemannSolverBase::rght_flux;
   using RiemannSolverBase::resv_flux;

   using FluxFunction = FluidSpecie<model, Form::flux, specie.id, passivity>;

protected:

//! Speed of the left extremal wave
   double SL;

//! Speed of the right extremal wave
   double SR;

//! Left flux in the moving frame
   FluxFunction prime_flux_left;

//! Right flux in the moving frame
   FluxFunction prime_flux_rght;

//! Iterate to find the extremal speeds SL and SR
   SPECTRUM_DEVICE_FUNC void IterateStarState(void);

//! Internal solver for the Riemann problem
   SPECTRUM_DEVICE_FUNC void SolveInternal(void);

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC RiemannSolverHLLE(void) = default;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RiemannSolverHLLC class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief HLLC Riemann solver class template
\author Vladimir Florinski
*/
template <typename RiemannSolverConfig_>
class RiemannSolverHLLC : public RiemannSolverHLLE<RiemannSolverConfig_>
{
   using RiemannSolverConfig = RiemannSolverConfig_;
   using RiemannSolverBase = RiemannSolverBase<RiemannSolverConfig>;
   using RiemannSolverBase::model;
   using RiemannSolverBase::specie;
   using RiemannSolverBase::passivity;
   using RiemannSolverBase::left_prim;
   using RiemannSolverBase::rght_prim;
   using RiemannSolverBase::resv_prim;
   using RiemannSolverBase::left_cons;
   using RiemannSolverBase::rght_cons;
   using RiemannSolverBase::resv_cons;
   using RiemannSolverBase::left_flux;
   using RiemannSolverBase::rght_flux;
   using RiemannSolverBase::resv_flux;

   using RiemannSolverHLLE = RiemannSolverHLLE<RiemannSolverConfig>;
   using RiemannSolverHLLE::SL;
   using RiemannSolverHLLE::SR;
   using RiemannSolverHLLE::prime_flux_left;
   using RiemannSolverHLLE::prime_flux_rght;
   using RiemannSolverHLLE::IterateStarState;

protected:

//! Speed of the centeral wave
   double SC;

//! Internal solver for the Riemann problem
   SPECTRUM_DEVICE_FUNC void SolveInternal(void);

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC RiemannSolverHLLC(void) = default;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RiemannSolverNoreflect class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Non-reflecting Riemann solver class template
\author Vladimir Florinski

Ref: Pogorelov, N. V., and Semenov, A. Y., Solar Wind Interaction with the Magnetized Interstellar Medium, Astronomy and Astrophysics, v. 321, p. 330 (1997).
*/
template <typename RiemannSolverConfig_>
class RiemannSolverNoreflect : public RiemannSolverHLLE<RiemannSolverConfig_>
{
   using RiemannSolverConfig = RiemannSolverConfig_;
   using RiemannSolverBase = RiemannSolverBase<RiemannSolverConfig>;
   using RiemannSolverBase::model;
   using RiemannSolverBase::specie;
   using RiemannSolverBase::passivity;
   using RiemannSolverBase::left_prim;
   using RiemannSolverBase::rght_prim;
   using RiemannSolverBase::resv_prim;
   using RiemannSolverBase::left_cons;
   using RiemannSolverBase::rght_cons;
   using RiemannSolverBase::resv_cons;
   using RiemannSolverBase::left_flux;
   using RiemannSolverBase::rght_flux;
   using RiemannSolverBase::resv_flux;

   using RiemannSolverHLLE = RiemannSolverHLLE<RiemannSolverConfig>;

protected:

//! Direction of the outflow (+1 or -1)
   int dir;

//! Swith the left and right states and the vector directions
//   SPECTRUM_DEVICE_FUNC void FlipStates(void);

//! Internal solver for the Riemann problem
   SPECTRUM_DEVICE_FUNC void SolveInternal(void);

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC RiemannSolverNoreflect(void) = default;
};

};

#endif
