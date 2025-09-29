/*!
\file riemann_solver.hh
\brief Declares some common one-dimensional Riemann solvers
\author Vladimir Florinski
\author XiaoCheng Guo
\author Dinshaw Balsara

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_RIEMANN_SOLVER_HH
#define SPECTRUM_RIEMANN_SOLVER_HH

#include "fluid/impl1/conservation_laws_gasdyn.hh"
#include "fluid/impl1/conservation_laws_mhd.hh"

namespace Spectrum {

//! Enables iteration to calculate maximum extent of the Riemann fan
#define ITERATE_HLL_WAVE false

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RiemannSolverBase class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Riemann solver base class template
\author Vladimir Florinski
*/
template<typename ConsLaw>
class RiemannSolverBase
{
protected:

//! Left primitive state
   ConsLaw::PrimitiveState left_prim;

//! Right primitive state
   ConsLaw::PrimitiveState rght_prim;

//! Resolved primitive state
   ConsLaw::PrimitiveState resv_prim;

//! Left conserved state
   ConsLaw::ConservedState left_cons;

//! Right conserved state
   ConsLaw::ConservedState rght_cons;

//! Resolved conserved state
   ConsLaw::ConservedState resv_cons;

//! Left flux function
   ConsLaw::FluxFunction left_flux;

//! Right flux function
   ConsLaw::FluxFunction rght_flux;

//! Resolved flux function
   ConsLaw::FluxFunction resv_flux;

//! Internal solver for the Riemann problem
   SPECTRUM_DEVICE_FUNC virtual void SolveInternal(void) {};

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC RiemannSolverBase(void) = default;

//! Solve the Riemann problem using primitive reconstruction
   SPECTRUM_DEVICE_FUNC void Solve(const ConsLaw::PrimitiveState& left_prim_in, const ConsLaw::PrimitiveState& rght_prim_in);

//! Solve the Riemann problem using conserved reconstruction
   SPECTRUM_DEVICE_FUNC void Solve(const ConsLaw::ConservedState& left_cons_in, const ConsLaw::ConservedState& rght_cons_in);

//! Return the resolved primitive state
   SPECTRUM_DEVICE_FUNC ConsLaw::PrimitiveState GetResolvedPrim(void) const;

//! Return the resolved conserved state
   SPECTRUM_DEVICE_FUNC ConsLaw::ConservedState GetResolvedCons(void) const;

//! Return the resolved flux function
   SPECTRUM_DEVICE_FUNC ConsLaw::FluxFunction GetResolvedFlux(void) const;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RiemannSolverRusanov class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Rusanov Riemann solver class template
\author Vladimir Florinski
*/
template<typename ConsLaw>
class RiemannSolverRusanov : public RiemannSolverBase<ConsLaw>
{
   using RiemannSolverBase<ConsLaw>::left_prim;
   using RiemannSolverBase<ConsLaw>::rght_prim;
   using RiemannSolverBase<ConsLaw>::resv_prim;
   using RiemannSolverBase<ConsLaw>::left_cons;
   using RiemannSolverBase<ConsLaw>::rght_cons;
   using RiemannSolverBase<ConsLaw>::resv_cons;
   using RiemannSolverBase<ConsLaw>::left_flux;
   using RiemannSolverBase<ConsLaw>::rght_flux;
   using RiemannSolverBase<ConsLaw>::resv_flux;

protected:

//! Speed of the extremal waves
   double S;

//! Internal solver for the Riemann problem
   SPECTRUM_DEVICE_FUNC void SolveInternal(void) override;

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
template<typename ConsLaw>
class RiemannSolverHLLE : public RiemannSolverBase<ConsLaw>
{
   using RiemannSolverBase<ConsLaw>::left_prim;
   using RiemannSolverBase<ConsLaw>::rght_prim;
   using RiemannSolverBase<ConsLaw>::resv_prim;
   using RiemannSolverBase<ConsLaw>::left_cons;
   using RiemannSolverBase<ConsLaw>::rght_cons;
   using RiemannSolverBase<ConsLaw>::resv_cons;
   using RiemannSolverBase<ConsLaw>::left_flux;
   using RiemannSolverBase<ConsLaw>::rght_flux;
   using RiemannSolverBase<ConsLaw>::resv_flux;

protected:

//! Speed of the left extremal wave
   double SL;

//! Speed of the right extremal wave
   double SR;

//! Left flux in the moving frame
   ConsLaw::FluxFunction prime_flux_left;

//! Right flux in the moving frame
   ConsLaw::FluxFunction prime_flux_rght;

//! Iterate to find the extremal speeds SL and SR
   SPECTRUM_DEVICE_FUNC void IterateStarState(void);

//! Internal solver for the Riemann problem
   SPECTRUM_DEVICE_FUNC void SolveInternal(void) override;

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
template<typename ConsLaw>
class RiemannSolverHLLC : public RiemannSolverHLLE<ConsLaw>
{
   using RiemannSolverBase<ConsLaw>::left_prim;
   using RiemannSolverBase<ConsLaw>::rght_prim;
   using RiemannSolverBase<ConsLaw>::resv_prim;
   using RiemannSolverBase<ConsLaw>::left_cons;
   using RiemannSolverBase<ConsLaw>::rght_cons;
   using RiemannSolverBase<ConsLaw>::resv_cons;
   using RiemannSolverBase<ConsLaw>::left_flux;
   using RiemannSolverBase<ConsLaw>::rght_flux;
   using RiemannSolverBase<ConsLaw>::resv_flux;

   using RiemannSolverHLLE<ConsLaw>::SL;
   using RiemannSolverHLLE<ConsLaw>::SR;
   using RiemannSolverHLLE<ConsLaw>::prime_flux_left;
   using RiemannSolverHLLE<ConsLaw>::prime_flux_rght;
   using RiemannSolverHLLE<ConsLaw>::IterateStarState;

protected:

//! Speed of the centeral wave
   double SC;

//! Internal solver for the Riemann problem
   SPECTRUM_DEVICE_FUNC void SolveInternal(void) override;

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
template<typename ConsLaw>
class RiemannSolverNoreflect : public RiemannSolverBase<ConsLaw>
{
   using RiemannSolverBase<ConsLaw>::left_prim;
   using RiemannSolverBase<ConsLaw>::rght_prim;
   using RiemannSolverBase<ConsLaw>::resv_prim;
   using RiemannSolverBase<ConsLaw>::left_cons;
   using RiemannSolverBase<ConsLaw>::rght_cons;
   using RiemannSolverBase<ConsLaw>::resv_cons;
   using RiemannSolverBase<ConsLaw>::left_flux;
   using RiemannSolverBase<ConsLaw>::rght_flux;
   using RiemannSolverBase<ConsLaw>::resv_flux;

protected:

//! Direction of the outflow (+1 or -1)
   int dir;

//! Swith the left and right states and the vector directions
//   SPECTRUM_DEVICE_FUNC void FlipStates(void);

//! Internal solver for the Riemann problem
   SPECTRUM_DEVICE_FUNC void SolveInternal(void) override;

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC RiemannSolverNoreflect(void) = default;
};

};

#endif
