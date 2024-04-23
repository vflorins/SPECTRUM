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

#include "fluid/conservation_laws_gasdyn.hh"
#include "fluid/conservation_laws_mhd.hh"
#ifdef GEO_DEBUG
#include "common/print_warn.hh"
#endif

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
template<typename cl_prim, typename cl_cons, typename cl_flux>
class RiemannSolverBase
{
protected:

//! Left primitive state
   cl_prim left_prim;

//! Right primitive state
   cl_cons rght_prim;

//! Resolved primitive state
   cl_prim resv_prim;

//! Left conserved state
   cl_cons left_cons;

//! Right conserved state
   cl_cons rght_cons;

//! Resolved conserved state
   cl_cons resv_cons;

//! Left flux function
   cl_flux left_flux;

//! Right flux function
   cl_flux rght_flux;

//! Resolved flux function
   cl_flux resv_flux;

//! Internal solver for the Riemann problem
   SPECTRUM_DEVICE_FUNC virtual void SolveInternal(void) {};

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC RiemannSolverBase(void) = default;

//! Solve the Riemann problem using primitive reconstruction
   SPECTRUM_DEVICE_FUNC void Solve(const cl_prim& left_prim_in, const cl_prim& rght_prim_in);

//! Solve the Riemann problem using conserved reconstruction
   SPECTRUM_DEVICE_FUNC void Solve(const cl_cons& left_cons_in, const cl_cons& rght_cons_in);

//! Return the resolved primitive state
   SPECTRUM_DEVICE_FUNC cl_prim GetResolvedPrim(void) const;

//! Return the resolved conserved state
   SPECTRUM_DEVICE_FUNC cl_cons GetResolvedCons(void) const;

//! Return the resolved flux function
   SPECTRUM_DEVICE_FUNC cl_flux GetResolvedFlux(void) const;
};

/*!
\author Vladimir Florinski
\date 03/18/2024
\params[in] left_prim_in Left primitive state
\params[in] rght_prim_in Right primitive state
*/
template<typename cl_prim, typename cl_cons, typename cl_flux>
SPECTRUM_DEVICE_FUNC inline void RiemannSolverBase<cl_prim, cl_cons, cl_flux>::Solve(const cl_prim& left_prim_in, const cl_prim& rght_prim_in)
{
   left_prim = left_prim_in;
   rght_prim = rght_prim_in;
   left_cons = left_prim.ToConserved(false);
   rght_cons = rght_prim.ToConserved(false);
   left_flux = left_prim.ToFlux(false);
   rght_flux = rght_prim.ToFlux(false);
   SolveInternal();
};

/*!
\author Vladimir Florinski
\date 03/18/2024
\params[in] left_cons_in Left conserved state
\params[in] rght_cons_in Right conserved state
*/
template<typename cl_prim, typename cl_cons, typename cl_flux>
SPECTRUM_DEVICE_FUNC inline void RiemannSolverBase<cl_prim, cl_cons, cl_flux>::Solve(const cl_cons& left_cons_in, const cl_cons& rght_cons_in)
{
   left_cons = left_cons_in;
   rght_cons = rght_cons_in;
   left_prim = left_cons.ToPrimitive(false);
   rght_prim = rght_cons.ToPrimitive(false);
   left_flux = left_prim.ToFlux(false);
   rght_flux = rght_prim.ToFlux(false);
   SolveInternal();
};

/*!
\author Vladimir Florinski
\date 03/07/2024
\return Resolved primitive state
*/
template<typename cl_prim, typename cl_cons, typename cl_flux>
SPECTRUM_DEVICE_FUNC inline cl_prim RiemannSolverBase<cl_prim, cl_cons, cl_flux>::GetResolvedPrim(void) const
{
   return resv_prim;
};

/*!
\author Vladimir Florinski
\date 03/18/2024
\return Resolved conserved state
*/
template<typename cl_prim, typename cl_cons, typename cl_flux>
SPECTRUM_DEVICE_FUNC inline cl_cons RiemannSolverBase<cl_prim, cl_cons, cl_flux>::GetResolvedCons(void) const
{
   return resv_cons;
};

/*!
\author Vladimir Florinski
\date 03/07/2024
\return Resolved flux function
*/
template<typename cl_prim, typename cl_cons, typename cl_flux>
SPECTRUM_DEVICE_FUNC inline cl_flux RiemannSolverBase<cl_prim, cl_cons, cl_flux>::GetResolvedFlux(void) const
{
   return resv_flux;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RiemannSolverRusanov class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Rusanov Riemann solver class template
\author Vladimir Florinski
*/
template<typename cl_prim, typename cl_cons, typename cl_flux>
class RiemannSolverRusanov : public RiemannSolverBase<cl_prim, cl_cons, cl_flux>
{
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::left_prim;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::rght_prim;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::resv_prim;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::left_cons;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::rght_cons;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::resv_cons;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::left_flux;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::rght_flux;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::resv_flux;

protected:

//! Speed of the extremal waves
   double S;

//! Internal solver for the Riemann problem
   SPECTRUM_DEVICE_FUNC void SolveInternal(void) override;

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC RiemannSolverRusanov(void) = default;
};

/*!
\author Vladimir Florinski
\date 03/18/2024
*/
template<typename cl_prim, typename cl_cons, typename cl_flux>
SPECTRUM_DEVICE_FUNC inline void RiemannSolverRusanov<cl_prim, cl_cons, cl_flux>::SolveInternal(void)
{
   double fastest_wave_left = left_prim.FastestWaveNormal();
   double fastest_wave_rght = rght_prim.FastestWaveNormal();
   S = std::max(std::abs(left_prim.vel[0]) + fastest_wave_left, std::abs(rght_prim.vel[0]) + fastest_wave_rght);

   resv_cons = 0.5 * (left_flux - rght_flux + S * (left_cons + rght_cons)) / S;
   resv_prim = resv_cons.ToPrimitive(false);
   resv_flux = 0.5 * (left_flux + rght_flux + S * (left_cons - rght_cons));
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RiemannSolverHLLE class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief HLLE Riemann solver class template
\author Vladimir Florinski
*/
template<typename cl_prim, typename cl_cons, typename cl_flux>
class RiemannSolverHLLE : public RiemannSolverBase<cl_prim, cl_cons, cl_flux>
{
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::left_prim;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::rght_prim;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::resv_prim;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::left_cons;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::rght_cons;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::resv_cons;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::left_flux;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::rght_flux;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::resv_flux;

protected:

//! Speed of the left extremal wave
   double SL;

//! Speed of the right extremal wave
   double SR;

//! Left flux in the moving frame
   cl_flux prime_flux_left;

//! Right flux in the moving frame
   cl_flux prime_flux_rght;

//! Iterate to find the extremal speeds SL and SR
   SPECTRUM_DEVICE_FUNC void IterateStarState(void);

//! Internal solver for the Riemann problem
   SPECTRUM_DEVICE_FUNC void SolveInternal(void) override;

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC RiemannSolverHLLE(void) = default;
};

/*!
\author Vladimir Florinski
\date 03/20/2024
*/
template<typename cl_prim, typename cl_cons, typename cl_flux>
SPECTRUM_DEVICE_FUNC inline void RiemannSolverHLLE<cl_prim, cl_cons, cl_flux>::IterateStarState(void)
{
#ifdef GEO_DEBUG
   int iter_count = 0;
#endif
   double SLStar, SRStar, fastest_wave_star;

   do {
      resv_cons = (SR * rght_cons - SL * left_cons - rght_flux + left_flux) / (SR - SL);
      resv_prim = resv_cons.ToPrimitive(false);
      fastest_wave_star = resv_prim.FastestWaveNormal();

      SLStar = resv_prim.vel[0] - fastest_wave_star;
      SRStar = resv_prim.vel[0] + fastest_wave_star;

      if((SLStar - SL >= 0.0) && (SR - SRStar >= 0.0)) break;
      SL = std::min(SL, SLStar);
      SR = std::max(SR, SRStar);

#ifdef GEO_DEBUG
      iter_count++;
#endif
   } while(true);
   
#ifdef GEO_DEBUG
   PrintMessage(__FILE__, __LINE__, std::to_string(iter_count) + " HLL iterations performed", true);
#endif
};

/*!
\author Vladimir Florinski
\date 03/18/2024
*/
template<typename cl_prim, typename cl_cons, typename cl_flux>
SPECTRUM_DEVICE_FUNC inline void RiemannSolverHLLE<cl_prim, cl_cons, cl_flux>::SolveInternal(void)
{
   double fastest_wave_left = left_prim.FastestWaveNormal();
   double fastest_wave_rght = rght_prim.FastestWaveNormal();

   SL = std::min(left_prim.vel()[0] - fastest_wave_left, rght_prim.vel()[0] - fastest_wave_rght);
   SR = std::max(left_prim.vel()[0] + fastest_wave_left, rght_prim.vel()[0] + fastest_wave_rght);

#if ITERATE_HLL_WAVE != 0
// Move the extremal waves as far as possible
   IterateStarState();
#endif

// Supersonic flow to the left, use the right state
   if(SR <= 0.0) {
      resv_cons = rght_cons;
      resv_flux = rght_prim.ToFlux(false);
      resv_prim = rght_prim;
   }

// Supersonic flow to the right, use the left state
   else if(SL >= 0.0) {
      resv_cons = left_cons;
      resv_flux = left_prim.ToFlux(false);
      resv_prim = left_prim;
   }

// Riemann fan is open - calculate the intermediate state and flux
   else {
      prime_flux_left = SL * left_cons - left_flux;
      prime_flux_rght = SR * rght_cons - rght_flux;
      resv_cons = (prime_flux_rght - prime_flux_left) / (SR - SL);
      resv_prim = resv_cons.ToPrimitive(false);
      resv_flux = (SR * left_flux - SL * rght_flux + SR * SL * (rght_cons - left_cons)) / (SR - SL);
   };   
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RiemannSolverHLLC class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief HLLC Riemann solver class template
\author Vladimir Florinski
*/
template<typename cl_prim, typename cl_cons, typename cl_flux>
class RiemannSolverHLLC : public RiemannSolverHLLE<cl_prim, cl_cons, cl_flux>
{
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::left_prim;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::rght_prim;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::resv_prim;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::left_cons;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::rght_cons;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::resv_cons;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::left_flux;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::rght_flux;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::resv_flux;

   using RiemannSolverHLLE<cl_prim, cl_cons, cl_flux>::SL;
   using RiemannSolverHLLE<cl_prim, cl_cons, cl_flux>::SR;
   using RiemannSolverHLLE<cl_prim, cl_cons, cl_flux>::prime_flux_left;
   using RiemannSolverHLLE<cl_prim, cl_cons, cl_flux>::prime_flux_rght;
   using RiemannSolverHLLE<cl_prim, cl_cons, cl_flux>::IterateStarState;

protected:

//! Speed of the centeral wave
   double SC;

//! Internal solver for the Riemann problem
   SPECTRUM_DEVICE_FUNC void SolveInternal(void) override;

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC RiemannSolverHLLC(void) = default;
};

/*!
\author Vladimir Florinski
\date 03/20/2024

Ref: Li, S., An HLLC Riemann Solver for Magneto-Hydrodynamics, Journal of Computational Physics, v. 203, p. 344 (2005).
*/
template<typename cl_prim, typename cl_cons, typename cl_flux>
SPECTRUM_DEVICE_FUNC inline void RiemannSolverHLLC<cl_prim, cl_cons, cl_flux>::SolveInternal(void)
{
   double resv_total_pre;

// Call the parent version
   RiemannSolverHLLE<cl_prim, cl_cons, cl_flux>::SolveInternal();

// Keep the HLLE result for supersonic flows
   if(SL * SR >= 0.0) return;

// For transsonic flows some work was already done in the parent function, namely V1,V5,V6,V7,U5,U6,U7
   SC = resv_prim.vel()[0];

// Intermediate speed to the right - compute the left intermediate state
   if(SC > 0.0) {
      resv_prim.den() = prime_flux_left.den() / (SL - SC); // V0
      resv_cons.mom()[0] = resv_prim.den() * SC; // U1
      resv_cons.FixEnergy(prime_flux_left, SL, SC); // U4
      resv_cons.den() = resv_prim.den(); // U0
      resv_cons.FixMomentum(prime_flux_left, SL, SC); // U2,U3
      resv_prim.vel()[1] = resv_cons.mom()[1] / resv_prim.den(); // V2
      resv_prim.vel()[2] = resv_cons.mom()[2] / resv_prim.den(); // V3
      resv_prim.pre() = resv_cons.GetPressure(); // V4
      resv_flux = left_flux - SL * (left_cons - resv_cons);
   }

// Intermediate speed to the left - compute the right intermediate state
   else {
      resv_prim.den() = prime_flux_rght.den() / (SR - SC); // V0
      resv_cons.mom()[0] = resv_prim.den() * SC; // U1
      resv_cons.FixEnergy(prime_flux_rght, SR, SC); // U4
      resv_cons.den() = resv_prim.den(); // U0
      resv_cons.FixMomentum(prime_flux_rght, SR, SC); // U2,U3
      resv_prim.vel()[1] = resv_cons.mom()[1] / resv_prim.den(); // V2
      resv_prim.vel()[2] = resv_cons.mom()[2] / resv_prim.den(); // V3
      resv_prim.pre() = resv_cons.GetPressure(); // V4
      resv_flux = rght_flux - SR * (rght_cons - resv_cons);
   };   

// Test for small pressure and possibly revert to HLLE
   if(resv_prim.pre() / resv_cons.enr() < tiny) {
      RiemannSolverHLLE<cl_prim, cl_cons, cl_flux>::SolveInternal();
#ifdef GEO_DEBUG
      PrintMessage(__FILE__, __LINE__, "HLLC solver pressure underflow, reverting to HLLE", true);
#endif
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RiemannSolverNoreflect class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Non-reflecting Riemann solver class template
\author Vladimir Florinski

Ref: Pogorelov, N. V., and Semenov, A. Y., Solar Wind Interaction with the Magnetized Interstellar Medium, Astronomy and Astrophysics, v. 321, p. 330 (1997).
*/
template<typename cl_prim, typename cl_cons, typename cl_flux>
class RiemannSolverNoreflect : public RiemannSolverBase<cl_prim, cl_cons, cl_flux>
{
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::left_prim;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::rght_prim;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::resv_prim;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::left_cons;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::rght_cons;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::resv_cons;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::left_flux;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::rght_flux;
   using RiemannSolverBase<cl_prim, cl_cons, cl_flux>::resv_flux;

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

/*!
\author Vladimir Florinski
\date 04/09/2024
*/
template<typename cl_prim, typename cl_cons, typename cl_flux>
SPECTRUM_DEVICE_FUNC inline void RiemannSolverNoreflect<cl_prim, cl_cons, cl_flux>::SolveInternal(void)
{
   double fastest_wave_left = left_prim.FastestWaveNormal();
   double fastest_wave_rght = rght_prim.FastestWaveNormal();

// Revert to the backup RS (HLLE) for inflow, super-luminal outflow, or sub-luminal BC
   if(dir > 0.0) {
      if((left_prim.vel()[0] < 0.0) || (fastest_wave_left - left_prim.vel()[0] < 0.0) || (fastest_wave_rght - rght_prim.vel()[0] > 0.0)) {
         RiemannSolverHLLE<cl_prim, cl_cons, cl_flux>::SolveInternal();
         return;
      };
   }
   else {
      if((rght_prim.vel()[0] > 0.0) || (fastest_wave_rght + rght_prim.vel()[0] < 0.0) || (fastest_wave_left + left_prim.vel()[0] > 0.0)) {
         RiemannSolverHLLE<cl_prim, cl_cons, cl_flux>::SolveInternal();
         return;
      };
   };

// Obtain the sonic state
   if(dir > 0.0) resv_prim = left_prim.ToSonic(false);
   else {
      rght_prim.Invert();
      resv_prim = rght_prim.ToSonic(false);
      rght_prim.Invert();
      resv_prim.Invert();
   };

   resv_cons = resv_prim.ToConserved(false);
   resv_flux = resv_prim.ToFlux(false);
};

};

#endif
