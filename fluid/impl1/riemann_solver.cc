/*!
\file riemann_solver.cc
\brief Implements several common one-dimensional Riemann solvers
\author Vladimir Florinski
\author XiaoCheng Guo
\author Dinshaw Balsara

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "fluid/impl1/riemann_solver.hh"
#ifdef GEO_DEBUG
#include "common/print_warn.hh"
#endif

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RiemannSolverBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/15/2025
\params[in] left_prim_in Left primitive state
\params[in] rght_prim_in Right primitive state
*/
template<typename ConsLaw>
SPECTRUM_DEVICE_FUNC void RiemannSolverBase<ConsLaw>::Solve(const ConsLaw::PrimitiveState& left_prim_in, const ConsLaw::PrimitiveState& rght_prim_in)
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
\date 01/15/2025
\params[in] left_cons_in Left conserved state
\params[in] rght_cons_in Right conserved state
*/
template<typename ConsLaw>
SPECTRUM_DEVICE_FUNC void RiemannSolverBase<ConsLaw>::Solve(const ConsLaw::ConservedState& left_cons_in, const ConsLaw::ConservedState& rght_cons_in)
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
\date 01/15/2025
\return Resolved primitive state
*/
template<typename ConsLaw>
SPECTRUM_DEVICE_FUNC ConsLaw::PrimitiveState RiemannSolverBase<ConsLaw>::GetResolvedPrim(void) const
{
   return resv_prim;
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Resolved conserved state
*/
template<typename ConsLaw>
SPECTRUM_DEVICE_FUNC ConsLaw::ConservedState RiemannSolverBase<ConsLaw>::GetResolvedCons(void) const
{
   return resv_cons;
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Resolved flux function
*/
template<typename ConsLaw>
SPECTRUM_DEVICE_FUNC ConsLaw::FluxFunction RiemannSolverBase<ConsLaw>::GetResolvedFlux(void) const
{
   return resv_flux;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RiemannSolverRusanov methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/15/2025
*/
template<typename ConsLaw>
SPECTRUM_DEVICE_FUNC void RiemannSolverRusanov<ConsLaw>::SolveInternal(void)
{
   double fastest_wave_left = left_prim.FastestWaveNormal();
   double fastest_wave_rght = rght_prim.FastestWaveNormal();
   S = std::max(std::abs(left_prim.vel[0]) + fastest_wave_left, std::abs(rght_prim.vel[0]) + fastest_wave_rght);

   resv_cons = 0.5 * (left_flux - rght_flux + S * (left_cons + rght_cons)) / S;
   resv_prim = resv_cons.ToPrimitive(false);
   resv_flux = 0.5 * (left_flux + rght_flux + S * (left_cons - rght_cons));
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RiemannSolverHLLE methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/15/2025
*/
template<typename ConsLaw>
SPECTRUM_DEVICE_FUNC void RiemannSolverHLLE<ConsLaw>::IterateStarState(void)
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
\date 01/15/2025

Ref: Einfeldt, B., On Godunov-type methods for gas dynamics, SIAM Journal on Numerical Analysis, v. 25, p. 294 (1988).
*/
template<typename ConsLaw>
SPECTRUM_DEVICE_FUNC void RiemannSolverHLLE<ConsLaw>::SolveInternal(void)
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
// RiemannSolverHLLC methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/15/2025

Ref: Li, S., An HLLC Riemann Solver for Magneto-Hydrodynamics, Journal of Computational Physics, v. 203, p. 344 (2005).
*/
template<typename ConsLaw>
SPECTRUM_DEVICE_FUNC void RiemannSolverHLLC<ConsLaw>::SolveInternal(void)
{
   double resv_total_pre;

// Call the parent version
   RiemannSolverHLLE<ConsLaw>::SolveInternal();

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
   if(resv_prim.pre() / resv_cons.enr() < sp_tiny) {
      RiemannSolverHLLE<ConsLaw>::SolveInternal();
#ifdef GEO_DEBUG
      PrintMessage(__FILE__, __LINE__, "HLLC solver pressure underflow, reverting to HLLE", true);
#endif
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RiemannSolverNoreflect methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/15/2025
*/
template<typename ConsLaw>
SPECTRUM_DEVICE_FUNC void RiemannSolverNoreflect<ConsLaw>::SolveInternal(void)
{
   double fastest_wave_left = left_prim.FastestWaveNormal();
   double fastest_wave_rght = rght_prim.FastestWaveNormal();

// Revert to the backup RS (HLLE) for inflow, super-luminal outflow, or sub-luminal BC
   if(dir > 0.0) {
      if((left_prim.vel()[0] < 0.0) || (fastest_wave_left - left_prim.vel()[0] < 0.0) || (fastest_wave_rght - rght_prim.vel()[0] > 0.0)) {
         RiemannSolverHLLE<ConsLaw>::SolveInternal();
         return;
      };
   }
   else {
      if((rght_prim.vel()[0] > 0.0) || (fastest_wave_rght + rght_prim.vel()[0] < 0.0) || (fastest_wave_left + left_prim.vel()[0] > 0.0)) {
         RiemannSolverHLLE<ConsLaw>::SolveInternal();
         return;
      };
   };

// Obtain the sonic state. The standard situation is the interior cell on the left and the ghost cell on the right. For the opposite case all vector quantities must be inverted on input to "ToSonic()" and then inverted back.
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
