/*!
\file variables.hh
\author Vladimir Florinski
\author Lucius Schoenbaum
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_VARIABLESSPECIES_HH
#define SPECTRUM_VARIABLESSPECIES_HH

#include "common/physics.hh"
#include "partially_generated/fields.hh"



// todo draft - deprecated



namespace Spectrum {

//const constexpr std::array<std::string_view, 12> SpeciesNames = {
//      std::string_view("MHD"),
//      std::string_view("GasDyn"),
//      std::string_view("Pickup-Ion"),
//};


/*!
\brief Cell-wise variables for a fluid or gas distribution
\author Lucius Schoenbaum
\author Vladimir Florinski
\date 03/25/2025
Some wrapper code has been added to make the very generic std::variant
interface more convenient for a physical/numerical application.
From the outside, it is as though all variables (primitive, and conservative)
are stored, but only one family is stored at a time. In some situations,
the user may wish to guarantee that one version is the one that is stored
(in a computation-heavy block that uses only one family primarily).
In that case, the user can call `ensureConserved` or `ensurePrimitive`.
The flux function for all conserved variables is obtained by calling
`getFluxFunction`.
 */
//template <unsigned int fluid_>
//struct VariablesGASDYN: ConservedStateGASDYN_t {
//
//   static constexpr const unsigned int fluid = fluid_;
////   static constexpr const std::string_view name = SpeciesNames[fluid];
//
///*!
//\brief velocity
//\author Lucius Schoenbaum
//\author Vladimir Florinski
//*/
//   Vel_t Vel() const {
//      return Mom() / Den();
//   }
//
//   void Vel_assign(Vel_t vel) {
//      Mom() = vel*Den();
//   }
//
///*!
//\brief pressure
//\author Lucius Schoenbaum
//\author Vladimir Florinski
//*/
//   Prs_t Prs() const {
//      return Pressure(Den(), Vel().Norm2(), 0.0, Enr(), fluid);
//   }
//
//   void Prs_assign(Prs_t prs) {
//      Enr() = Energy(Den(), Vel().Norm2(), 0.0, prs, fluid);;
//   }
//
//
//public: // ctor:
//
//   /*
//    * Note: copy, move, and default constructors are provided by std::variant.
//    *
//    */
//
//
//public: // getters:
//
//   ConservedStateGASDYN_t getConservedState() {
//      return static_cast<ConservedStateGASDYN_t>(*this);
//   }
//
//   PrimitiveStateGASDYN_t getPrimitiveState() {
//      return PrimitiveStateGASDYNT(Den(),Vel(),Prs());
//   }
//
//   FluxFunctionGASDYN_t getFluxFunction() {
//      return FluxFunctionGASDYN_t();
//   }
//
//
//public: // derived values:
//
///*!
// * Fastest wave speed (sound)
// *
// */
//   double FastestWave() const {
//      return SoundSpeed(Den(), Prs(), fluid);
//   }
//
//
///*!
// * Fastest normal wave speed (sound)
// *
// */
//   double FastestWaveNormal() const {
//      //return SoundSpeed(den(), prs(), fluid);
//      return 0.0;
//   }
//
///*!
// * Vector of primitive variables in a sonic state
// *
// */
//   VariablesGASDYN toSonic() const {
//      auto variables = VariablesGASDYN();
//      return variables;
//   }
//
//   void Invert() {
//      Mom() = -Mom();
//   }
//
//
//public: // helpers:
//
//   ;
//
//};


};


#endif
