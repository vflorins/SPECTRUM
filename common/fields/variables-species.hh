/*!
\file variables.hh
\author Vladimir Florinski
\author Lucius Schoenbaum
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_FIELDS_HH
#define SPECTRUM_FIELDS_HH

#include "common/physics.hh"
#include "partially_generated/fields.hh"



// todo this used to be variables but probably should become species but I'm not sure yet with all the changes happening



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
template <unsigned int fluid_>
struct VariablesGASDYN: ConservedStateGASDYN_t {

   static constexpr const unsigned int fluid = fluid_;
   // todo __________________________________________________________________
//   static constexpr const std::string_view name = SpeciesNames[fluid];

/*!
\brief velocity
\author Lucius Schoenbaum
\author Vladimir Florinski
*/
   vel_t vel() const {
      return mom() / den();
   }

   void vel_assign(vel_t vel) {
      mom() = vel*den();
   }

/*!
\brief pressure
\author Lucius Schoenbaum
\author Vladimir Florinski
*/
   prs_t prs() const {
      return Pressure(den(), vel().Norm2(), 0.0, enr(), fluid);
   }

   void prs_assign(prs_t prs) {
      enr() = Energy(den(), vel().Norm2(), 0.0, prs, fluid);;
   }


public: // ctor:

   /*
    * Note: copy, move, and default constructors are provided by std::variant.
    *
    */


public: // getters:

   ConservedStateGASDYN_t getConservedState() {
      return static_cast<ConservedStateGASDYN_t>(*this);
   }

   PrimitiveStateGASDYN_t getPrimitiveState() {
      return PrimitiveStateGASDYNT(den(),vel(),prs());
   }

   FluxFunctionGASDYN_t getFluxFunction() {
      // todo __________________________________________________________________
      return FluxFunctionGASDYN_t();
   }


public: // derived values:

/*!
 * Fastest wave speed (sound)
 *
 */
   double FastestWave() const {
      return SoundSpeed(den(), prs(), fluid);
   }


/*!
 * Fastest normal wave speed (sound)
 *
 */
   double FastestWaveNormal() const {
      // todo __________________________________________________________________
      //return SoundSpeed(den(), prs(), fluid);
      return 0.0;
   }

/*!
 * Vector of primitive variables in a sonic state
 *
 */
   VariablesGASDYN toSonic() const {
      // todo __________________________________________________________________
      auto variables = VariablesGASDYN();
      return variables;
   }

   void Invert() {
      mom() = -mom();
   }


public: // helpers:

   ;

};


};


#endif
