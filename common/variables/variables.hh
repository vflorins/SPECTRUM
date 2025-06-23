/*!
\file variables.hh
\author Vladimir Florinski
\author Lucius Schoenbaum
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_VARIABLES_HH
#define SPECTRUM_VARIABLES_HH

#include <variant>
#include "mhdtuple/mhdtuple.hh"

#include <iostream>
using std::cout; using std::endl;


namespace Spectrum {

// QUASI-REALISTIC STUB - set up all fluid species that you wish to be handled by the solver

#define SPECIES_MHD 0
#define SPECIES_GASDYN 1
#define SPECIES_PICKUP_ION 2

const constexpr std::array<std::string_view, 12> SpeciesNames = {
      std::string_view("MHD"),
      std::string_view("GasDyn"),
      std::string_view("Pickup-Ion"),
};


// COPY - FOR THE SAKE OF THE TEST VERSION - don't include physics.hh or link gsl

#define M_8PI     2.513274122871834590770114706623602307E+1

//! Polytropic indices
// Indexed by fluid type number
constexpr double gamma_eos[] = {4.0 / 3.0, 5.0 / 3.0, 6.0 / 3.0, 7.0 / 3.0};

/*!
\brief Calculate the MHD energy density
\author Vladimir Florinski
\date 07/31/2019
\param[in] den Mass density
\param[in] u2  Square of velocity
\param[in] B2  Square of the magnetic field
\param[in] pre Gas pressure
\param[in] ifl Index of the fluid
\return Energy per unit volume
*/
inline double Energy(double den, double u2, double B2, double pre, unsigned int ifl = 0)
{
   return den * u2 / 2.0 + pre / (gamma_eos[ifl] - 1.0) + B2 / M_8PI;
};

/*!
\brief Calculate gas pressure from MHD energy density
\author Vladimir Florinski
\date 07/31/2019
\param[in] den Mass density
\param[in] u2  Square of velocity
\param[in] B2  Square of the magnetic field
\param[in] enr MHD energy density
\param[in] ifl Index of the fluid
\return gas pressure
*/
SPECTRUM_DEVICE_FUNC inline double Pressure(double den, double u2, double B2, double enr, unsigned int ifl = 0)
{
   return (enr - den * u2 / 2.0 - B2 / M_8PI) * (gamma_eos[ifl] - 1.0);
};



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
   static constexpr const std::string_view name = SpeciesNames[fluid];

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
