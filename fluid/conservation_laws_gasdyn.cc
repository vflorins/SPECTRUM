/*!
\file conservation_laws_gasdyn.cc
\brief Implements rules for compressible gas-dynamic equations
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "common/physics.hh"
#include "fluid/conservation_laws_gasdyn.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Gasdyn::PrimitiveState methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] den_in Density
\param[in] vel_in Velocity
\param[in] pre_in Pressure
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC Gasdyn<fluid>::PrimitiveState::PrimitiveState(double den_in, const GeoVector&vel_in, double pre_in)
{
   den() = den_in;
   vel() = vel_in;
   pre() = pre_in;
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] other Object to initialize from
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC Gasdyn<fluid>::PrimitiveState::PrimitiveState(const SimpleArray<double, CL_GASDYN_NVARS>& other)
                                                  : SimpleArray<double, CL_GASDYN_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Fastest wave speed (sound)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC double Gasdyn<fluid>::PrimitiveState::FastestWave(void) const
{
   return SoundSpeed(den(), pre(), fluid);
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Fastest normal wave speed (sound)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC double Gasdyn<fluid>::PrimitiveState::FastestWaveNormal(void) const
{
   return SoundSpeed(den(), pre(), fluid);
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Vector of conserved variables
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC Gasdyn<fluid>::ConservedState Gasdyn<fluid>::PrimitiveState::ToConserved(void) const
{
   Gasdyn<fluid>::ConservedState cons;

   cons.den() = den();
   cons.mom() = den() * vel();
   cons.enr() = Energy(den(), vel().Norm2(), 0.0, pre(), fluid);

   return cons;
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Flux vector
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC Gasdyn<fluid>::FluxFunction Gasdyn<fluid>::PrimitiveState::ToFlux(void) const
{
   double enr = Energy(den(), vel().Norm2(), 0.0, pre(), fluid);
   Gasdyn<fluid>::FluxFunction flux;

   flux.denf() = den() * vel()[0];
   flux.momf() = den() * vel()[0] * vel() + pre() * gv_nx;
   flux.enrf() = (enr + pre()) * vel()[0];

   return flux;
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Vector of primitive variables in a sonic state
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC Gasdyn<fluid>::PrimitiveState Gasdyn<fluid>::PrimitiveState::ToSonic(void) const
{
   double cs = Sound(den(), pre(), fluid);
   double del_den = 2.0 * den() / (gamma_eos[fluid] + 1.0) * (vel()[0] / cs - 1.0);
   double cs_star = vel()[0] - cs * del_den / den();

   Gasdyn<fluid>::PrimitiveState prim;
   
   prim.den() = den() + del_den;
   prim.vel()[0] = cs_star;
   prim.vel()[1] = vel()[1];
   prim.vel()[2] = vel()[2];
   prim.pre() = prim.den() * Sqr(cs_star) / gamma_eos[fluid];
   return prim;
};

/*!
\author Vladimir Florinski
\date 01/15/2025
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC void Gasdyn<fluid>::PrimitiveState::Invert(void)
{
   vel()[0] = -vel()[0];
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Gasdyn::ConservedState methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] den_in Density
\param[in] mom_in Momentum
\param[in] enr_in Energy
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC Gasdyn<fluid>::ConservedState::ConservedState(double den_in, const GeoVector& mom_in, double enr_in)
{
   den() = den_in;
   mom() = mom_in;
   enr() = enr_in;
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] other Object to initialize from
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC Gasdyn<fluid>::ConservedState::ConservedState(const SimpleArray<double, CL_GASDYN_NVARS>& other)
                                                  : SimpleArray<double, CL_GASDYN_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Fastest wave speed (sound)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC double Gasdyn<fluid>::ConservedState::FastestWave(void) const
{
   GeoVector vel = mom() / den();
   double pre = Pressure(den(), vel.Norm2(), 0.0, enr(), fluid);
   return SoundSpeed(den(), pre, fluid);
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Fastest normal wave speed (sound)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC double Gasdyn<fluid>::ConservedState::FastestWaveNormal(void) const
{
   GeoVector vel = mom() / den();
   double pre = Pressure(den(), vel.Norm2(), 0.0, enr(), fluid);
   return SoundSpeed(den(), pre, fluid);
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Vector of primitive variables
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC Gasdyn<fluid>::PrimitiveState Gasdyn<fluid>::ConservedState::ToPrimitive(void) const
{
   Gasdyn<fluid>::PrimitiveState prim;

   prim.den() = den();
   prim.vel() = mom() / den();
   prim.pre() = Pressure(den(), prim.vel().Norm2(), 0.0, enr(), fluid);

   return prim;
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Flux vector
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC Gasdyn<fluid>::FluxFunction Gasdyn<fluid>::ConservedState::ToFlux(void) const
{
   double pre = Pressure(den(), mom().Norm2() / Sqr(den()), 0.0, enr(), fluid);

   Gasdyn<fluid>::FluxFunction flux;
   flux.denf() = mom()[0];
   flux.momf() = mom()[0] * mom() / den() + pre * gv_nx;
   flux.enrf() = (enr() + pre) * mom()[0] / den();

   return flux;
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Gas pressure
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC double Gasdyn<fluid>::ConservedState::GetPressure(void) const
{
   return Pressure(den(), mom().Norm2() / Sqr(den()), 0.0, enr(), fluid);
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] flux Flux vector
\param[in] S1   First wave speed
\param[in] S2   Second wave speed
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC void Gasdyn<fluid>::ConservedState::FixMomentum(Gasdyn<fluid>::FluxFunction& flux, double S1, double S2)
{
   mom()[1] = flux.mom()[1] / (S1 - S2);
   mom()[2] = flux.mom()[2] / (S1 - S2);
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] flux Flux vector
\param[in] S1   First wave speed
\param[in] S2   Second wave speed
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC void Gasdyn<fluid>::ConservedState::FixEnergy(Gasdyn<fluid>::FluxFunction& flux, double S1, double S2)
{
   double pre_total = S1 * mom()[0] - flux.mom()[0];
   enr() = (flux.enr() + pre_total * S2) / (S1 - S2);
};

/*!
\author Vladimir Florinski
\date 01/15/2025
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC void Gasdyn<fluid>::ConservedState::Invert(void)
{
   mom()[0] = -mom()[0];
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Gasdyn::FluxFunction methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] denf_in Density flux
\param[in] momf_in Momentum flux
\param[in] enrf_in Energy flux
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC Gasdyn<fluid>::FluxFunction::FluxFunction(double denf_in, const GeoVector& momf_in, double enrf_in)
{
   denf() = denf_in;
   momf() = momf_in;
   enrf() = enrf_in;
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] other Object to initialize from
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC Gasdyn<fluid>::FluxFunction::FluxFunction(const SimpleArray<double, CL_GASDYN_NVARS>& other)
                                                : SimpleArray<double, CL_GASDYN_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 01/15/2025
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC void Gasdyn<fluid>::FluxFunction::Invert(void)
{
   momf()[0] = -momf()[0];
};

};
