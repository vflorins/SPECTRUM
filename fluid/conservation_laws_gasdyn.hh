/*!
\file conservation_laws_gasdyn.hh
\brief Declares rules for compressible gas-dynamic equations
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_CONSERVATION_LAWS_GASDYN_HH
#define SPECTRUM_CONSERVATION_LAWS_GASDYN_HH

#include "common/physics.hh"

namespace Spectrum {

//! Number of fluid variables (DEN,VELx3,PRE)
#define CL_GASDYN_NVARS 5

template <int fluid> struct PrimitiveStateGasdyn;
template <int fluid> struct ConservedStateGasdyn;
template <int fluid> struct FluxFunctionGasdyn;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// PrimitiveStateGasdyn class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Codifies physics for compressible gas dynamics
\author Vladimir Florinski
*/
template <int fluid>
struct PrimitiveStateGasdyn : public SimpleArray<double, CL_GASDYN_NVARS>
{
   using SimpleArray<double, CL_GASDYN_NVARS>::data;

//! A trait to be used in template specializations
   static constexpr bool is_cons_law_active = true;

//! Alias for density (RO)
   const double& den(void) const {return data[0];};

//! Alias for density (RW)
   double& den(void) {return data[0];};

//! Alias for velocity (RO)
   const GeoVector& vel(void) const {return (const GeoVector&)data[1];};

//! Alias for velocity (RW)
   GeoVector& vel(void) {return (GeoVector&)data[1];};

//! Alias for pressure (RO)
   const double& pre(void) const {return data[4];};

//! Alias for pressure (RW)
   double& pre(void) {return data[4];};

//! Return the "fluid" template parameter
   SPECTRUM_DEVICE_FUNC static constexpr int Fluid(void) {return fluid;};

//! Return the number of variables
   SPECTRUM_DEVICE_FUNC static constexpr int Nvars(void) {return CL_GASDYN_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC PrimitiveStateGasdyn(void) = default;

//! Constructor from components
   SPECTRUM_DEVICE_FUNC PrimitiveStateGasdyn(double den_in, const GeoVector& vel_in, double pre_in);

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC PrimitiveStateGasdyn(const SimpleArray<double, CL_GASDYN_NVARS>& other);

//! Calculate the fastest wave speed
   SPECTRUM_DEVICE_FUNC double FastestWave(void) const;

//! Calculate the fastest normal wave speed
   SPECTRUM_DEVICE_FUNC double FastestWaveNormal(void) const;

//! Calculate conserved state
   SPECTRUM_DEVICE_FUNC ConservedStateGasdyn<fluid> ToConserved(void) const;

//! Calculate flux function
   SPECTRUM_DEVICE_FUNC FluxFunctionGasdyn<fluid> ToFlux(void) const;

//! Calculate sonic state
   SPECTRUM_DEVICE_FUNC PrimitiveStateGasdyn<fluid> ToSonic(void) const;

//! Invert vector variables
   SPECTRUM_DEVICE_FUNC void Invert(void);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ConservedStateGasdyn class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Codifies physics for compressible gas dynamics
\author Vladimir Florinski
*/
template <int fluid>
struct ConservedStateGasdyn : public SimpleArray<double, CL_GASDYN_NVARS>
{
   using SimpleArray<double, CL_GASDYN_NVARS>::data;

//! A trait to be used in template specializations
   static constexpr bool is_cons_law_active = true;

//! Alias for density (RO)
   const double& den(void) const {return data[0];};

//! Alias for density (RW)
   double& den(void) {return data[0];};

//! Alias for momentum (RO)
   const GeoVector& mom(void) const {return (const GeoVector&)data[1];};

//! Alias for momentum (RW)
   GeoVector& mom(void) {return (GeoVector&)data[1];};

//! Alias for energy (RO)
   const double& enr(void) const {return data[4];};

//! Alias for energy (RW)
   double& enr(void) {return data[4];};

//! Return the "fluid" template parameter
   SPECTRUM_DEVICE_FUNC static constexpr int Fluid(void) {return fluid;};

//! Return the number of variables
   SPECTRUM_DEVICE_FUNC static constexpr int Nvars(void) {return CL_GASDYN_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC ConservedStateGasdyn(void) = default;

//! Constructor from components
   SPECTRUM_DEVICE_FUNC ConservedStateGasdyn(double den_in, const GeoVector& mom_in, double enr_in);

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC ConservedStateGasdyn(const SimpleArray<double, CL_GASDYN_NVARS>& other);

//! Calculate the fastest wave speed
   SPECTRUM_DEVICE_FUNC double FastestWave(void) const;

//! Calculate the fastest normal wave speed
   SPECTRUM_DEVICE_FUNC double FastestWaveNormal(void) const;

//! Calculate primitive state
   SPECTRUM_DEVICE_FUNC PrimitiveStateGasdyn<fluid> ToPrimitive(void) const;

//! Calculate flux function
   SPECTRUM_DEVICE_FUNC FluxFunctionGasdyn<fluid> ToFlux(void) const;

//! Calculate _gas_ pressure
   SPECTRUM_DEVICE_FUNC double GetPressure(void) const;

//! Calculate y and z components of momentum from the externally provided flux in a moving frame
   SPECTRUM_DEVICE_FUNC void FixMomentum(FluxFunctionGasdyn<fluid>& flux, double S1, double S2);

//! Calculate energy from the externally provided flux in a moving frame
   SPECTRUM_DEVICE_FUNC void FixEnergy(FluxFunctionGasdyn<fluid>& flux, double S1, double S2);

//! Invert vector variables
   SPECTRUM_DEVICE_FUNC void Invert(void);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// FluxFunctionGasdyn class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Codifies physics for compressible gas dynamics
\author Vladimir Florinski
*/
template <int fluid>
struct FluxFunctionGasdyn : public SimpleArray<double, CL_GASDYN_NVARS>
{
   using SimpleArray<double, CL_GASDYN_NVARS>::data;

//! A trait to be used in template specializations
   static constexpr bool is_cons_law_active = true;

//! Alias for density flux (RO)
   const double& denf(void) const {return data[0];};

//! Alias for density flux (RW)
   double& denf(void) {return data[0];};

//! Alias for momentum flux (RO)
   const GeoVector& momf(void) const {return (const GeoVector&)data[1];};

//! Alias for momentum flux (RW)
   GeoVector& momf(void) {return (GeoVector&)data[1];};

//! Alias for energy flux (RO)
   const double& enrf(void) const {return data[4];};

//! Alias for energy flux (RW)
   double& enrf(void) {return data[4];};

//! Return the "fluid" template parameter
   SPECTRUM_DEVICE_FUNC static constexpr int Fluid(void) {return fluid;};

//! Return the number of variables
   SPECTRUM_DEVICE_FUNC static constexpr int Nvars(void) {return CL_GASDYN_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC FluxFunctionGasdyn(void) = default;

//! Constructor from components
   SPECTRUM_DEVICE_FUNC FluxFunctionGasdyn(double denf_in, const GeoVector& momf_in, double enrf_in);

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC FluxFunctionGasdyn(const SimpleArray<double, CL_GASDYN_NVARS>& other);

//! Invert vector variables
   SPECTRUM_DEVICE_FUNC void Invert(void);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// PrimitiveStateGasdyn inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 03/15/2024
\param[in] den_in Density
\param[in] vel_in Velocity
\param[in] pre_in Pressure
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline PrimitiveStateGasdyn<fluid>::PrimitiveStateGasdyn(double den_in, const GeoVector&vel_in, double pre_in)
{
   den() = den_in;
   vel() = vel_in;
   pre() = pre_in;
};

/*!
\author Vladimir Florinski
\date 03/14/2024
\param[in] other Object to initialize from
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline PrimitiveStateGasdyn<fluid>::PrimitiveStateGasdyn(const SimpleArray<double, CL_GASDYN_NVARS>& other)
                                                       : SimpleArray<double, CL_GASDYN_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 03/16/2024
\return Fastest wave speed (sound)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline double PrimitiveStateGasdyn<fluid>::FastestWave(void) const
{
   return SoundSpeed(den(), pre(), fluid);
};

/*!
\author Vladimir Florinski
\date 03/16/2024
\return Fastest normal wave speed (sound)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline double PrimitiveStateGasdyn<fluid>::FastestWaveNormal(void) const
{
   return SoundSpeed(den(), pre(), fluid);
};

/*!
\author Vladimir Florinski
\date 03/06/2024
\return Vector of conserved variables
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline ConservedStateGasdyn<fluid> PrimitiveStateGasdyn<fluid>::ToConserved(void) const
{
   ConservedStateGasdyn<fluid> cons;

   cons.den() = den();
   cons.mom() = den() * vel();
   cons.enr() = Energy(den(), vel().Norm2(), 0.0, pre(), fluid);

   return cons;
};

/*!
\author Vladimir Florinski
\date 03/06/2024
\return Flux vector
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline FluxFunctionGasdyn<fluid> PrimitiveStateGasdyn<fluid>::ToFlux(void) const
{
   double enr = Energy(den(), vel().Norm2(), 0.0, pre(), fluid);
   FluxFunctionGasdyn<fluid> flux;

   flux.denf() = den() * vel()[0];
   flux.momf() = den() * vel()[0] * vel() + pre() * gv_nx;
   flux.enrf() = (enr + pre()) * vel()[0];

   return flux;
};

/*!
\author Vladimir Florinski
\date 04/24/2024
\return Vector of primitive variables in a sonic state
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline PrimitiveStateGasdyn<fluid> PrimitiveStateGasdyn<fluid>::ToSonic(void) const
{
   double cs = Sound(den(), pre(), fluid);
   double del_den = 2.0 * den() / (gamma_eos[fluid] + 1.0) * (vel()[0] / cs - 1.0);
   double cs_star = vel()[0] - cs * del_den / den();

   PrimitiveStateGasdyn<fluid> prim;
   
   prim.den() = den() + del_den;
   prim.vel()[0] = cs_star;
   prim.vel()[1] = vel()[1];
   prim.vel()[2] = vel()[2];
   prim.pre() = prim.den() * Sqr(cs_star) / gamma_eos[fluid];
   return prim;
};

/*!
\author Vladimir Florinski
\date 04/09/2024
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline void PrimitiveStateGasdyn<fluid>::Invert(void)
{
   vel()[0] = -vel()[0];
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ConservedStateGasdyn inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 03/14/2024
\param[in] den_in Density
\param[in] mom_in Momentum
\param[in] enr_in Energy
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline ConservedStateGasdyn<fluid>::ConservedStateGasdyn(double den_in, const GeoVector& mom_in, double enr_in)
{
   den() = den_in;
   mom() = mom_in;
   enr() = enr_in;
};

/*!
\author Vladimir Florinski
\date 03/14/2024
\param[in] other Object to initialize from
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline ConservedStateGasdyn<fluid>::ConservedStateGasdyn(const SimpleArray<double, CL_GASDYN_NVARS>& other)
                                                       : SimpleArray<double, CL_GASDYN_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 03/16/2024
\return Fastest wave speed (sound)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline double ConservedStateGasdyn<fluid>::FastestWave(void) const
{
   GeoVector vel = mom() / den();
   double pre = Pressure(den(), vel.Norm2(), 0.0, enr(), fluid);
   return SoundSpeed(den(), pre, fluid);
};

/*!
\author Vladimir Florinski
\date 03/16/2024
\return Fastest normal wave speed (sound)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline double ConservedStateGasdyn<fluid>::FastestWaveNormal(void) const
{
   GeoVector vel = mom() / den();
   double pre = Pressure(den(), vel.Norm2(), 0.0, enr(), fluid);
   return SoundSpeed(den(), pre, fluid);
};

/*!
\author Vladimir Florinski
\date 03/06/2024
\return Vector of primitive variables
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline PrimitiveStateGasdyn<fluid> ConservedStateGasdyn<fluid>::ToPrimitive(void) const
{
   PrimitiveStateGasdyn<fluid> prim;

   prim.den() = den();
   prim.vel() = mom() / den();
   prim.pre() = Pressure(den(), prim.vel().Norm2(), 0.0, enr(), fluid);

   return prim;
};

/*!
\author Vladimir Florinski
\date 03/06/2024
\return Flux vector
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline FluxFunctionGasdyn<fluid> ConservedStateGasdyn<fluid>::ToFlux(void) const
{
   double pre = Pressure(den(), mom().Norm2() / Sqr(den()), 0.0, enr(), fluid);

   FluxFunctionGasdyn<fluid> flux;
   flux.denf() = mom()[0];
   flux.momf() = mom()[0] * mom() / den() + pre * gv_nx;
   flux.enrf() = (enr() + pre) * mom()[0] / den();

   return flux;
};

/*!
\author Vladimir Florinski
\date 03/20/2024
\return Gas pressure
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline double ConservedStateGasdyn<fluid>::GetPressure(void) const
{
   return Pressure(den(), mom().Norm2() / Sqr(den()), 0.0, enr(), fluid);
};

/*!
\author Vladimir Florinski
\date 03/18/2024
\param[in] flux Flux vector
\param[in] S1   First wave speed
\param[in] S2   Second wave speed
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline void ConservedStateGasdyn<fluid>::FixMomentum(FluxFunctionGasdyn<fluid>& flux, double S1, double S2)
{
   mom()[1] = flux.mom()[1] / (S1 - S2);
   mom()[2] = flux.mom()[2] / (S1 - S2);
};

/*!
\author Vladimir Florinski
\date 03/20/2024
\param[in] flux Flux vector
\param[in] S1   First wave speed
\param[in] S2   Second wave speed
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline void ConservedStateGasdyn<fluid>::FixEnergy(FluxFunctionGasdyn<fluid>& flux, double S1, double S2)
{
   double pre_total = S1 * mom()[0] - flux.mom()[0];
   enr() = (flux.enr() + pre_total * S2) / (S1 - S2);
};

/*!
\author Vladimir Florinski
\date 04/09/2024
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline void ConservedStateGasdyn<fluid>::Invert(void)
{
   mom()[0] = -mom()[0];
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// FluxFunctionGasdyn inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 03/15/2024
\param[in] denf_in Density flux
\param[in] momf_in Momentum flux
\param[in] enrf_in Energy flux
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline FluxFunctionGasdyn<fluid>::FluxFunctionGasdyn(double denf_in, const GeoVector& momf_in, double enrf_in)
{
   denf() = denf_in;
   momf() = momf_in;
   enrf() = enrf_in;
};

/*!
\author Vladimir Florinski
\date 03/14/2024
\param[in] other Object to initialize from
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline FluxFunctionGasdyn<fluid>::FluxFunctionGasdyn(const SimpleArray<double, CL_GASDYN_NVARS>& other)
                                                     : SimpleArray<double, CL_GASDYN_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 04/09/2024
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline void FluxFunctionGasdyn<fluid>::Invert(void)
{
   momf()[0] = -momf()[0];
};

};

#endif
