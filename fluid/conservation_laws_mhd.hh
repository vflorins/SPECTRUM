/*!
\file conservation_laws_mhd.hh
\brief Declares rules for ideal MHD equations
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_CONSERVATION_LAWS_MHD_HH
#define SPECTRUM_CONSERVATION_LAWS_MHD_HH

#include "common/physics.hh"

namespace Spectrum {

//! Number of GLM variables
#define CL_MHD_GLM

#ifdef CL_MHD_GLM
//! Number of fluid variables (DEN,VELx3,PRE,MAGx3, GLM)
#define CL_MHD_NVARS 9
//! Safety factor for GLM transport
#define CL_MHD_GLMSAFETY 1.05

#else
//! Number of fluid variables (DEN,VELx3,PRE,MAGx3)
#define CL_MHD_NVARS 8
#endif

template <int fluid> struct PrimitiveStateMHD;
template <int fluid> struct ConservedStateMHD;
template <int fluid> struct FluxFunctionMHD;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// PrimitiveStateMHD class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Codifies physics for ideal single fluid MHD
\author Vladimir Florinski
*/
template <int fluid>
struct PrimitiveStateMHD : public SimpleArray<double, CL_MHD_NVARS>
{
   using SimpleArray<double, CL_MHD_NVARS>::data;

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

//! Alias for magnetic field (RO)
   const GeoVector& mag(void) const {return (const GeoVector&)data[5];};

//! Alias for magnetic field (RW)
   GeoVector& mag(void) {return (GeoVector&)data[5];};

#ifdef CL_MHD_GLM

//! Alias for Largange multiplier (RO)
   const double& glm(void) const {return data[8];};

//! Alias for Largange multiplier (RW)
   double& glm(void) {return data[8];};

#endif

//! Return the "fluid" template parameter
   SPECTRUM_DEVICE_FUNC static constexpr int Fluid(void) {return fluid;};

//! Return the number of variables
   SPECTRUM_DEVICE_FUNC static constexpr int Nvars(void) {return CL_MHD_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC PrimitiveStateMHD(void) = default;

#ifdef CL_MHD_GLM
//! Constructor from components
   SPECTRUM_DEVICE_FUNC PrimitiveStateMHD(double den_in, const GeoVector& vel_in, double pre_in, const GeoVector& mag_in, double glm_in);
#else
//! Constructor from components
   SPECTRUM_DEVICE_FUNC PrimitiveStateMHD(double den_in, const GeoVector& vel_in, double pre_in, const GeoVector& mag_in);
#endif

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC PrimitiveStateMHD(const SimpleArray<double, CL_MHD_NVARS>& other);

//! Calculate the fastest wave speed
   SPECTRUM_DEVICE_FUNC double FastestWave(void) const;

//! Calculate the fastest normal wave speed
   SPECTRUM_DEVICE_FUNC double FastestWaveNormal(void) const;

//! Calculate primitive state
   SPECTRUM_DEVICE_FUNC ConservedStateMHD<fluid> ToConserved(void) const;

//! Calculate flux function
   SPECTRUM_DEVICE_FUNC FluxFunctionMHD<fluid> ToFlux(void) const;

//! Calculate sonic state
   SPECTRUM_DEVICE_FUNC PrimitiveStateMHD<fluid> ToSonic(void) const;

//! Invert vector variables
   SPECTRUM_DEVICE_FUNC void Invert(void);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ConservedStateMHD class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Codifies physics for ideal single fluid MHD
\author Vladimir Florinski
*/
template <int fluid>
struct ConservedStateMHD : public SimpleArray<double, CL_MHD_NVARS>
{
   using SimpleArray<double, CL_MHD_NVARS>::data;

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

//! Alias for magnetic field (RO)
   const GeoVector& mag(void) const {return (const GeoVector&)data[5];};

//! Alias for magnetic field (RW)
   GeoVector& mag(void) {return (GeoVector&)data[5];};

#ifdef CL_MHD_GLM

//! Alias for Largange multiplier (RO)
   const double& glm(void) const {return data[8];};

//! Alias for Largange multiplier (RW)
   double& glm(void) {return data[8];};

#endif

//! Return the "fluid" template parameter
   SPECTRUM_DEVICE_FUNC static constexpr int Fluid(void) {return fluid;};

//! Return the number of variables
   SPECTRUM_DEVICE_FUNC static constexpr int Nvars(void) {return CL_MHD_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC ConservedStateMHD(void) = default;

#ifdef CL_MHD_GLM
//! Constructor from components
   SPECTRUM_DEVICE_FUNC ConservedStateMHD(double den_in, const GeoVector& mom_in, double enr_in, const GeoVector& mag_in, double glm_in);
#else
//! Constructor from components
   SPECTRUM_DEVICE_FUNC ConservedStateMHD(double den_in, const GeoVector& mom_in, double enr_in, const GeoVector& mag_in);
#endif

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC ConservedStateMHD(const SimpleArray<double, CL_MHD_NVARS>& other);

//! Calculate the fastest wave speed
   SPECTRUM_DEVICE_FUNC double FastestWave(void) const;

//! Calculate the fastest normal wave speed
   SPECTRUM_DEVICE_FUNC double FastestWaveNormal(void) const;

//! Calculate primitive state
   SPECTRUM_DEVICE_FUNC PrimitiveStateMHD<fluid> ToPrimitive(void) const;

//! Calculate flux function
   SPECTRUM_DEVICE_FUNC FluxFunctionMHD<fluid> ToFlux(void) const;

//! Calculate _gas_ pressure
   SPECTRUM_DEVICE_FUNC double GetPressure(void) const;

//! Calculate y and z components of momentum from the externally provided flux in a moving frame
   SPECTRUM_DEVICE_FUNC void FixMomentum(FluxFunctionMHD<fluid>& flux, double S1, double S2);

//! Calculate energy from the externally provided flux in a moving frame
   SPECTRUM_DEVICE_FUNC void FixEnergy(FluxFunctionMHD<fluid>& flux, double S1, double S2);

//! Invert vector variables
   SPECTRUM_DEVICE_FUNC void Invert(void);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// FluxFunctionMHD class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Codifies physics for ideal single fluid MHD
\author Vladimir Florinski
*/
template <int fluid>
struct FluxFunctionMHD : public SimpleArray<double, CL_MHD_NVARS>
{
   using SimpleArray<double, CL_MHD_NVARS>::data;

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

//! Alias for magnetic flux (RO)
   const GeoVector& magf(void) const {return (const GeoVector&)data[5];};

//! Alias for magnetic flux (RW)
   GeoVector& magf(void) {return (GeoVector&)data[5];};

#ifdef CL_MHD_GLM

//! Alias for Largange multiplier flux (RO)
   const double& glmf(void) const {return data[8];};

//! Alias for Largange multiplier flux (RO)
   double& glmf(void) {return data[8];};

#endif

//! Return the "fluid" template parameter
   SPECTRUM_DEVICE_FUNC static constexpr int Fluid(void) {return fluid;};

//! Return the number of variables
   SPECTRUM_DEVICE_FUNC static constexpr int Nvars(void) {return CL_MHD_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC FluxFunctionMHD(void) = default;

#ifdef CL_MHD_GLM
//! Constructor from components
   SPECTRUM_DEVICE_FUNC FluxFunctionMHD(double denf_in, const GeoVector& momf_in, double enrf_in, const GeoVector& magf_in, double glmf_in);
#else
//! Constructor from components
   SPECTRUM_DEVICE_FUNC FluxFunctionMHD(double denf_in, const GeoVector& momf_in, double enrf_in, const GeoVector& magf_in);
#endif

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC FluxFunctionMHD(const SimpleArray<double, CL_MHD_NVARS>& other);

//! Invert vector variables
   SPECTRUM_DEVICE_FUNC void Invert(void);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// PrimitiveStateMHD inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 03/14/2024
\param[in] den_in Density
\param[in] vel_in Velocity
\param[in] pre_in Pressure
\param[in] mag_in Magnetic field
\param[in] glm_in Largange multiplier
*/
template <int fluid>
#ifdef CL_MHD_GLM
SPECTRUM_DEVICE_FUNC inline PrimitiveStateMHD<fluid>::PrimitiveStateMHD(double den_in, const GeoVector& vel_in, double pre_in, const GeoVector& mag_in, double glm_in)
#else
SPECTRUM_DEVICE_FUNC inline PrimitiveStateMHD<fluid>::PrimitiveStateMHD(double den_in, const GeoVector& vel_in, double pre_in, const GeoVector& mag_in)
#endif
{
   den() = den_in;
   vel() = vel_in;
   pre() = pre_in;
   mag() = mag_in;

#ifdef CL_MHD_GLM
   glm() = glm_in;
#endif
};

/*!
\author Vladimir Florinski
\date 03/14/2024
\param[in] other Object to initialize from
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline PrimitiveStateMHD<fluid>::PrimitiveStateMHD(const SimpleArray<double, CL_MHD_NVARS>& other)
                                                    : SimpleArray<double, CL_MHD_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 03/16/2024
\return Fastest wave speed (max fast with GLM safety)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline double PrimitiveStateMHD<fluid>::FastestWave(void) const
{
   double Va2 = Alfven2(den(), mag().Norm2());
   double Cs2 = Sound2(den(), pre(), fluid);

// FIXME decide where to use "CL_MHD_GLMSAFETY"
//#ifdef CL_MHD_GLM
   return sqrt(Va2 + Cs2);
//#else
//   return CL_MHD_GLMSAFETY * sqrt(Va2 + Cs2);
//#endif
};

/*!
\author Vladimir Florinski
\date 03/16/2024
\return Fastest normal wave speed (fast_n)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline double PrimitiveStateMHD<fluid>::FastestWaveNormal(void) const
{
   double Va2 = Alfven2(den(), mag().Norm2());
   double Vax2 = Alfven2(den(), Sqr(mag()[0]));
   double Cs2 = Sound2(den(), pre(), fluid);

   return FastMagnetosonicSpeed(Va2, Vax2, Cs2);
};

/*!
\author Vladimir Florinski
\date 03/06/2024
\return Vector of conserved variables
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline ConservedStateMHD<fluid> PrimitiveStateMHD<fluid>::ToConserved(void) const
{
   ConservedStateMHD<fluid> cons;

   cons.den() = den();
   cons.mom() = vel() * den();
   cons.enr() = Energy(den(), vel().Norm2(), mag().Norm2(), pre(), fluid);
   cons.mag() = mag();

#ifdef CL_MHD_GLM
   cons.glm() = glm();
#endif

   return cons;
};

/*!
\author Vladimir Florinski
\date 03/06/2024
\return Flux vector
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline FluxFunctionMHD<fluid> PrimitiveStateMHD<fluid>::ToFlux(void) const
{
   double enr = Energy(den(), vel().Norm2(), mag().Norm2(), pre(), fluid);
   double pre_tot = pre() + mag().Norm2() / M_8PI;
   FluxFunctionMHD<fluid> flux;

   flux.denf() = den() * vel()[0];
   flux.momf() = den() * vel()[0] * vel() + pre_tot * gv_nx - (mag()[0] / M_4PI) * mag();
   flux.enrf() = (enr + pre_tot) * vel()[0] - (vel() * mag()) * mag()[0] / M_4PI;
   flux.magf() = vel()[0] * mag() - mag()[0] * vel();

#ifdef CL_MHD_GLM
// The GLM flux is proportional to the square of the fastest speed. The safety factor "CL_MHD_GLMSAFETY" is used to make sure it is slightly outside of the "normal" Riemann fan.
   flux.magf() += glm() * gv_nx;

// FIXME should this flux be calculated in the RS?
//   flux.glmf() = Sqr(CL_MHD_GLMSAFETY * FastestWaveNormal()) * mag()[0];
#endif

   return flux;
};

/*!
\author Vladimir Florinski
\date 04/24/2024
\return Vector of primitive variables in a sonic state
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline PrimitiveStateMHD<fluid> PrimitiveStateMHD<fluid>::ToSonic(void) const
{
   GeoVector aa = mag() / sqrt(M_4PI * den());
   double ax2 = Sqr(aa[0]);
   double at2 = Sqr(aa[1]) + Sqr(aa[2]);
   double c2 = Sound2(den(), pre(), fluid);
   double c = sqrt(c2);
   double s2 = c2 + ax2 + at2;
   double d2 = sqrt(s2 - 4.0 * c2 * ax2);
   double af2 = 0.5 * (s2 + d2);
   double af = sqrt(af2);
   double as2 = af2 - d2;
   double den_prime = den() * fmax(1.0 - as2 / c2, small);

   double q2 = 0.25 * den() / af2 * ((gamma_eos[fluid] * c2 - s2) / den() + 2.0 * at2 * (d2 + s2) / d2 / den_prime
      + (gamma_eos[fluid] * c2 * (s2 - 2.0 * ax2) - s2 * s2) / d2 / den());
   double del_den = den() / (1.0 + q2) * (vel()[0] / af - 1.0);
   double af_star = vel()[0] - af * del_den / den();
   double cs_star = vel()[0] - c  * del_den / den();

   PrimitiveStateMHD<fluid> prim;
   
   prim.den() = den() + del_den;
   prim.vel()[0] = af_star;
   prim.vel()[1] = vel()[1] + sqrt(as2) * aa[1] * sign(mag()[0]) * del_den / c / den_prime;
   prim.vel()[2] = vel()[2] + sqrt(as2) * aa[2] * sign(mag()[0]) * del_den / c / den_prime;
   prim.mag()[0] = mag()[0];
   prim.mag()[1] = mag()[1] * (1.0 + del_den / den_prime);
   prim.mag()[2] = mag()[2] * (1.0 + del_den / den_prime);
   prim.pre() = prim.den() * Sqr(cs_star) / gamma_eos[fluid];
   return prim;
};

/*!
\author Vladimir Florinski
\date 04/09/2024
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline void PrimitiveStateMHD<fluid>::Invert(void)
{
   vel()[0] = -vel()[0];
   mag()[0] = -mag()[0];
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ConservedStateMHD inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 03/14/2024
\param[in] den_in Density
\param[in] mom_in Momentum
\param[in] enr_in Energy
\param[in] mag_in Magnetic field
\param[in] glm_in Largange multiplier
*/
template <int fluid>
#ifdef CL_MHD_GLM
SPECTRUM_DEVICE_FUNC inline ConservedStateMHD<fluid>::ConservedStateMHD(double den_in, const GeoVector& mom_in, double enr_in, const GeoVector& mag_in, double glm_in)
#else
SPECTRUM_DEVICE_FUNC inline ConservedStateMHD<fluid>::ConservedStateMHD(double den_in, const GeoVector& mom_in, double enr_in, const GeoVector& mag_in)
#endif
{
   den() = den_in;
   mom() = mom_in;
   enr() = enr_in;
   mag() = mag_in;

#ifdef CL_MHD_GLM
   glm() = glm_in;
#endif
};

/*!
\author Vladimir Florinski
\date 03/14/2024
\param[in] other Object to initialize from
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline ConservedStateMHD<fluid>::ConservedStateMHD(const SimpleArray<double, CL_MHD_NVARS>& other)
                                                    : SimpleArray<double, CL_MHD_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 03/16/2024
\return Fastest wave speed (max fast with GLM safety)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline double ConservedStateMHD<fluid>::FastestWave(void) const
{
   GeoVector vel = mom() / den();
   double pre = Pressure(den(), vel.Norm2(), mag().Norm2(), enr(), fluid);
   double Va2 = Alfven2(den(), mag().Norm2());
   double Cs2 = Sound2(den(), pre, fluid);

// FIXME decide where to use "CL_MHD_GLMSAFETY"
//#ifdef CL_MHD_GLM
   return sqrt(Va2 + Cs2);
//#else
//   return CL_MHD_GLMSAFETY * sqrt(Va2 + Cs2);
//#endif
};

/*!
\author Vladimir Florinski
\date 03/16/2024
\return Fastest normal wave speed (fast_n)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline double ConservedStateMHD<fluid>::FastestWaveNormal(void) const
{
   GeoVector vel = mom() / den();
   double pre = Pressure(den(), vel.Norm2(), mag().Norm2(), enr(), fluid);
   double Va2 = Alfven2(den(), mag().Norm2());
   double Vax2 = Alfven2(den(), Sqr(mag()[0]));
   double Cs2 = Sound2(den(), pre, fluid);

   return FastMagnetosonicSpeed(Va2, Vax2, Cs2);
};

/*!
\author Vladimir Florinski
\date 03/06/2024
\return Vector of primitive variables
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline PrimitiveStateMHD<fluid> ConservedStateMHD<fluid>::ToPrimitive(void) const
{
   PrimitiveStateMHD<fluid> prim;

   prim.den() = den();
   prim.vel() = mom() / den();
   prim.pre() = Pressure(den(), prim.vel().Norm2(), mag().Norm2(), enr(), fluid);
   prim.mag() = mag();

#ifdef CL_MHD_LAGR
   prim.glm() = glm();
#endif

   return prim;
};

/*!
\author Vladimir Florinski
\date 03/14/2024
\return Flux vector
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline FluxFunctionMHD<fluid> ConservedStateMHD<fluid>::ToFlux(void) const
{
   GeoVector vel = mom() / den();
   double pre = Pressure(den(), vel.Norm2(), mag().Norm2(), enr(), fluid);
   double pre_tot = pre + mag().Norm2() / M_8PI;
   FluxFunctionMHD<fluid> flux;

   flux.denf() = mom()[0];
   flux.momf() = mom()[0] * vel + pre_tot * gv_nx - (mag()[0] / M_4PI) * mag();
   flux.enrf() = (enr() + pre_tot) * vel[0] - (vel * mag()) * mag()[0] / M_4PI;
   flux.magf() = vel[0] * mag() - mag()[0] * vel;

#ifdef CL_MHD_GLM
// The GLM flux is proportional to the square of the fastest speed. The safety factor "CL_MHD_GLMSAFETY" is used to make sure it is slightly outside of the "normal" Riemann fan.
   flux.magf() += glm() * gv_nx;

// FIXME should this flux be calculated in the RS?
//   flux.glmf() = Sqr(CL_MHD_GLMSAFETY * FastestWaveNormal()) * mag()[0];
#endif

   return flux;
};

/*!
\author Vladimir Florinski
\date 03/20/2024
\return Gas pressure
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline double ConservedStateMHD<fluid>::GetPressure(void) const
{
   return Pressure(den(), mom().Norm2() / Sqr(den()), mag().Norm2(), enr(), fluid);
};

/*!
\author Vladimir Florinski
\date 03/20/2024
\param[in] flux Flux vector
\param[in] S1   First wave speed
\param[in] S2   Second wave speed
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline void ConservedStateMHD<fluid>::FixMomentum(FluxFunctionMHD<fluid>& flux, double S1, double S2)
{
   mom()[1] = (flux.mom()[1] - mag()[0] * mag()[1] / M_4PI) / (S1 - S2);
   mom()[2] = (flux.mom()[2] - mag()[0] * mag()[2] / M_4PI) / (S1 - S2);
};

/*!
\author Vladimir Florinski
\date 03/20/2024
\param[in] flux Flux vector
\param[in] S1   First wave speed
\param[in] S2   Second wave speed
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline void ConservedStateMHD<fluid>::FixEnergy(FluxFunctionMHD<fluid>& flux, double S1, double S2)
{
   double pre_total = S1 * mom()[0] + mag()[0] * mag()[0] / M_4PI - flux.mom()[0];
   enr() = (flux.enr() - (mag()[0] * S2 + (mag()[1] * mom()[1] + mag()[2] * mom()[2]) / den()) * mag()[0] + pre_total * S2) / (S1 - S2);
};

/*!
\author Vladimir Florinski
\date 04/09/2024
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline void ConservedStateMHD<fluid>::Invert(void)
{
   mom()[0] = -mom()[0];
   mag()[0] = -mag()[0];
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// FluxFunctionMHD inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 03/15/2024
\param[in] denf_in Density flux
\param[in] momf_in Momentum flux
\param[in] enrf_in Energy flux
\param[in] magf_in Magnetic field
\param[in] glmf_in Largange multiplier flux
*/
template <int fluid>
#ifdef CL_MHD_GLM
SPECTRUM_DEVICE_FUNC inline FluxFunctionMHD<fluid>::FluxFunctionMHD(double denf_in, const GeoVector& momf_in, double enrf_in, const GeoVector& magf_in, double glmf_in)
#else
SPECTRUM_DEVICE_FUNC inline FluxFunctionMHD<fluid>::FluxFunctionMHD(double denf_in, const GeoVector& momf_in, double enrf_in, const GeoVector& magf_in)
#endif
{
   denf() = denf_in;
   momf() = momf_in;
   enrf() = enrf_in;
   magf() = magf_in;

#ifdef CL_MHD_GLM
   glmf() = glmf_in;
#endif
};

/*!
\author Vladimir Florinski
\date 03/14/2024
\param[in] other Object to initialize from
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline FluxFunctionMHD<fluid>::FluxFunctionMHD(const SimpleArray<double, CL_MHD_NVARS>& other)
                                                  : SimpleArray<double, CL_MHD_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 04/09/2024
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC inline void FluxFunctionMHD<fluid>::Invert(void)
{
   momf()[0] = -momf()[0];
   magf()[0] = -magf()[0];
};

};

#endif
