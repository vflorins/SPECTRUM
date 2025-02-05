/*!
\file conservation_laws_mhd.cc
\brief Implements rules for ideal MHD equations
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "common/physics.hh"
#include "fluid/conservation_laws_mhd.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// MHD::PrimitiveState methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] den_in Density
\param[in] vel_in Velocity
\param[in] pre_in Pressure
\param[in] mag_in Magnetic field
\param[in] glm_in Largange multiplier
*/
template <int fluid>
#ifdef CL_MHD_GLM
SPECTRUM_DEVICE_FUNC MHD<fluid>::PrimitiveState::PrimitiveState(double den_in, const GeoVector& vel_in, double pre_in, const GeoVector& mag_in, double glm_in)
#else
SPECTRUM_DEVICE_FUNC MHD<fluid>::PrimitiveState::PrimitiveState(double den_in, const GeoVector& vel_in, double pre_in, const GeoVector& mag_in)
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
\date 01/15/2025
\param[in] other Object to initialize from
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC MHD<fluid>::PrimitiveState::PrimitiveState(const SimpleArray<double, CL_MHD_NVARS>& other)
                                               : SimpleArray<double, CL_MHD_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Fastest wave speed (max fast with GLM safety)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC double MHD<fluid>::PrimitiveState::FastestWave(void) const
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
\date 01/15/2025
\return Fastest normal wave speed (fast_n)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC double MHD<fluid>::PrimitiveState::FastestWaveNormal(void) const
{
   double Va2 = Alfven2(den(), mag().Norm2());
   double Vax2 = Alfven2(den(), Sqr(mag()[0]));
   double Cs2 = Sound2(den(), pre(), fluid);

   return FastMagnetosonicSpeed(Va2, Vax2, Cs2);
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Vector of conserved variables
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC MHD<fluid>::ConservedState MHD<fluid>::PrimitiveState::ToConserved(void) const
{
   MHD<fluid>::ConservedState cons;

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
\date 01/15/2025
\return Flux vector
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC MHD<fluid>::FluxFunction MHD<fluid>::PrimitiveState::ToFlux(void) const
{
   double enr = Energy(den(), vel().Norm2(), mag().Norm2(), pre(), fluid);
   double pre_tot = pre() + mag().Norm2() / M_8PI;
   MHD<fluid>::FluxFunction flux;

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
\date 01/15/2025
\return Vector of primitive variables in a sonic state
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC MHD<fluid>::PrimitiveState MHD<fluid>::PrimitiveState::ToSonic(void) const
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
   double den_prime = den() * fmax(1.0 - as2 / c2, sp_small);

   double q2 = 0.25 * den() / af2 * ((gamma_eos[fluid] * c2 - s2) / den() + 2.0 * at2 * (d2 + s2) / d2 / den_prime
      + (gamma_eos[fluid] * c2 * (s2 - 2.0 * ax2) - s2 * s2) / d2 / den());
   double del_den = den() / (1.0 + q2) * (vel()[0] / af - 1.0);
   double af_star = vel()[0] - af * del_den / den();
   double cs_star = vel()[0] - c  * del_den / den();

   MHD<fluid>::PrimitiveState prim;
   
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
\date 01/15/2025
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC void MHD<fluid>::PrimitiveState::Invert(void)
{
   vel()[0] = -vel()[0];
   mag()[0] = -mag()[0];
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// MHD<fluid>::ConservedState methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] den_in Density
\param[in] mom_in Momentum
\param[in] enr_in Energy
\param[in] mag_in Magnetic field
\param[in] glm_in Largange multiplier
*/
template <int fluid>
#ifdef CL_MHD_GLM
SPECTRUM_DEVICE_FUNC MHD<fluid>::ConservedState::ConservedState(double den_in, const GeoVector& mom_in, double enr_in, const GeoVector& mag_in, double glm_in)
#else
SPECTRUM_DEVICE_FUNC MHD<fluid>::ConservedState::ConservedState(double den_in, const GeoVector& mom_in, double enr_in, const GeoVector& mag_in)
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
\date 01/15/2025
\param[in] other Object to initialize from
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC MHD<fluid>::ConservedState::ConservedState(const SimpleArray<double, CL_MHD_NVARS>& other)
                                               : SimpleArray<double, CL_MHD_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\return Fastest wave speed (max fast with GLM safety)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC double MHD<fluid>::ConservedState::FastestWave(void) const
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
\date 01/15/2025
\return Fastest normal wave speed (fast_n)
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC double MHD<fluid>::ConservedState::FastestWaveNormal(void) const
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
\date 01/15/2025
\return Vector of primitive variables
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC MHD<fluid>::PrimitiveState MHD<fluid>::ConservedState::ToPrimitive(void) const
{
   MHD<fluid>::PrimitiveState prim;

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
\date 01/15/2025
\return Flux vector
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC MHD<fluid>::FluxFunction MHD<fluid>::ConservedState::ToFlux(void) const
{
   GeoVector vel = mom() / den();
   double pre = Pressure(den(), vel.Norm2(), mag().Norm2(), enr(), fluid);
   double pre_tot = pre + mag().Norm2() / M_8PI;
   MHD<fluid>::FluxFunction flux;

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
\date 01/15/2025
\return Gas pressure
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC double MHD<fluid>::ConservedState::GetPressure(void) const
{
   return Pressure(den(), mom().Norm2() / Sqr(den()), mag().Norm2(), enr(), fluid);
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] flux Flux vector
\param[in] S1   First wave speed
\param[in] S2   Second wave speed
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC void MHD<fluid>::ConservedState::FixMomentum(MHD<fluid>::FluxFunction& flux, double S1, double S2)
{
   mom()[1] = (flux.mom()[1] - mag()[0] * mag()[1] / M_4PI) / (S1 - S2);
   mom()[2] = (flux.mom()[2] - mag()[0] * mag()[2] / M_4PI) / (S1 - S2);
};

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] flux Flux vector
\param[in] S1   First wave speed
\param[in] S2   Second wave speed
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC void MHD<fluid>::ConservedState::FixEnergy(MHD<fluid>::FluxFunction& flux, double S1, double S2)
{
   double pre_total = S1 * mom()[0] + mag()[0] * mag()[0] / M_4PI - flux.mom()[0];
   enr() = (flux.enr() - (mag()[0] * S2 + (mag()[1] * mom()[1] + mag()[2] * mom()[2]) / den()) * mag()[0] + pre_total * S2) / (S1 - S2);
};

/*!
\author Vladimir Florinski
\date 01/15/2025
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC void MHD<fluid>::ConservedState::Invert(void)
{
   mom()[0] = -mom()[0];
   mag()[0] = -mag()[0];
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// FluxFunctionMHD methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 01/15/2025
\param[in] denf_in Density flux
\param[in] momf_in Momentum flux
\param[in] enrf_in Energy flux
\param[in] magf_in Magnetic field
\param[in] glmf_in Largange multiplier flux
*/
template <int fluid>
#ifdef CL_MHD_GLM
SPECTRUM_DEVICE_FUNC MHD<fluid>::FluxFunction::FluxFunction(double denf_in, const GeoVector& momf_in, double enrf_in, const GeoVector& magf_in, double glmf_in)
#else
SPECTRUM_DEVICE_FUNC MHD<fluid>::FluxFunction::FluxFunction(double denf_in, const GeoVector& momf_in, double enrf_in, const GeoVector& magf_in)
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
\date 01/15/2025
\param[in] other Object to initialize from
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC MHD<fluid>::FluxFunction::FluxFunction(const SimpleArray<double, CL_MHD_NVARS>& other)
                                             : SimpleArray<double, CL_MHD_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 01/15/2025
*/
template <int fluid>
SPECTRUM_DEVICE_FUNC void MHD<fluid>::FluxFunction::Invert(void)
{
   momf()[0] = -momf()[0];
   magf()[0] = -magf()[0];
};

};
