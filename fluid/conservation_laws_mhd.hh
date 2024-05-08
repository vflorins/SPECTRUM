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

//! Number of fluid variables (DEN,VELx3,PRE,MAGx3)
#define CL_MHD_NVARS 8

//! Number of GLM variables
#define CL_MHD_GLM 1
#if (CL_MHD_GLM != 0) && (CL_MHD_GLM != 1)
#error The number of GLM variables must be 0 or 1
#endif

//! Safety factor for GLM transport
#define CL_MHD_GLMSAFETY 1.05

//! Total number of variables
#define CL_MHD_TOTAL (CL_MHD_NVARS + CL_MHD_GLM + n_ind)

template <int fluid, int n_ind> struct PrimitiveStateMHD;
template <int fluid, int n_ind> struct ConservedStateMHD;
template <int fluid, int n_ind> struct FluxFunctionMHD;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// PrimitiveStateMHD class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

template <int fluid, int n_ind>
struct PrimitiveStateMHD : public SimpleArray<double, CL_MHD_TOTAL>
{
   using SimpleArray<double, CL_MHD_TOTAL>::data;

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

#if CL_MHD_GLM != 0

//! Alias for Largange multiplier (RO)
   const double& glm(void) const {return data[8];};

//! Alias for Largange multiplier (RW)
   double& glm(void) {return data[8];};

#endif

//! Return the "fluid" template parameter
   SPECTRUM_DEVICE_FUNC int Fluid(void) const {return fluid;};

//! Return the "n_ind" template parameter
   SPECTRUM_DEVICE_FUNC int Nind(void) const {return n_ind;};

//! Return the number of main variables
   SPECTRUM_DEVICE_FUNC int Nmain(void) const {return CL_MHD_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC PrimitiveStateMHD(void) = default;

#if CL_MHD_GLM == 0
//! Constructor from components
   SPECTRUM_DEVICE_FUNC PrimitiveStateMHD(double den_in, const GeoVector& vel_in, double pre_in, const GeoVector& mag_in, double* ind_in);
#else
//! Constructor from components
   SPECTRUM_DEVICE_FUNC PrimitiveStateMHD(double den_in, const GeoVector& vel_in, double pre_in, const GeoVector& mag_in, double glm_in, double* ind_in);
#endif

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC PrimitiveStateMHD(const SimpleArray<double, CL_MHD_TOTAL>& other);

//! Calculate the fastest wave speed
   SPECTRUM_DEVICE_FUNC double FastestWave(void) const;

//! Calculate the fastest normal wave speed
   SPECTRUM_DEVICE_FUNC double FastestWaveNormal(void) const;

//! Calculate primitive state
   SPECTRUM_DEVICE_FUNC ConservedStateMHD<fluid, n_ind> ToConserved(bool ind_ok) const;

//! Calculate flux function
   SPECTRUM_DEVICE_FUNC FluxFunctionMHD<fluid, n_ind> ToFlux(bool ind_ok) const;

//! Invert vector variables
   SPECTRUM_DEVICE_FUNC void Invert(void);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ConservedStateMHD class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

template <int fluid, int n_ind>
struct ConservedStateMHD : public SimpleArray<double, CL_MHD_TOTAL>
{
   using SimpleArray<double, CL_MHD_TOTAL>::data;

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

#if CL_MHD_GLM != 0

//! Alias for Largange multiplier (RO)
   const double& glm(void) const {return data[8];};

//! Alias for Largange multiplier (RW)
   double& glm(void) {return data[8];};

#endif

//! Return the "fluid" template parameter
   SPECTRUM_DEVICE_FUNC int Fluid(void) const {return fluid;};

//! Return the "n_ind" template parameter
   SPECTRUM_DEVICE_FUNC int Nind(void) const {return n_ind;};

//! Return the number of main variables
   SPECTRUM_DEVICE_FUNC int Nmain(void) const {return CL_MHD_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC ConservedStateMHD(void) = default;

#if CL_MHD_GLM == 0
//! Constructor from components
   SPECTRUM_DEVICE_FUNC ConservedStateMHD(double den_in, const GeoVector& mom_in, double enr_in, const GeoVector& mag_in, double* ind_in);
#else
//! Constructor from components
   SPECTRUM_DEVICE_FUNC ConservedStateMHD(double den_in, const GeoVector& mom_in, double enr_in, const GeoVector& mag_in, double glm_in, double* ind_in);
#endif

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC ConservedStateMHD(const SimpleArray<double, CL_MHD_TOTAL>& other);

//! Calculate the fastest wave speed
   SPECTRUM_DEVICE_FUNC double FastestWave(void) const;

//! Calculate the fastest normal wave speed
   SPECTRUM_DEVICE_FUNC double FastestWaveNormal(void) const;

//! Calculate primitive state
   SPECTRUM_DEVICE_FUNC PrimitiveStateMHD<fluid, n_ind> ToPrimitive(bool ind_ok) const;

//! Calculate flux function
   SPECTRUM_DEVICE_FUNC FluxFunctionMHD<fluid, n_ind> ToFlux(bool ind_ok) const;

//! Calculate _gas_ pressure
   SPECTRUM_DEVICE_FUNC double GetPressure(void) const;

//! Calculate y and z components of momentum from the externally provided flux in a moving frame
   SPECTRUM_DEVICE_FUNC void FixMomentum(FluxFunctionMHD<fluid, n_ind>& flux, double S1, double S2);

//! Calculate energy from the externally provided flux in a moving frame
   SPECTRUM_DEVICE_FUNC void FixEnergy(FluxFunctionMHD<fluid, n_ind>& flux, double S1, double S2);

//! Invert vector variables
   SPECTRUM_DEVICE_FUNC void Invert(void);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// FluxFunctionMHD class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

template <int fluid, int n_ind>
struct FluxFunctionMHD : public SimpleArray<double, CL_MHD_TOTAL>
{
   using SimpleArray<double, CL_MHD_TOTAL>::data;

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

#if CL_MHD_GLM != 0

//! Alias for Largange multiplier flux (RO)
   const double& glmf(void) const {return data[8];};

//! Alias for Largange multiplier flux (RO)
   double& glmf(void) {return data[8];};

#endif

//! Return the "fluid" template parameter
   SPECTRUM_DEVICE_FUNC int Fluid(void) const {return fluid;};

//! Return the "n_ind" template parameter
   SPECTRUM_DEVICE_FUNC int Nind(void) const {return n_ind;};

//! Return the number of main variables
   SPECTRUM_DEVICE_FUNC int Nmain(void) const {return CL_MHD_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC FluxFunctionMHD(void) = default;

#if CL_MHD_GLM == 0
//! Constructor from components
   SPECTRUM_DEVICE_FUNC FluxFunctionMHD(double denf_in, const GeoVector& momf_in, double enrf_in, const GeoVector& magf_in, double* indf_in);
#else
//! Constructor from components
   SPECTRUM_DEVICE_FUNC FluxFunctionMHD(double denf_in, const GeoVector& momf_in, double enrf_in, const GeoVector& magf_in, double glmf_in, double* indf_in);
#endif

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC FluxFunctionMHD(const SimpleArray<double, CL_MHD_TOTAL>& other);

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
\param[in] ind_in Indicator variables
*/
template <int fluid, int n_ind>
#if CL_MHD_GLM == 0
SPECTRUM_DEVICE_FUNC inline PrimitiveStateMHD<fluid, n_ind>::PrimitiveStateMHD(double den_in, const GeoVector& vel_in, double pre_in,
                                                                               const GeoVector& mag_in, double* ind_in)
#else
SPECTRUM_DEVICE_FUNC inline PrimitiveStateMHD<fluid, n_ind>::PrimitiveStateMHD(double den_in, const GeoVector& vel_in, double pre_in,
                                                                               const GeoVector& mag_in, double glm_in, double* ind_in)
#endif
{
   den() = den_in;
   vel() = vel_in;
   pre() = pre_in;
   mag() = mag_in;

#if CL_MHD_GLM != 0
   glm() = glm_in;
#endif

   if(n_ind != 0) memcpy(data + CL_MHD_NVARS + CL_MHD_GLM, ind_in, n_ind * sizeof(double));
};

/*!
\author Vladimir Florinski
\date 03/14/2024
\param[in] other Object to initialize from
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline PrimitiveStateMHD<fluid, n_ind>::PrimitiveStateMHD(const SimpleArray<double, CL_MHD_TOTAL>& other)
   : SimpleArray<double, CL_MHD_TOTAL>(other)
{
};

/*!
\author Vladimir Florinski
\date 03/16/2024
\return Fastest wave speed (max fast with GLM safety)
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline double PrimitiveStateMHD<fluid, n_ind>::FastestWave(void) const
{
   double Va2 = Alfven2(den(), mag().Norm2());
   double Cs2 = Sound2(den(), pre(), fluid);

// FIXME decide where to use "CL_MHD_GLMSAFETY"
//#if CL_MHD_GLM == 0
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
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline double PrimitiveStateMHD<fluid, n_ind>::FastestWaveNormal(void) const
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
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline ConservedStateMHD<fluid, n_ind> PrimitiveStateMHD<fluid, n_ind>::ToConserved(bool ind_ok) const
{
   ConservedStateMHD<fluid, n_ind> cons;

   cons.den() = den();
   cons.mom() = vel() * den();
   cons.enr() = Energy(den(), vel().Norm2(), mag().Norm2(), pre(), fluid);
   cons.mag() = mag();

#if CL_MHD_GLM != 0
   cons.glm() = glm();
#endif

// Optionally convert the indicator variables
   if(ind_ok) {
      for(auto i = CL_MHD_NVARS + CL_MHD_GLM; i < CL_MHD_NVARS + CL_MHD_GLM + n_ind; i++) {
         cons.data[i] = den() * data[i];
      };
   };
   return cons;
};

/*!
\author Vladimir Florinski
\date 03/06/2024
\return Flux vector
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline FluxFunctionMHD<fluid, n_ind> PrimitiveStateMHD<fluid, n_ind>::ToFlux(bool ind_ok) const
{
   double enr = Energy(den(), vel().Norm2(), mag().Norm2(), pre(), fluid);
   double pre_tot = pre() + mag().Norm2() / eightpi;
   FluxFunctionMHD<fluid, n_ind> flux;

   flux.denf() = den() * vel()[0];
   flux.momf() = den() * vel()[0] * vel() + pre_tot * gv_nx - (mag()[0] / fourpi) * mag();
   flux.enrf() = (enr + pre_tot) * vel()[0] - (vel() * mag()) * mag()[0] / fourpi;
   flux.magf() = vel()[0] * mag() - mag()[0] * vel();

#if CL_MHD_GLM != 0
// The GLM flux is proportional to the square of the fastest speed. The safety factor "CL_MHD_GLMSAFETY" is used to make sure it is slightly outside of the "normal" Riemann fan.
   flux.magf() += glm() * gv_nx;

// FIXME should this flux be calculated in the RS?
//   flux.glmf() = Sqr(CL_MHD_GLMSAFETY * FastestWaveNormal()) * mag()[0];
#endif

// Optionally convert the indicator variables
   if(ind_ok) {
      for(auto i = CL_MHD_NVARS + CL_MHD_GLM; i < CL_MHD_NVARS + CL_MHD_GLM + n_ind; i++) {
         flux.data[i] = den() * data[i] * vel()[0];
      };
   };
   return flux;
};

/*!
\author Vladimir Florinski
\date 04/09/2024
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline void PrimitiveStateMHD<fluid, n_ind>::Invert(void)
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
\param[in] ind_in Indicator variables
*/
template <int fluid, int n_ind>
#if CL_MHD_GLM == 0
SPECTRUM_DEVICE_FUNC inline ConservedStateMHD<fluid, n_ind>::ConservedStateMHD(double den_in, const GeoVector& mom_in, double enr_in,
                                                                               const GeoVector& mag_in, double* ind_in)
#else
SPECTRUM_DEVICE_FUNC inline ConservedStateMHD<fluid, n_ind>::ConservedStateMHD(double den_in, const GeoVector& mom_in, double enr_in,
                                                                               const GeoVector& mag_in, double glm_in, double* ind_in)
#endif
{
   den() = den_in;
   mom() = mom_in;
   enr() = enr_in;
   mag() = mag_in;

#if CL_MHD_GLM != 0
   glm() = glm_in;
#endif

   if(n_ind != 0) memcpy(data + CL_MHD_NVARS + CL_MHD_GLM, ind_in, n_ind * sizeof(double));
};

/*!
\author Vladimir Florinski
\date 03/14/2024
\param[in] other Object to initialize from
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline ConservedStateMHD<fluid, n_ind>::ConservedStateMHD(const SimpleArray<double, CL_MHD_TOTAL>& other)
   : SimpleArray<double, CL_MHD_TOTAL>(other)
{
};

/*!
\author Vladimir Florinski
\date 03/16/2024
\return Fastest wave speed (max fast with GLM safety)
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline double ConservedStateMHD<fluid, n_ind>::FastestWave(void) const
{
   GeoVector vel = mom() / den();
   double pre = Pressure(den(), vel.Norm2(), mag().Norm2(), enr(), fluid);
   double Va2 = Alfven2(den(), mag().Norm2());
   double Cs2 = Sound2(den(), pre, fluid);

// FIXME decide where to use "CL_MHD_GLMSAFETY"
//#if CL_MHD_GLM == 0
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
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline double ConservedStateMHD<fluid, n_ind>::FastestWaveNormal(void) const
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
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline PrimitiveStateMHD<fluid, n_ind> ConservedStateMHD<fluid, n_ind>::ToPrimitive(bool ind_ok) const
{
   PrimitiveStateMHD<fluid, n_ind> prim;

   prim.den() = den();
   prim.vel() = mom() / den();
   prim.pre() = Pressure(den(), prim.vel().Norm2(), mag().Norm2(), enr(), fluid);
   prim.mag() = mag();

#if CL_MHD_LAGR != 0
   prim.glm() = glm();
#endif

// Optionally convert the indicator variables
   if(ind_ok) {
      for(auto i = CL_MHD_NVARS + CL_MHD_GLM; i < CL_MHD_NVARS + CL_MHD_GLM + n_ind; i++) {
         prim.data[i] = data[i] / den();
      };
   };
   return prim;
};

/*!
\author Vladimir Florinski
\date 03/14/2024
\return Flux vector
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline FluxFunctionMHD<fluid, n_ind> ConservedStateMHD<fluid, n_ind>::ToFlux(bool ind_ok) const
{
   GeoVector vel = mom() / den();
   double pre = Pressure(den(), vel.Norm2(), mag().Norm2(), enr(), fluid);
   double pre_tot = pre + mag().Norm2() / eightpi;
   FluxFunctionMHD<fluid, n_ind> flux;

   flux.denf() = mom()[0];
   flux.momf() = mom()[0] * vel + pre_tot * gv_nx - (mag()[0] / fourpi) * mag();
   flux.enrf() = (enr() + pre_tot) * vel[0] - (vel * mag()) * mag()[0] / fourpi;
   flux.magf() = vel[0] * mag() - mag()[0] * vel;

#if CL_MHD_GLM != 0
// The GLM flux is proportional to the square of the fastest speed. The safety factor "CL_MHD_GLMSAFETY" is used to make sure it is slightly outside of the "normal" Riemann fan.
   flux.magf() += glm() * gv_nx;

// FIXME should this flux be calculated in the RS?
//   flux.glmf() = Sqr(CL_MHD_GLMSAFETY * FastestWaveNormal()) * mag()[0];
#endif

// Optionally convert the indicator variables
   if(ind_ok) {
      for(auto i = CL_MHD_NVARS; i < CL_MHD_NVARS + n_ind; i++) {
         flux.data[i] = data[i] * mom()[0];
      };
   };
   return flux;
};

/*!
\author Vladimir Florinski
\date 03/20/2024
\return Gas pressure
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline double ConservedStateMHD<fluid, n_ind>::GetPressure(void) const
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
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline void ConservedStateMHD<fluid, n_ind>::FixMomentum(FluxFunctionMHD<fluid, n_ind>& flux, double S1, double S2)
{
   mom()[1] = (flux.mom()[1] - mag()[0] * mag()[1] / fourpi) / (S1 - S2);
   mom()[2] = (flux.mom()[2] - mag()[0] * mag()[2] / fourpi) / (S1 - S2);
};

/*!
\author Vladimir Florinski
\date 03/20/2024
\param[in] flux Flux vector
\param[in] S1   First wave speed
\param[in] S2   Second wave speed
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline void ConservedStateMHD<fluid, n_ind>::FixEnergy(FluxFunctionMHD<fluid, n_ind>& flux, double S1, double S2)
{
   double pre_total = S1 * mom()[0] + mag()[0] * mag()[0] / fourpi - flux.mom()[0];
   enr() = (flux.enr() - (mag()[0] * S2 + (mag()[1] * mom()[1] + mag()[2] * mom()[2]) / den()) * mag()[0] + pre_total * S2) / (S1 - S2);
};

/*!
\author Vladimir Florinski
\date 04/09/2024
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline void ConservedStateMHD<fluid, n_ind>::Invert(void)
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
\param[in] indf_in Indicator variables flux
*/
template <int fluid, int n_ind>
#if CL_MHD_GLM == 0
SPECTRUM_DEVICE_FUNC inline FluxFunctionMHD<fluid, n_ind>::FluxFunctionMHD(double denf_in, const GeoVector& momf_in, double enrf_in,
                                                                           const GeoVector& magf_in, double* indf_in)
#else
SPECTRUM_DEVICE_FUNC inline FluxFunctionMHD<fluid, n_ind>::FluxFunctionMHD(double denf_in, const GeoVector& momf_in, double enrf_in,
                                                                           const GeoVector& magf_in, double glmf_in, double* indf_in)
#endif
{
   denf() = denf_in;
   momf() = momf_in;
   enrf() = enrf_in;
   magf() = magf_in;

#if CL_MHD_GLM != 0
   glmf() = glmf_in;
#endif

   if(n_ind != 0) memcpy(data + CL_MHD_NVARS + CL_MHD_GLM, indf_in, n_ind * sizeof(double));
};

/*!
\author Vladimir Florinski
\date 03/14/2024
\param[in] other Object to initialize from
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline FluxFunctionMHD<fluid, n_ind>::FluxFunctionMHD(const SimpleArray<double, CL_MHD_TOTAL>& other)
   : SimpleArray<double, CL_MHD_TOTAL>(other)
{
};

/*!
\author Vladimir Florinski
\date 04/09/2024
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline void FluxFunctionMHD<fluid, n_ind>::Invert(void)
{
   momf()[0] = -momf()[0];
   magf()[0] = -magf()[0];
};

};

#endif
