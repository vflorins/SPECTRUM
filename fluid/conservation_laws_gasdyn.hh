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

//! Total number of variables
#define CL_GASDYN_TOTAL (CL_GASDYN_NVARS + n_ind)

template <int fluid, int n_ind> struct PrimitiveStateGasdyn;
template <int fluid, int n_ind> struct ConservedStateGasdyn;
template <int fluid, int n_ind> struct FluxFunctionGasdyn;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// PrimitiveStateGasdyn class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

template <int fluid, int n_ind>
struct PrimitiveStateGasdyn : public SimpleArray<double, CL_GASDYN_TOTAL>
{
   using SimpleArray<double, CL_GASDYN_TOTAL>::data;

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
   SPECTRUM_DEVICE_FUNC int Fluid(void) const {return fluid;};

//! Return the "n_ind" template parameter
   SPECTRUM_DEVICE_FUNC int Nind(void) const {return n_ind;};

//! Return the number of main variables
   SPECTRUM_DEVICE_FUNC int Nmain(void) const {return CL_GASDYN_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC PrimitiveStateGasdyn(void) = default;

//! Constructor from components
   SPECTRUM_DEVICE_FUNC PrimitiveStateGasdyn(double den_in, const GeoVector& vel_in, double pre_in, double* ind_in);

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC PrimitiveStateGasdyn(const SimpleArray<double, CL_GASDYN_TOTAL>& other);

//! Calculate the fastest wave speed
   SPECTRUM_DEVICE_FUNC double FastestWave(void) const;

//! Calculate the fastest normal wave speed
   SPECTRUM_DEVICE_FUNC double FastestWaveNormal(void) const;

//! Calculate conserved state
   SPECTRUM_DEVICE_FUNC ConservedStateGasdyn<fluid, n_ind> ToConserved(bool ind_ok) const;

//! Calculate flux function
   SPECTRUM_DEVICE_FUNC FluxFunctionGasdyn<fluid, n_ind> ToFlux(bool ind_ok) const;

//! Calculate sonic state
   SPECTRUM_DEVICE_FUNC PrimitiveStateGasdyn<fluid, n_ind> ToSonic(bool ind_ok) const;

//! Invert vector variables
   SPECTRUM_DEVICE_FUNC void Invert(void);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ConservedStateGasdyn class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

template <int fluid, int n_ind>
struct ConservedStateGasdyn : public SimpleArray<double, CL_GASDYN_TOTAL>
{
   using SimpleArray<double, CL_GASDYN_TOTAL>::data;

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
   SPECTRUM_DEVICE_FUNC int Fluid(void) const {return fluid;};

//! Return the "n_ind" template parameter
   SPECTRUM_DEVICE_FUNC int Nind(void) const {return n_ind;};

//! Return the number of main variables
   SPECTRUM_DEVICE_FUNC int Nmain(void) const {return CL_GASDYN_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC ConservedStateGasdyn(void) = default;

//! Constructor from components
   SPECTRUM_DEVICE_FUNC ConservedStateGasdyn(double den_in, const GeoVector& mom_in, double enr_in, double* ind_in);

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC ConservedStateGasdyn(const SimpleArray<double, CL_GASDYN_TOTAL>& other);

//! Calculate the fastest wave speed
   SPECTRUM_DEVICE_FUNC double FastestWave(void) const;

//! Calculate the fastest normal wave speed
   SPECTRUM_DEVICE_FUNC double FastestWaveNormal(void) const;

//! Calculate primitive state
   SPECTRUM_DEVICE_FUNC PrimitiveStateGasdyn<fluid, n_ind> ToPrimitive(bool ind_ok) const;

//! Calculate flux function
   SPECTRUM_DEVICE_FUNC FluxFunctionGasdyn<fluid, n_ind> ToFlux(bool ind_ok) const;

//! Calculate _gas_ pressure
   SPECTRUM_DEVICE_FUNC double GetPressure(void) const;

//! Calculate y and z components of momentum from the externally provided flux in a moving frame
   SPECTRUM_DEVICE_FUNC void FixMomentum(FluxFunctionGasdyn<fluid, n_ind>& flux, double S1, double S2);

//! Calculate energy from the externally provided flux in a moving frame
   SPECTRUM_DEVICE_FUNC void FixEnergy(FluxFunctionGasdyn<fluid, n_ind>& flux, double S1, double S2);

//! Invert vector variables
   SPECTRUM_DEVICE_FUNC void Invert(void);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// FluxFunctionGasdyn class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

template <int fluid, int n_ind>
struct FluxFunctionGasdyn : public SimpleArray<double, CL_GASDYN_TOTAL>
{
   using SimpleArray<double, CL_GASDYN_TOTAL>::data;

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
   SPECTRUM_DEVICE_FUNC int Fluid(void) const {return fluid;};

//! Return the "n_ind" template parameter
   SPECTRUM_DEVICE_FUNC int Nind(void) const {return n_ind;};

//! Return the number of main variables
   SPECTRUM_DEVICE_FUNC int Nmain(void) const {return CL_GASDYN_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC FluxFunctionGasdyn(void) = default;

//! Constructor from components
   SPECTRUM_DEVICE_FUNC FluxFunctionGasdyn(double denf_in, const GeoVector& momf_in, double enrf_in, double* indf_in);

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC FluxFunctionGasdyn(const SimpleArray<double, CL_GASDYN_TOTAL>& other);

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
\param[in] ind_in Indicator variables
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline PrimitiveStateGasdyn<fluid, n_ind>::PrimitiveStateGasdyn(double den_in, const GeoVector& vel_in,
                                                                                     double pre_in, double* ind_in)
{
   den() = den_in;
   vel() = vel_in;
   pre() = pre_in;
   if(n_ind != 0) memcpy(data + CL_GASDYN_NVARS, ind_in, n_ind * sizeof(double));
};

/*!
\author Vladimir Florinski
\date 03/14/2024
\param[in] other Object to initialize from
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline PrimitiveStateGasdyn<fluid, n_ind>::PrimitiveStateGasdyn(const SimpleArray<double, CL_GASDYN_TOTAL>& other)
   : SimpleArray<double, CL_GASDYN_TOTAL>(other)
{
};

/*!
\author Vladimir Florinski
\date 03/16/2024
\return Fastest wave speed (sound)
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline double PrimitiveStateGasdyn<fluid, n_ind>::FastestWave(void) const
{
   return SoundSpeed(den(), pre(), fluid);
};

/*!
\author Vladimir Florinski
\date 03/16/2024
\return Fastest normal wave speed (sound)
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline double PrimitiveStateGasdyn<fluid, n_ind>::FastestWaveNormal(void) const
{
   return SoundSpeed(den(), pre(), fluid);
};

/*!
\author Vladimir Florinski
\date 03/06/2024
\return Vector of conserved variables
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline ConservedStateGasdyn<fluid, n_ind> PrimitiveStateGasdyn<fluid, n_ind>::ToConserved(bool ind_ok) const
{
   ConservedStateGasdyn<fluid, n_ind> cons;

   cons.den() = den();
   cons.mom() = den() * vel();
   cons.enr() = Energy(den(), vel().Norm2(), 0.0, pre(), fluid);

// Optionally convert the indicator variables
   if(ind_ok) {
      for(auto i = CL_GASDYN_NVARS; i < CL_GASDYN_NVARS + n_ind; i++) {
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
SPECTRUM_DEVICE_FUNC inline FluxFunctionGasdyn<fluid, n_ind> PrimitiveStateGasdyn<fluid, n_ind>::ToFlux(bool ind_ok) const
{
   double enr = Energy(den(), vel().Norm2(), 0.0, pre(), fluid);
   FluxFunctionGasdyn<fluid, n_ind> flux;

   flux.denf() = den() * vel()[0];
   flux.momf() = den() * vel()[0] * vel() + pre() * gv_nx;
   flux.enrf() = (enr + pre()) * vel()[0];

// Optionally convert the indicator variables
   if(ind_ok) {
      for(auto i = CL_GASDYN_NVARS; i < CL_GASDYN_NVARS + n_ind; i++) {
         flux.data[i] = den() * data[i] * vel()[0];
      };
   };
   return flux;
};

/*!
\author Vladimir Florinski
\date 04/10/2024
\return Vector of primitive variables in a sonic state
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC PrimitiveStateGasdyn<fluid, n_ind> PrimitiveStateGasdyn<fluid, n_ind>::ToSonic(bool ind_ok) const
{
   double cs = Sound(den(), pre(), fluid);
   double gamm1 = gamma_eos[fluid] - 1.0;
   double gamp1 = gamma_eos[fluid] + 1.0;
   double del_den = 2.8 * den() / gamp1 * (vel()[0] / cs - 1.0);
   double cs_star = cs * (1.0 + 0.5 * gamm1 * cs * del_den / den());
   PrimitiveStateGasdyn<fluid, n_ind> prim;
   
   prim.den() = den() + del_den;
   prim.vel()[0] = vel()[0] - cs * del_den / den();
   prim.vel()[1] = vel()[1];
   prim.vel()[2] = vel()[2];
   prim.pre() = prim.den() * Sqr(cs_star) / gamma_eos[fluid];
};

/*!
\author Vladimir Florinski
\date 04/09/2024
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline void PrimitiveStateGasdyn<fluid, n_ind>::Invert(void)
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
\param[in] ind_in Indicator variables
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline ConservedStateGasdyn<fluid, n_ind>::ConservedStateGasdyn(double den_in, const GeoVector& mom_in,
                                                                                     double enr_in, double* ind_in)
{
   den() = den_in;
   mom() = mom_in;
   enr() = enr_in;
   if(n_ind != 0) memcpy(data + CL_GASDYN_NVARS, ind_in, n_ind * sizeof(double));
};

/*!
\author Vladimir Florinski
\date 03/14/2024
\param[in] other Object to initialize from
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline ConservedStateGasdyn<fluid, n_ind>::ConservedStateGasdyn(const SimpleArray<double, CL_GASDYN_TOTAL>& other)
   : SimpleArray<double, CL_GASDYN_TOTAL>(other)
{
};

/*!
\author Vladimir Florinski
\date 03/16/2024
\return Fastest wave speed (sound)
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline double ConservedStateGasdyn<fluid, n_ind>::FastestWave(void) const
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
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline double ConservedStateGasdyn<fluid, n_ind>::FastestWaveNormal(void) const
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
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline PrimitiveStateGasdyn<fluid, n_ind> ConservedStateGasdyn<fluid, n_ind>::ToPrimitive(bool ind_ok) const
{
   PrimitiveStateGasdyn<fluid, n_ind> prim;

   prim.den() = den();
   prim.vel() = mom() / den();
   prim.pre() = Pressure(den(), prim.vel().Norm2(), 0.0, enr(), fluid);

// Optionally convert the indicator variables
   if(ind_ok) {
      for(auto i = CL_GASDYN_NVARS; i < CL_GASDYN_NVARS + n_ind; i++) {
         prim.data[i] = data[i] / den();
      };
   };
   return prim;
};

/*!
\author Vladimir Florinski
\date 03/06/2024
\return Flux vector
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline FluxFunctionGasdyn<fluid, n_ind> ConservedStateGasdyn<fluid, n_ind>::ToFlux(bool ind_ok) const
{
   double pre = Pressure(den(), mom().Norm2() / Sqr(den()), 0.0, enr(), fluid);

   FluxFunctionGasdyn<fluid, n_ind> flux;
   flux.denf() = mom()[0];
   flux.momf() = mom()[0] * mom() / den() + pre * gv_nx;
   flux.enrf() = (enr() + pre) * mom()[0] / den();

// Optionally convert the indicator variables
   if(ind_ok) {
      for(auto i = CL_GASDYN_NVARS; i < CL_GASDYN_NVARS + n_ind; i++) {
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
SPECTRUM_DEVICE_FUNC inline double ConservedStateGasdyn<fluid, n_ind>::GetPressure(void) const
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
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline void ConservedStateGasdyn<fluid, n_ind>::FixMomentum(FluxFunctionGasdyn<fluid, n_ind>& flux, double S1, double S2)
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
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline void ConservedStateGasdyn<fluid, n_ind>::FixEnergy(FluxFunctionGasdyn<fluid, n_ind>& flux, double S1, double S2)
{
   double pre_total = S1 * mom()[0] - flux.mom()[0];
   enr() = (flux.enr() + pre_total * S2) / (S1 - S2);
};

/*!
\author Vladimir Florinski
\date 04/09/2024
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline void ConservedStateGasdyn<fluid, n_ind>::Invert(void)
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
\param[in] indf_in Indicator variables flux
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline FluxFunctionGasdyn<fluid, n_ind>::FluxFunctionGasdyn(double denf_in, const GeoVector& momf_in,
                                                                                 double enrf_in, double* indf_in)
{
   denf() = denf_in;
   momf() = momf_in;
   enrf() = enrf_in;
   if(n_ind != 0) memcpy(data + CL_GASDYN_NVARS, indf_in, n_ind * sizeof(double));
};

/*!
\author Vladimir Florinski
\date 03/14/2024
\param[in] other Object to initialize from
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline FluxFunctionGasdyn<fluid, n_ind>::FluxFunctionGasdyn(const SimpleArray<double, CL_GASDYN_TOTAL>& other)
   : SimpleArray<double, CL_GASDYN_TOTAL>(other)
{
};

/*!
\author Vladimir Florinski
\date 04/09/2024
*/
template <int fluid, int n_ind>
SPECTRUM_DEVICE_FUNC inline void FluxFunctionGasdyn<fluid, n_ind>::Invert(void)
{
   momf()[0] = -momf()[0];
};

};

#endif
