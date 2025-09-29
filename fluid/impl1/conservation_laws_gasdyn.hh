/*!
\file conservation_laws_gasdyn.hh
\brief Declares rules for compressible gas-dynamic equations
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_CONSERVATION_LAWS_GASDYN_HH
#define SPECTRUM_CONSERVATION_LAWS_GASDYN_HH

#include "common/vectors.hh"

namespace Spectrum {

//! Number of fluid variables (DEN,VELx3,PRE)
#define CL_GASDYN_NVARS 5

/*!
\brief Codifies physics of compressible gas dynamics
\author Vladimir Florinski
*/
template <int fluid>
class Gasdyn
{
public:

   struct ConservedState;
   struct FluxFunction;

//! A trait to be used in template specializations
   static constexpr bool is_cons_law_active = true;

//! Return the "fluid" template parameter
   SPECTRUM_DEVICE_FUNC static constexpr int Fluid(void) {return fluid;};

//! Return the number of variables
   SPECTRUM_DEVICE_FUNC static constexpr int Nvars(void) {return CL_GASDYN_NVARS;};

/*!
\brief Codifies physics of compressible gas dynamics for primitive variables
\author Vladimir Florinski
*/
   struct PrimitiveState : public SimpleArray<double, CL_GASDYN_NVARS>
   {
      using SimpleArray<double, CL_GASDYN_NVARS>::data;

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

   //! Default constructor
      SPECTRUM_DEVICE_FUNC PrimitiveState(void) = default;

   //! Constructor from components
      SPECTRUM_DEVICE_FUNC PrimitiveState(double den_in, const GeoVector& vel_in, double pre_in);

   //! Constructor from SimpleArray
      SPECTRUM_DEVICE_FUNC PrimitiveState(const SimpleArray<double, CL_GASDYN_NVARS>& other);

   //! Calculate the fastest wave speed
      SPECTRUM_DEVICE_FUNC double FastestWave(void) const;

   //! Calculate the fastest normal wave speed
      SPECTRUM_DEVICE_FUNC double FastestWaveNormal(void) const;

   //! Calculate conserved state
      SPECTRUM_DEVICE_FUNC ConservedState ToConserved(void) const;

   //! Calculate flux function
      SPECTRUM_DEVICE_FUNC FluxFunction ToFlux(void) const;

   //! Calculate sonic state
      SPECTRUM_DEVICE_FUNC PrimitiveState ToSonic(void) const;

   //! Invert vector variables about the x-axis
      SPECTRUM_DEVICE_FUNC void Invert(void);
   };

/*!
\brief Codifies physics of compressible gas dynamics for conserved variables
\author Vladimir Florinski
*/
   struct ConservedState : public SimpleArray<double, CL_GASDYN_NVARS>
   {
      using SimpleArray<double, CL_GASDYN_NVARS>::data;

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

   //! Default constructor
      SPECTRUM_DEVICE_FUNC ConservedState(void) = default;

   //! Constructor from components
      SPECTRUM_DEVICE_FUNC ConservedState(double den_in, const GeoVector& mom_in, double enr_in);

   //! Constructor from SimpleArray
      SPECTRUM_DEVICE_FUNC ConservedState(const SimpleArray<double, CL_GASDYN_NVARS>& other);

   //! Calculate the fastest wave speed
      SPECTRUM_DEVICE_FUNC double FastestWave(void) const;

   //! Calculate the fastest normal wave speed
      SPECTRUM_DEVICE_FUNC double FastestWaveNormal(void) const;

   //! Calculate primitive state
      SPECTRUM_DEVICE_FUNC PrimitiveState ToPrimitive(void) const;

   //! Calculate flux function
      SPECTRUM_DEVICE_FUNC FluxFunction ToFlux(void) const;

   //! Calculate _gas_ pressure
      SPECTRUM_DEVICE_FUNC double GetPressure(void) const;

   //! Calculate y and z components of momentum from the externally provided flux in a moving frame
      SPECTRUM_DEVICE_FUNC void FixMomentum(FluxFunction& flux, double S1, double S2);

   //! Calculate energy from the externally provided flux in a moving frame
      SPECTRUM_DEVICE_FUNC void FixEnergy(FluxFunction& flux, double S1, double S2);

   //! Invert vector variables about the x-axis
      SPECTRUM_DEVICE_FUNC void Invert(void);
   };

/*!
\brief Codifies physics of compressible gas dynamics for the fluxes
\author Vladimir Florinski
*/
   struct FluxFunction : public SimpleArray<double, CL_GASDYN_NVARS>
   {
      using SimpleArray<double, CL_GASDYN_NVARS>::data;

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

   //! Default constructor
      SPECTRUM_DEVICE_FUNC FluxFunction(void) = default;

   //! Constructor from components
      SPECTRUM_DEVICE_FUNC FluxFunction(double denf_in, const GeoVector& momf_in, double enrf_in);

   //! Constructor from SimpleArray
      SPECTRUM_DEVICE_FUNC FluxFunction(const SimpleArray<double, CL_GASDYN_NVARS>& other);

   //! Invert vector variables about the x-axis
      SPECTRUM_DEVICE_FUNC void Invert(void);
   };

};

};

#endif
