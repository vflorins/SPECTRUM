/*!
\file diffusion_base.hh
\brief Declares a base class to compute difusion coefficients
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DIFFUSION_BASE_HH
#define SPECTRUM_DIFFUSION_BASE_HH

// This includes (algorithm, cmath, cstdint, cstring, exception, fstream, vector), data_container, definitions, multi_index, vectors
#include "config.h"
#include "common/params.hh"
#include "common/physics.hh"
#include "common/derivativedata.hh"
#include <memory>

namespace Spectrum {

using namespace DiffusionOptions;

//! The diffusion model is independent of the background
constexpr uint16_t DIFF_NOBACKGROUND = 0x0010;

//! Clone function pattern
#define CloneFunctionDiffusion(T) std::unique_ptr<DiffusionBase> Clone(void) const override {return std::make_unique<T>();};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionBase class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Base class to calculate diffusion coefficients
\author Vladimir Florinski

The "DiffusionXXXXX" classes calculate diffusion coefficients in a broad sense. These could be pitch-angle scattering, spatial diffusion, gyrophase diffusion, momentum diffusion, etc. This class is required to be fast, so many of the checks performed by larger classes are not done. This class uses the momentum vector in the form (p,mu,phi).

Parameters:
*/
template <typename HConfig_>
class DiffusionBase : public Params {
public:

   using HConfig = HConfig_;
   using DiffusionCoordinates = HConfig::DiffusionConfig::Coordinates;
   using DiffusionFields = HConfig::DiffusionConfig::Fields;
   using HConfig::specie;

protected:

//! Coordinates (transient)
   DiffusionCoordinates _coords;

//! Fields data (transient)
   DiffusionFields _fields;

//! Square of the pitch angle sine (transient)
   double st2;

//! Cyclotron frequency (transient)
   double Omega;

//! Diffusion coefficient packed as a vector (transient)
   GeoVector Kappa;

//! Default constructor (protected, class not designed to be instantiated)
   DiffusionBase(void);

//! Constructor with arguments (to speed up construction of derived classes)
   DiffusionBase(const std::string& name_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   DiffusionBase(const DiffusionBase& other);

//! Set up the diffusion model based on "params"
   virtual void SetupDiffusion(bool construct);

//! Compute the diffusion coefficients
   virtual void EvaluateDiffusion(Component comp);

public:

//! Destructor
   virtual ~DiffusionBase() = default;

//! Clone function (stub)
   virtual std::unique_ptr<DiffusionBase> Clone(void) const = 0;

//! Set up the class parameters
   void SetupObject(const DataContainer& cont_in);

//! Stage at target coordinates for any kind of diffusion evaluation
   void Stage(const DiffusionCoordinates& coords, const DiffusionFields& fields);

//! Evaluate and return one diffusion component
   double Get(Component comp);

//! Compute derivative of diffusion coefficient in position or time. By default, it is computed numerically, but specific classes can override with analytic expressions.
   virtual double GetDirectionalDerivative(Component comp, int xyz, const DerivativeData& ddata);

//! Compute derivative of diffusion coefficient in mu. By default, it is computed numerically, but specific classes can override with analytic expressions.
   virtual double GetMuDerivative(Component comp);

};

};

// Something like this is needed for templated classes
#include "diffusion_base.cc"

#endif
