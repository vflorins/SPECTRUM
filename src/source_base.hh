/*!
\file source_base.hh
\brief Declares a base class to compute source terms
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SOURCE_BASE_HH
#define SPECTRUM_SOURCE_BASE_HH

// This includes (algorithm, cmath, cstdint, cstring, exception, fstream, vector), data_container, definitions, multi_index, vectors
#include "config.h"
#include "common/params.hh"
#include "common/physics.hh"
#include "common/spatial_data.hh"
#include <memory>

namespace Spectrum {

//! The source model is independent of the background
const uint16_t SOURCE_NOBACKGROUND = 0x0010;

//! Clone function pattern
#define CloneFunctionSource(T) std::unique_ptr<SourceBase> Clone(void) const override {return std::make_unique<T>();};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SourceBase class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Base class to calculate source terms
\author Juan G Alonso Guzman

The "SourceXXXXX" classes calculate source coefficients in a broad sense.

Parameters:
*/
class SourceBase : public Params {

protected:

//! Spatial data (transient)
   SpatialData _spdata;

//! Value of source term (transient)
   double SourceTerm;

//! Value of time step for point sources (transient)
   double _dt;

//! Default constructor (protected, class not designed to be instantiated)
   SourceBase(void);

//! Constructor with arguments (to speed up construction of derived classes)
   SourceBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   SourceBase(const SourceBase& other);

//! Set up the source model based on "params"
   virtual void SetupSource(bool construct);

//! Compute the source coefficients
   virtual void EvaluateSource(void);

public:

//! Destructor
   virtual ~SourceBase() = default;

//! Clone function (stub)
   virtual std::unique_ptr<SourceBase> Clone(void) const = 0;

//! Set up the class parameters
   void SetupObject(const DataContainer& cont_in);

//! Evaluate and return one source term
   double GetSource(double t_in, const GeoVector& pos_in, const GeoVector& mom_in, const SpatialData& spdata_in, double dt_in);

};

};

#endif
