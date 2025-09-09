/*!
\file distance_map.hh
\brief Mapping between reference and physical distances
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific nThis file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DISTANCE_MAP
#define SPECTRUM_DISTANCE_MAP

#include <memory>

#include "common/status_class.hh"

namespace Spectrum {

//! Clone function pattern
#define CloneFunctionDistance(T) std::unique_ptr<DistanceBase> Clone(void) const override {return std::make_unique<T>();};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistanceBase class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Base class for a mapping function between a reference uniform coordinate and a stretched physical coordinate
\author Vladimir Florinski

Parameters: rmin, rmax
*/
class DistanceBase : public Params {

protected:

//! Largest physical distance (persistent)
   double rmin;

//! Smallest physical distance (persistent)
   double rmax;

//! Tells the evaluator which map is needed
   int selector;

//! Default constructor (protected, class not designed to be instantiated)
   DistanceBase(void);

//! Constructor with arguments (to speed up construction of derived classes)
   DistanceBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   DistanceBase(const DistanceBase& other);

//! Set up the distance map based on "params"
   virtual void SetupDistance(bool construct);

//! Compute the result of applying the map
   virtual void EvaluateDistance(void);

public:

//! Destructor
   virtual ~DistanceBase() = default;

//! Clone function (stub)
   virtual std::unique_ptr<DistanceBase> Clone(void) const = 0;

//! Set up the object's persistent class data members (generic)
   void SetupObject(const DataContainer& cont_in);

//! Return the physical coordinate (strored in _pos[0])
   double GetPhysical(double ref);

//! Return the reference coordinate (strored in _pos[1])
   double GetReference(double phys);

//! Return the derivative of the physical coordinate (strored in _pos[2])
   double GetDerivative(double ref);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistanceExponential class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DistanceExponential class
const std::string distance_name_exponential = "DistanceExponential";

/*!
\brief An exponential map
\author Vladimir Florinski

This map gives exponential scaling, where the density of points is inversely proportional to the physical distance.

Parameters: (DistanceBase)
*/
class DistanceExponential : public DistanceBase {

protected:

//! Ratio of distances
   double ratio;

//! Logarithmic ratio of distances
   double log_ratio;

//! Set up the radial map based on "params"
   void SetupDistance(bool construct) override;

//! Compute the result of applying the map
   void EvaluateDistance(void) override;

public:

//! Default constructor
   DistanceExponential(void);

//! Copy constructor
   DistanceExponential(const DistanceExponential& other);

//! Destructor
   ~DistanceExponential() override = default;

//! Clone function
   CloneFunctionDistance(DistanceExponential);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistancePowerLaw class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the RadialPowerLaw class
const std::string distance_name_power_law = "DistancePowerLaw";

/*!
\brief A power law map
\author Vladimir Florinski

This map generates a scaling controlled by the pareameter "powl". A value of 1 gives linear scaling, and a large number gives an exponential scaling. Small "powl" (<1) yield a map with inverted point density.

Parameters: (DistanceBase), powl
*/
class DistancePowerLaw : public DistanceBase {

protected:

//! Power law index (persistent)
   double powl;

//! Scaling for the reference coordinate (persistent)
   double ref0;

//! Inverse of the power law index (persistent)
   double pow_inv;

//! Coefficient for derivative calculation (persistent)
   double c;

//! Set up the radial map based on "params"
   void SetupDistance(bool construct) override;

//! Compute the result of applying the map
   void EvaluateDistance(void) override;

public:

//! Default constructor
   DistancePowerLaw(void);

//! Copy constructor
   DistancePowerLaw(const DistancePowerLaw& other);

//! Destructor
   ~DistancePowerLaw() override = default;

//! Clone function
   CloneFunctionDistance(DistancePowerLaw);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistanceLinearExp class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DistanceLinearExp class
const std::string distance_name_linear_exp = "DistanceLinearExp";

/*!
\brief A linear map smoothly transitioning into an exponential map
\author Vladimir Florinski

This map consists of two regions meeting at "rmed". The argument "frac" is a fraction of the points that fall within the inner region. The argument "logy" specifies the smoothness of the transition (the larger it is the sharper the transition). For large "logy" (>5-10) the inner region is linear and the outer region is exponential. Small "logy" (<=1) can yield a map that is (almost) exponential everywhere or a map with inverted point density.

Parameters: (DistanceBase), rmed, frac, logy
*/
class DistanceLinearExp : public DistanceBase {

protected:

//! Base for exponentiation (persistent)
   double base;

//! Logarithm of the base (persistent)
   double logy;

//! Scaling for the reference coordinate (persistent)
   double ref0;

//! Coefficient for the exponential (persistent)
   double C;

//! Coefficient for derivative calculation (persistent)
   double c;

//! Set up the radial map based on "params"
   void SetupDistance(bool construct) override;

//! Compute the result of applying the map
   void EvaluateDistance(void) override;

public:

//! Default constructor
   DistanceLinearExp(void);

//! Copy constructor
   DistanceLinearExp(const DistanceLinearExp& other);

//! Destructor
   ~DistanceLinearExp() override = default;

//! Clone function
   CloneFunctionDistance(DistanceLinearExp);
};

#ifdef GEO_DEBUG
//! Testing routine
void PrintDistanceMaps(DistanceBase* dist_map, double powl, bool as_phys);
#endif

};

#endif
