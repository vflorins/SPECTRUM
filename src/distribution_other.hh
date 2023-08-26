/*!
\file distribution_other.hh
\brief Declares several event distribution classes
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DISTRIBUTION_OTHER_HH
#define SPECTRUM_DISTRIBUTION_OTHER_HH

#include "distribution_templated.hh"

namespace Spectrum {

//! Clone function pattern
#define CloneFunctionDistribution(T) std::shared_ptr<DistributionBase> Clone(void) const override {return std::make_shared<T>();};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionUniform class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A parent class for uniform distributions
\author Juan G Alonso Guzman
\author Vladimir Florinski

Type: Unspecified
Parameters: (DistributionTemplated), double val_hot, double val_cold
*/
template <class distroClass>
class DistributionUniform : public DistributionTemplated<distroClass> {

protected:

//! Constant value for the "hot" condition (persistent)
   distroClass val_hot;

//! Constant value for the "cold" condition (persistent)
   distroClass val_cold;

//! Set up the distribution accumulator based on "params"
   void SetupDistribution(bool construct) override;

//! Weight from a "hot" boundary
   void UniformHot(void);

//! Weight from a "cold" boundary
   void UniformCold(void);

public:

//! Default constructor (never meant to be used)
   DistributionUniform(void) = delete;

//! Constructor with arguments (to speed up construction of derived classes)
   DistributionUniform(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor
   DistributionUniform(const DistributionUniform& other);

//! Destructor
   ~DistributionUniform() override = default;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionTimeUniform class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DistributionTimeUniform class
const std::string dist_name_time_uniform = "DistributionTimeUniform";

/*!
\brief A uniform function of the simulated time
\author Juan G Alonso Guzman

Type: 1D time
Parameters: (DistributionUniform)
*/
class DistributionTimeUniform : public DistributionUniform<double> {

protected:

//! Set up the distribution accumulator based on "params"
   void SetupDistribution(bool construct) override;

//! Determine the value to be binned from a phase space position and other arguments
   void EvaluateValue(void) override;

public:

//! Default constructor
   DistributionTimeUniform(void);

//! Copy constructor
   DistributionTimeUniform(const DistributionTimeUniform& other);

//! Destructor
   ~DistributionTimeUniform() override = default;

//! Clone function
   CloneFunctionDistribution(DistributionTimeUniform);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionPositionUniform class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DistributionPositionUniform class
const std::string dist_name_position_uniform = "DistributionPositionUniform";

/*!
\brief An isotropic function of the spatial position
\author Juan G Alonso Guzman

Type: 3D position
Parameters: (DistributionUniform), int val_coord, int val_time
*/
class DistributionPositionUniform : public DistributionUniform<double> {

protected:

//! Which coordinates to use for value: 0 initial, 1 final (persistent)
   int val_time;

//! Which coordinate representation to use for value: 0 cartesian, 1 spherical (persisent)
   int val_coord;

//! Set up the distribution accumulator based on "params"
   void SetupDistribution(bool construct) override;

//! Determine the value to be binned from a phase space position and other arguments
   void EvaluateValue(void) override;

public:

//! Default constructor
   DistributionPositionUniform(void);

//! Copy constructor
   DistributionPositionUniform(const DistributionPositionUniform& other);

//! Destructor
   ~DistributionPositionUniform() override = default;

//! Clone function
   CloneFunctionDistribution(DistributionPositionUniform);
};

#if TRAJ_TYPE != TRAJ_PARKER

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionPitchUniform class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DistributionPitchUniform class
const std::string dist_name_pitch_uniform = "DistributionPitchUniform";

/*!
\brief An isotropic function of the pitch angle cosine (mu)
\author Vladimir Florinski

Type: 1D momentum
Parameters: (DistributionUniform), int val_time
*/
class DistributionPitchUniform : public DistributionUniform<double> {

protected:

//! Which coordinates to use for value: 0 initial, 1 final (persistent)
   int val_time;

//! Set up the distribution accumulator based on "params"
   void SetupDistribution(bool construct) override;

//! Determine the value to be binned from a phase space position and other arguments
   void EvaluateValue(void) override;

public:

//! Default constructor
   DistributionPitchUniform(void);

//! Copy constructor
   DistributionPitchUniform(const DistributionPitchUniform& other);

//! Destructor
   ~DistributionPitchUniform() override = default;

//! Clone function
   CloneFunctionDistribution(DistributionPitchUniform);
};

#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionSpectrumLISM class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DistributionSpectrumLISM class
const std::string dist_name_spectrum_lism = "DistributionSpectrumLISM";

/*!
\brief Differential intensity as a specified function of kinetic energy J(T)
\author Vladimir Florinski

Type: 1D momentum
Parameters: (DistributionTemplated), double J0, double T0, double pow_law, double val_cold
*/
class DistributionSpectrumLISM : public DistributionTemplated<double> {

protected:

//! Normalization for the "hot" boundary (persistent)
   double J0;

//! Characteristic energy (persistent)
   double T0;

//! Spectral power law (persistent)
   double pow_law;

//! Constant value for the "cold" condition (persistent)
   double val_cold;

//! Set up the distribution accumulator based on "params"
   void SetupDistribution(bool construct) override;

//! Determine the value to be binned from a phase space position and other arguments
   void EvaluateValue(void) override;

//! Weight from a "hot" boundary
   void SpectrumLISMHot(void);

//! Weight from a "cold" boundary
   void SpectrumLISMCold(void);

public:

//! Default constructor
   DistributionSpectrumLISM(void);

//! Copy constructor
   DistributionSpectrumLISM(const DistributionSpectrumLISM& other);

//! Destructor
   ~DistributionSpectrumLISM() override = default;

//! Clone function
   CloneFunctionDistribution(DistributionSpectrumLISM);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionPositionCumulativeOrder1 class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DistributionPositionCumulativeOrder1 class
const std::string dist_name_pos_cumulative_O1 = "DistributionPositionCumulativeOrder1";

/*!
\brief First order spatial cumulatives <x>, <y>, <z>
\author Juan G Alonso Guzman
\author Vladimir Florinski

Type: 1D time
Parameters: (DistributionTemplated)
*/
class DistributionPositionCumulativeOrder1 : public DistributionTemplated<GeoVector> {

protected:

//! Set up the distribution accumulator based on "params"
   void SetupDistribution(bool construct) override;

//! Determine the value to be binned from a phase space position and other arguments
   void EvaluateValue(void) override;

//! Record position for cumulative
   void RecordPosition(void);

public:

//! Default constructor
   DistributionPositionCumulativeOrder1(void);

//! Copy constructor
   DistributionPositionCumulativeOrder1(const DistributionPositionCumulativeOrder1& other);

//! Destructor
   ~DistributionPositionCumulativeOrder1() override = default;

//! Clone function
   CloneFunctionDistribution(DistributionPositionCumulativeOrder1);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionPositionCumulativeOrder2 class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DistributionPositionCumulativeOrder2 class
const std::string dist_name_pos_cumulative_O2 = "DistributionPositionCumulativeOrder2";

/*!
\brief Second order spatial cumulatives <xx>, <xy>, <xz>, etc.
\author Vladimir Florinski
\author Juan G Alonso Guzman

Type: 1D time
Parameters: (DistributionTemplated)
*/
class DistributionPositionCumulativeOrder2 : public DistributionTemplated<GeoMatrix> {

protected:

//! Set up the distribution accumulator based on "params"
   void SetupDistribution(bool construct) override;

//! Determine the value to be binned from a phase space position and other arguments
   void EvaluateValue(void) override;

//! Record position for cumulative
   void RecordPosition(void);

public:

//! Default constructor
   DistributionPositionCumulativeOrder2(void);

//! Copy constructor
   DistributionPositionCumulativeOrder2(const DistributionPositionCumulativeOrder2& other);

//! Destructor
   ~DistributionPositionCumulativeOrder2() override = default;

//! Clone function
   CloneFunctionDistribution(DistributionPositionCumulativeOrder2);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionLossCone class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DistributionLossCone class
const std::string dist_name_loss_cone = "DistributionLossCone";

/*!
\brief Minimum |B|, maximum |B|, and minimum loss cone along a trajectory
\author Juan G Alonso Guzman

Type: 3D position
Parameters: (DistributionTemplated), int val_time, int val_coord
*/
class DistributionLossCone : public DistributionTemplated<GeoVector> {

protected:

//! Which coordinates to use for value: 0 initial, 1 final (persistent)
   int val_time;

//! Which coordinate representation to use for value: 0 cartesian, 1 spherical (persisent)
   int val_coord;

//! Set up the distribution accumulator based on "params"
   void SetupDistribution(bool construct) override;

//! Determine the value to be binned from a phase space position and other arguments
   void EvaluateValue(void) override;

//! Record minimum |B|, maximum |B|, and cumulative loss cone half-angle, asin(sqrt(|B|_min/|B|_max)), along trajectory
   void RecordLossCone(void);

public:

//! Default constructor
   DistributionLossCone(void);

//! Copy constructor
   DistributionLossCone(const DistributionLossCone& other);

//! Destructor
   ~DistributionLossCone() override = default;

//! Clone function
   CloneFunctionDistribution(DistributionLossCone);
};

};

#endif
