/*!
\file distribution_other.hh
\brief Declares several event distribution classes
\author Vladimir Florinski
\author Juan G Alonso Guzman

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
template <typename Trajectory_, class distroClass>
class DistributionUniform : public DistributionTemplated<Trajectory_, distroClass> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using DistributionBase = DistributionBase<Trajectory>;
   using DistributionTemplated = DistributionTemplated<Trajectory, distroClass>;

protected:

   using DistributionTemplated::_status;
   using DistributionTemplated::container;
   using DistributionTemplated::ActionTable;
   using DistributionTemplated::_weight;
   using DistributionTemplated::dims;
   using DistributionTemplated::_value;
   using DistributionTemplated::_t;
   using DistributionTemplated::_t2;
   using DistributionTemplated::_mom;
   using DistributionTemplated::_mom2;
   using DistributionTemplated::_fields;
   using DistributionTemplated::_fields2;

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

//! Default constructor (protected, class not designed to be instantiated)
   DistributionUniform(void);

//! Constructor with arguments (to speed up construction of derived classes)
   DistributionUniform(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   DistributionUniform(const DistributionUniform& other);

public:

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
Parameters: (DistributionUniform), int val_time
*/
template <typename Trajectory_>
class DistributionTimeUniform : public DistributionUniform<Trajectory_, double> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using DistributionBase = DistributionBase<Trajectory>;
   using DistributionUniform = DistributionUniform<Trajectory, double>;

protected:

   using DistributionUniform::_status;
   using DistributionUniform::container;
   using DistributionUniform::ActionTable;
   using DistributionUniform::_weight;
   using DistributionUniform::dims;
   using DistributionUniform::_value;
   using DistributionUniform::_t;
   using DistributionUniform::_t2;

protected:

//! Which coordinates to use for value: 0 initial, 1 final (persistent)
   int val_time;

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
\brief A uniform function of the spatial position
\author Juan G Alonso Guzman

Type: 3D position
Parameters: (DistributionUniform), int val_time, int val_coord
*/
template <typename Trajectory_>
class DistributionPositionUniform : public DistributionUniform<Trajectory_, double> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using DistributionBase = DistributionBase<Trajectory>;
   using DistributionUniform = DistributionUniform<Trajectory, double>;

protected:

   using DistributionUniform::_status;
   using DistributionUniform::container;
   using DistributionUniform::ActionTable;
   using DistributionUniform::_weight;
   using DistributionUniform::dims;
   using DistributionUniform::_value;
   using DistributionUniform::_pos;
   using DistributionUniform::_pos2;

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

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionMomentumUniform class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DistributionMomentumUniform class
const std::string dist_name_momentum_uniform = "DistributionMomentumUniform";

/*!
\brief A uniform function of the momentum
\author Juan G Alonso Guzman

Type: 3D momentum
Parameters: (DistributionUniform), int val_time, int val_coord
*/
template <typename Trajectory_>
class DistributionMomentumUniform : public DistributionUniform<Trajectory_, double> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using DistributionBase = DistributionBase<Trajectory>;
   using DistributionUniform = DistributionUniform<Trajectory, double>;

protected:

   using DistributionUniform::_status;
   using DistributionUniform::container;
   using DistributionUniform::ActionTable;
   using DistributionUniform::_weight;
   using DistributionUniform::dims;
   using DistributionUniform::_value;
   using DistributionUniform::_mom;
   using DistributionUniform::_mom2;
   using DistributionUniform::_fields;
   using DistributionUniform::_fields2;

protected:

//   using DistributionUniform::_

protected:

//! Which coordinates to use for value: 0 initial, 1 final (persistent)
   int val_time;

//! Which coordinate representation to use for value: 0 "native coordinates", 1 locally spherical with B || z (persistent)
   int val_coord;

//! Set up the distribution accumulator based on "params"
   void SetupDistribution(bool construct) override;

//! Determine the value to be binned from a phase space position and other arguments
   void EvaluateValue(void) override;

public:

//! Default constructor
   DistributionMomentumUniform(void);

//! Copy constructor
   DistributionMomentumUniform(const DistributionMomentumUniform& other);

//! Destructor
   ~DistributionMomentumUniform() override = default;

//! Clone function
   CloneFunctionDistribution(DistributionMomentumUniform);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionPositionMomentumUniform class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DistributionPositionMomentumUniform class
const std::string dist_name_position_momentum_uniform = "DistributionPositionMomentumUniform";

/*!
\brief A uniform function of the position and momentum
\author Juan G Alonso Guzman

Type: 2D position-momentum
Parameters: (DistributionUniform), int val_time, int pos_idx, int mom_idx
*/
template <typename Trajectory_>
class DistributionPositionMomentumUniform : public DistributionUniform<Trajectory_, double> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using DistributionBase = DistributionBase<Trajectory>;
   using DistributionUniform = DistributionUniform<Trajectory, double>;

protected:

   using DistributionUniform::_status;
   using DistributionUniform::container;
   using DistributionUniform::ActionTable;
   using DistributionUniform::_weight;
   using DistributionUniform::dims;
   using DistributionUniform::_value;
   using DistributionUniform::_pos;
   using DistributionUniform::_pos2;
   using DistributionUniform::_mom;
   using DistributionUniform::_mom2;
   using DistributionUniform::_fields;
   using DistributionUniform::_fields2;

protected:

//! Which coordinates to use for value: 0 initial, 1 final (persistent)
   int val_time;

//! Which component of position to use for value[0] (persistent)
   int pos_idx;

//! Which component of momentum to use for value[1] (persistent)
   int mom_idx;

//! Set up the distribution accumulator based on "params"
   void SetupDistribution(bool construct) override;

//! Determine the value to be binned from a phase space position and other arguments
   void EvaluateValue(void) override;

public:

//! Default constructor
   DistributionPositionMomentumUniform(void);

//! Copy constructor
   DistributionPositionMomentumUniform(const DistributionPositionMomentumUniform& other);

//! Destructor
   ~DistributionPositionMomentumUniform() override = default;

//! Clone function
   CloneFunctionDistribution(DistributionPositionMomentumUniform);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionAnisotropyLISM class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DistributionAnisotropyLISM class
const std::string dist_name_anisotropy_LISM = "DistributionAnisotropyLISM";

/*!
\brief Prescribed functions of spatial position and momentum to characterize LISM anisotropy of GCRs (See Zhang et al. 2014)
\author Juan G Alonso Guzman

Type: 3D momentum
Parameters: (DistributionUniform), GeoVector[3] rot_matrix, GeoVector U_LISM, double mom_pow_law
*/
template <typename Trajectory_>
class DistributionAnisotropyLISM : public DistributionTemplated<Trajectory_, double> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using DistributionBase = DistributionBase<Trajectory>;
   using DistributionTemplated = DistributionTemplated<Trajectory, double>;

   static_assert(std::same_as<Trajectory, TrajectoryLorentz<Fields>>, "DistributionAnisotropyLISM requires Lorentz Trajectory.");

protected:

   using DistributionTemplated::_status;
   using DistributionTemplated::container;
   using DistributionTemplated::ActionTable;
   using DistributionTemplated::_weight;
   using DistributionTemplated::dims;
   using DistributionTemplated::_value;
//   using DistributionTemplated::_t;
//   using DistributionTemplated::_t2;
   using DistributionTemplated::_mom;
   using DistributionTemplated::_mom2;
   using DistributionTemplated::_pos;
   using DistributionTemplated::_pos2;
   using DistributionTemplated::_fields;
   using DistributionTemplated::_fields2;
   using DistributionTemplated::specie;

protected:

//! Rotation matrix for binning coordinate system (persistent)
   GeoVector rot_matrix[3];

//! Relative velocity of LISM (persistent)
   GeoVector U_LISM;

//! Relative momentum in the LISM frame (transient)
   GeoVector mom_rel;

//! Slope of momentum power law (persistent)
   double mom_pow_law;

//! Gradient of density perpendicular to LISM magnetic field (persistent)
   GeoVector grad_perp_dens;

//! Set up the distribution accumulator based on "params"
   void SetupDistribution(bool construct) override;

//! Determine the value to be binned from a phase space position and other arguments
   void EvaluateValue(void) override;

//! Evaluate Compton-Getting factor
   void ComptonGettingFactor(void);

//! Momentum power law anisotropy
   void MomPowerLawAnisotropy(void);

//! First Legendre polynomial of pitch angle cosine anisotropy (dipole)
   void FirstLegendreAnisotropy(void);

//! Second Legendre polynomial of pitch angle cosine anisotropy (quadrupole)
   void SecondLegendreAnisotropy(void);

//! b x gradient anisotropy (guiding center dotted with perpendicular density gradient)
   void bCrossGradientAnisotropy(void);

public:

//! Default constructor
   DistributionAnisotropyLISM(void);

//! Copy constructor
   DistributionAnisotropyLISM(const DistributionAnisotropyLISM& other);

//! Destructor
   ~DistributionAnisotropyLISM() override = default;

//! Clone function
   CloneFunctionDistribution(DistributionAnisotropyLISM);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionSpectrumKineticEnergyPowerLaw class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Flag to define power law spectrum as differential density 0, differential intensity 1, or distribution function 2
#define DISTRO_KINETIC_ENERGY_POWER_LAW_TYPE 0

//! Readable name of the DistributionSpectrumKineticEnergyPowerLaw class
const std::string dist_name_spectrum_kinetic_energy_power_law = "DistributionSpectrumKineticEnergyPowerLaw";

/*!
\brief Differential intensity as a specified function of kinetic energy J(T)
\author Vladimir Florinski

Type: 1D momentum
Parameters: (DistributionTemplated), double J0, double T0, double pow_law, double val_cold
*/
template <typename Trajectory_>
class DistributionSpectrumKineticEnergyPowerLaw : public DistributionTemplated<Trajectory_, double> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using DistributionBase = DistributionBase<Trajectory>;
   using DistributionTemplated = DistributionTemplated<Trajectory, double>;

protected:

   using DistributionTemplated::_status;
   using DistributionTemplated::container;
   using DistributionTemplated::ActionTable;
   using DistributionTemplated::_weight;
   using DistributionTemplated::dims;
   using DistributionTemplated::_value;
//   using DistributionTemplated::_t;
//   using DistributionTemplated::_t2;
   using DistributionTemplated::_mom;
   using DistributionTemplated::_mom2;
   using DistributionTemplated::_pos;
   using DistributionTemplated::_pos2;
   using DistributionTemplated::_fields;
   using DistributionTemplated::_fields2;
   using DistributionTemplated::specie;

protected:

//! Normalization for the "hot" boundary (persistent)
   double J0;

//! Characteristic energy (persistent)
   double T0;

//! Spectral power law (persistent)
   double pow_law;

//! Constant value for the "cold" condition (persistent)
   double val_cold;

//! Kinetic energy at binning (transient)
   double kin_energy;

//! Set up the distribution accumulator based on "params"
   void SetupDistribution(bool construct) override;

//! Determine the value to be binned from a phase space position and other arguments
   void EvaluateValue(void) override;

//! Weight from a "hot" boundary
   virtual void SpectrumKineticEnergyPowerLawHot(void);

//! Weight from a "cold" boundary
   void SpectrumKineticEnergyPowerLawCold(void);

public:

//! Default constructor
   DistributionSpectrumKineticEnergyPowerLaw(void);

//! Constructor with arguments (to speed up construction of derived classes)
   DistributionSpectrumKineticEnergyPowerLaw(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor
   DistributionSpectrumKineticEnergyPowerLaw(const DistributionSpectrumKineticEnergyPowerLaw& other);

//! Destructor
   ~DistributionSpectrumKineticEnergyPowerLaw() override = default;

//! Clone function
   CloneFunctionDistribution(DistributionSpectrumKineticEnergyPowerLaw);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistributionSpectrumKineticEnergyBentPowerLaw class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DistributionSpectrumKineticEnergyBentPowerLaw class
const std::string dist_name_spectrum_kinetic_energy_bent_power_law = "DistributionSpectrumKineticEnergyBentPowerLaw";

/*!
\brief Differential intensity as a specified function of kinetic energy J(T)
\author Juan G Alonso Guzman

Type: 1D momentum
Parameters: (DistributionSpectrumKineticEnergyPowerLaw), double T_b, double pow_law_b
*/
template <typename Trajectory_>
class DistributionSpectrumKineticEnergyBentPowerLaw : public DistributionSpectrumKineticEnergyPowerLaw<Trajectory_> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using DistributionBase = DistributionBase<Trajectory>;
   using DistributionSpectrumKineticEnergyPowerLaw = DistributionSpectrumKineticEnergyPowerLaw<Trajectory>;

   using DistributionSpectrumKineticEnergyPowerLaw::pow_law;
   using DistributionSpectrumKineticEnergyPowerLaw::kin_energy;
protected:

   using DistributionSpectrumKineticEnergyPowerLaw::_status;
   using DistributionSpectrumKineticEnergyPowerLaw::container;
   using DistributionSpectrumKineticEnergyPowerLaw::ActionTable;
   using DistributionSpectrumKineticEnergyPowerLaw::_weight;
   using DistributionSpectrumKineticEnergyPowerLaw::dims;
   using DistributionSpectrumKineticEnergyPowerLaw::_value;
   using DistributionSpectrumKineticEnergyPowerLaw::_mom;
   using DistributionSpectrumKineticEnergyPowerLaw::_mom2;
   using DistributionSpectrumKineticEnergyPowerLaw::_pos;
   using DistributionSpectrumKineticEnergyPowerLaw::_pos2;
   using DistributionSpectrumKineticEnergyPowerLaw::_fields;
   using DistributionSpectrumKineticEnergyPowerLaw::_fields2;
   using DistributionSpectrumKineticEnergyPowerLaw::specie;


protected:

//! Bendover energy (persistent)
   double T_b;

//! Combined power law in the denominator energy ratio (persistent)
   double pow_law_comb;

//! Factor to control bend smoothness (persistent)
   double bend_smoothness;

//! Set up the distribution accumulator based on "params"
   void SetupDistribution(bool construct) override;

//! Weight from a "hot" boundary
   void SpectrumKineticEnergyPowerLawHot(void) override;

public:

//! Default constructor
   DistributionSpectrumKineticEnergyBentPowerLaw(void);

//! Copy constructor
   DistributionSpectrumKineticEnergyBentPowerLaw(const DistributionSpectrumKineticEnergyBentPowerLaw& other);

//! Destructor
   ~DistributionSpectrumKineticEnergyBentPowerLaw() override = default;

//! Clone function
   CloneFunctionDistribution(DistributionSpectrumKineticEnergyBentPowerLaw);
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
template <typename Trajectory_>
class DistributionPositionCumulativeOrder1 : public DistributionTemplated<Trajectory_, GeoVector> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using DistributionBase = DistributionBase<Trajectory>;
   using DistributionTemplated = DistributionTemplated<Trajectory, GeoVector>;

protected:

   using DistributionTemplated::_status;
   using DistributionTemplated::container;
   using DistributionTemplated::ActionTable;
   using DistributionTemplated::_weight;
   using DistributionTemplated::dims;
   using DistributionTemplated::_value;
   using DistributionTemplated::_t;
   using DistributionTemplated::_t2;
   using DistributionTemplated::_mom;
   using DistributionTemplated::_mom2;
   using DistributionTemplated::_pos;
   using DistributionTemplated::_pos2;
   using DistributionTemplated::_fields;
   using DistributionTemplated::_fields2;
   using DistributionTemplated::specie;

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
template <typename Trajectory_>
class DistributionPositionCumulativeOrder2 : public DistributionTemplated<Trajectory_, GeoMatrix> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using DistributionBase = DistributionBase<Trajectory>;
   using DistributionTemplated = DistributionTemplated<Trajectory, GeoMatrix>;

protected:

   using DistributionTemplated::_status;
   using DistributionTemplated::container;
   using DistributionTemplated::ActionTable;
   using DistributionTemplated::_weight;
   using DistributionTemplated::dims;
   using DistributionTemplated::_value;
   using DistributionTemplated::_t;
   using DistributionTemplated::_t2;
   using DistributionTemplated::_mom;
   using DistributionTemplated::_mom2;
   using DistributionTemplated::_pos;
   using DistributionTemplated::_pos2;
   using DistributionTemplated::_fields;
   using DistributionTemplated::_fields2;
   using DistributionTemplated::specie;

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
template <typename Trajectory_>
class DistributionLossCone : public DistributionTemplated<Trajectory_, GeoVector> {
public:

   using Trajectory = Trajectory_;
   using Fields = Trajectory::Fields;
   using DistributionBase = DistributionBase<Trajectory>;
   using DistributionTemplated = DistributionTemplated<Trajectory, GeoVector>;

protected:

   using DistributionTemplated::_status;
   using DistributionTemplated::container;
   using DistributionTemplated::ActionTable;
   using DistributionTemplated::_weight;
   using DistributionTemplated::dims;
   using DistributionTemplated::_value;
   using DistributionTemplated::_mom;
   using DistributionTemplated::_mom2;
   using DistributionTemplated::_pos;
   using DistributionTemplated::_pos2;
   using DistributionTemplated::_fields;
   using DistributionTemplated::_fields2;
   using DistributionTemplated::specie;
   using DistributionTemplated::_edata;
   using DistributionTemplated::_edata2;

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

// Something like this is needed for templated classes
#include "distribution_other.cc"

#endif
