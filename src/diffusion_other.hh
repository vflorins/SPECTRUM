/*!
\file diffusion_other.hh
\brief Declares several classes to compute difusion coefficients
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DIFFUSION_OTHER_HH
#define SPECTRUM_DIFFUSION_OTHER_HH

#include "diffusion_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionNone class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief No diffusion class, for reduced physics
\author Lucius Schoenbaum

Parameters: (DiffusionBase) None
*/
template <typename HConfig_>
class DiffusionNone : public DiffusionBase<HConfig_> {

//! Readable name of the DiffusionIsotropicConstant class
   static constexpr std::string_view diff_name = "DiffusionNone";

public:

   using HConfig = HConfig_;
   using Config = HConfig::DiffusionConfig;

   using DiffusionBase = DiffusionBase<HConfig>;
   using Coordinates = DiffusionBase::Coordinates;
   using Fields = DiffusionBase::Fields;

   using DiffusionBase::_status;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::Kappa;
//   using DiffusionBase::mu;
   using DiffusionBase::Stage;

protected:

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(Component comp) override;

public:

//! Default constructor
   DiffusionNone(void);

//! Copy constructor
   DiffusionNone(const DiffusionNone& other);

//! Destructor
   ~DiffusionNone() override = default;

//! Clone function
//   CloneFunctionDiffusion(DiffusionNone);

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(Component comp) override;

};


//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionIsotropicConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Isotropic scattering in pitch angle only, uniform in space
\author Vladimir Florinski

Parameters: (DiffusionBase), double D0
*/
template <typename HConfig_>
class DiffusionIsotropicConstant : public DiffusionBase<HConfig_> {

//! Readable name of the DiffusionIsotropicConstant class
   static constexpr std::string_view diff_name = "DiffusionIsotropicConstant";

public:

   using HConfig = HConfig_;
   static_assert(HConfig::trajectory != Config::Trajectory::Parker, "DiffusionIsotropicConstant diffusion type cannot be applied to the Parker Trajectory type.");

   using Config = HConfig::DiffusionConfig;

   using DiffusionBase = DiffusionBase<HConfig>;
   using Coordinates = DiffusionBase::Coordinates;
   using Fields = DiffusionBase::Fields;

   using DiffusionBase::_status;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::Kappa;
//   using DiffusionBase::mu;
   using DiffusionBase::Stage;


protected:

//! Scattering frequency (persistent)
   static constexpr double D0 = Config::D0;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(Component comp) override;

public:

//! Default constructor
   DiffusionIsotropicConstant(void);

//! Copy constructor
   DiffusionIsotropicConstant(const DiffusionIsotropicConstant& other);

//! Destructor
   ~DiffusionIsotropicConstant() override = default;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(Component comp) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionQLTConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Scattering in pitch angle according to quasi-linear theory for slab turbulence, uniform in space
\author Vladimir Florinski

Parameters: (DiffusionBase), double A2A, double l_max, double ps_index
*/
template <typename HConfig_>
class DiffusionQLTConstant : public DiffusionBase<HConfig_> {

//! Readable name of the DiffusionQLTConstant class
   static constexpr std::string_view diff_name = "DiffusionQLTConstant";

public:

   using HConfig = HConfig_;
   static_assert(HConfig::trajectory != Config::Trajectory::Parker, "DiffusionIsotropicConstant diffusion type cannot be applied to the Parker Trajectory type.");
   using Config = HConfig::DiffusionConfig;

   using DiffusionBase = DiffusionBase<HConfig>;
   using Coordinates = DiffusionBase::Coordinates;
   using Fields = DiffusionBase::Fields;

//   using DiffusionBase::_status;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::Kappa;
   using DiffusionBase::Omega;
   using DiffusionBase::st2;
   using DiffusionBase::Stage;

protected:

//! Alfven turbulence relative variance (persistent)
   static constexpr double A2A = Config::A2A;

//! Maximum turbulent lengthscale (persistent)
   static constexpr double l_max = Config::l_max;

//! Characteristic wavenumber (persistent)
   static constexpr double k_min = M_2PI / l_max;

//! Power spectral index (persistent)
   static constexpr double ps_index = Config::ps_index;

//! Power spectral index minus one (persistent)
   static constexpr double ps_minus = ps_index - 1.0;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(Component comp) override;

public:

//! Default constructor
   DiffusionQLTConstant(void);

//! Constructor with arguments (to speed up construction of derived classes)
   DiffusionQLTConstant(const std::string_view& name_in, status_t status_in);

//! Copy constructor
   DiffusionQLTConstant(const DiffusionQLTConstant& other);

//! Destructor
   ~DiffusionQLTConstant() override = default;

};


//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionWNLTConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A weakly nonlinear model with pitch angle scattering and perpendicular diffusion with constant turbulent ratios.
\author Vladimir Florinski

Parameters: (DiffusionQLTConstant), double A2T, double A2L
*/
template <typename HConfig_>
class DiffusionWNLTConstant : public DiffusionQLTConstant<HConfig_> {

//! Readable name of the DiffusionWNLTConstant class
   static constexpr std::string_view diff_name = "DiffusionWNLTConstant";

public:

   using HConfig = HConfig_;
   static_assert(HConfig::trajectory != Config::Trajectory::Parker, "DiffusionIsotropicConstant diffusion type cannot be applied to the Parker Trajectory type.");

   using Config = HConfig::DiffusionConfig;

   using DiffusionBase = DiffusionBase<HConfig>;
   using Coordinates = DiffusionBase::Coordinates;
   using Fields = DiffusionBase::Fields;

   using DiffusionQLTConstant = DiffusionQLTConstant<HConfig>;

   //   using DiffusionBase::_status;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::Kappa;
   using DiffusionBase::Omega;
   using DiffusionBase::st2;
   using DiffusionBase::Stage;

   using DiffusionQLTConstant::ps_index;
   using DiffusionQLTConstant::k_min;
   using DiffusionQLTConstant::ps_minus;

protected:

//! Whether to use QLT scattering
   static constexpr auto use_qlt_scatt = Config::use_qlt_scatt;

//! Transverse turbulence relative variance (persistent)
   static constexpr double A2T = Config::A2T;

//! Longitudinal turbulence relative variance (persistent)
   static constexpr double A2L = Config::A2L;

//! Power spectral index plus one (persistent)
   static constexpr double ps_plus = ps_index - 1.0;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(Component comp) override;

public:

//! Default constructor
   DiffusionWNLTConstant(void);

//! Constructor with arguments (to speed up construction of derived classes)
   DiffusionWNLTConstant(const std::string_view& name_in, status_t status_in);

//! Copy constructor
   DiffusionWNLTConstant(const DiffusionWNLTConstant& other);

//! Destructor
   ~DiffusionWNLTConstant() override = default;

};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionWNLTRampVLISM class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief WNLT diffusion in VLISM where the largest scale to turbulence increases linearly (ramp) between the HP (Rankine half body) and the sheath surface (also a Rankine half body) and is constant elsewhere.
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionWNLTConstant)
*/
template <typename HConfig_>
class DiffusionWNLTRampVLISM : public DiffusionWNLTConstant<HConfig_> {

//! Readable name of the DiffusionWNLTRampVLISM class
   static constexpr std::string_view diff_name = "DiffusionWNLTRampVLISM";

public:

   using HConfig = HConfig_;
   static_assert(HConfig::trajectory != Config::Trajectory::Parker, "DiffusionIsotropicConstant diffusion type cannot be applied to the Parker Trajectory type.");

   using Config = HConfig::DiffusionConfig;

   using DiffusionBase = DiffusionBase<HConfig>;
   using Coordinates = DiffusionBase::Coordinates;
   using Fields = DiffusionBase::Fields;

   using DiffusionWNLTConstant = DiffusionWNLTConstant<HConfig>;

   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::Stage;

   using DiffusionWNLTConstant::ps_index;
   using DiffusionWNLTConstant::ps_minus;
   using DiffusionWNLTConstant::use_qlt_scatt;
// parameters from WNLTConstant that become transient
//   using DiffusionWNLTConstant::l_max;
//   using DiffusionWNLTConstant::k_min;
//   using DiffusionWNLTConstant::A2A;
//   using DiffusionWNLTConstant::A2T;
//   using DiffusionWNLTConstant::A2L;

protected:

   //! Maximum turbulent lengthscale (persistent)
   static constexpr double l_max_ref = Config::l_max_ref;

//! Reference characteristic wavenumber (persistent)
   static constexpr double k_min_ref = M_2PI / l_max_ref;

//! Alfven turbulence relative variance (persistent)
   static constexpr double A2A_ref = Config::A2A_ref;

//! Transverse turbulence relative variance (persistent)
   static constexpr double A2T_ref = Config::A2T_ref;

//! Longitudinal turbulence relative variance (persistent)
   static constexpr double A2L_ref = Config::A2L_ref;

//! Largest scale of turbulence at HP (persistent)
   static constexpr double l_max_HP = Config::l_max_HP;

//! Difference between largest scales (persistent)
   static constexpr double dl_max = l_max_ref - l_max_HP;

//! Extent of the HP in the nose direction (persistent)
   static constexpr double z_nose = Config::z_nose;

//! Extent of the sheath in the nose direction (persistent)
   static constexpr double z_sheath = Config::z_sheath;

//! Difference between nose distance (persistent)
   static constexpr double dz = z_sheath - z_nose;

//! parameters from WNLTConstant that become transient
   double l_max = l_max_ref;
   double k_min = k_min_ref;
   double A2A = A2A_ref;
   double A2T = A2T_ref;
   double A2L = A2L_ref;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(Component comp) override;

public:

//! Default constructor
   DiffusionWNLTRampVLISM(void);

//! Copy constructor
   DiffusionWNLTRampVLISM(const DiffusionWNLTRampVLISM& other);

//! Destructor
   ~DiffusionWNLTRampVLISM() override = default;

};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionParaConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Parallel diffusion uniform in space
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionBase), double D0
*/
template <typename HConfig_>
class DiffusionParaConstant : public DiffusionBase<HConfig_> {

//! Readable name of the DiffusionParaConstant class
   static constexpr std::string_view diff_name = "DiffusionParaConstant";

public:

   using HConfig = HConfig_;
   using Config = HConfig::DiffusionConfig;

   using DiffusionBase = DiffusionBase<HConfig>;
   using Coordinates = DiffusionBase::Coordinates;
   using Fields = DiffusionBase::Fields;

   using DiffusionBase::_status;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Diffusion coefficient (persistent)
   static constexpr double D0 = Config::D0;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(Component comp) override;

public:

//! Default constructor
   DiffusionParaConstant(void);

//! Copy constructor
   DiffusionParaConstant(const DiffusionParaConstant& other);

//! Destructor
   ~DiffusionParaConstant() override = default;

//! Compute derivative of diffusion coefficient in position or time
   double GetDirectionalDerivative(Component comp, int xyz, const DerivativeData& ddata) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(Component comp) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionPerpConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Perpendicular diffusion uniform in space
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionBase), double D0
*/
template <typename HConfig_>
class DiffusionPerpConstant : public DiffusionBase<HConfig_> {

//! Readable name of the DiffusionPerpConstant class
   static constexpr std::string_view diff_name = "DiffusionPerpConstant";

public:

   using HConfig = HConfig_;
   using Config = HConfig::DiffusionConfig;

   using DiffusionBase = DiffusionBase<HConfig>;
   using Coordinates = DiffusionBase::Coordinates;
   using Fields = DiffusionBase::Fields;

   using DiffusionBase::_status;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Diffusion coefficient (persistent)
   static constexpr double D0 = Config::D0;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(Component comp) override;

public:

//! Default constructor
   DiffusionPerpConstant(void);

//! Copy constructor
   DiffusionPerpConstant(const DiffusionPerpConstant& other);

//! Destructor
   ~DiffusionPerpConstant() override = default;

//! Compute derivative of diffusion coefficient in position or time
   double GetDirectionalDerivative(Component comp, int xyz, const DerivativeData& ddata) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(Component comp) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionFullConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Full (perpendicular + parallel) diffusion uniform in space
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionBase), double Dperp, double Dpara
*/
template <typename HConfig_>
class DiffusionFullConstant : public DiffusionBase<HConfig_> {

//! Readable name of the DiffusionFullConstant class
   static constexpr std::string_view diff_name = "DiffusionFullConstant";

public:

   using HConfig = HConfig_;
   using Config = HConfig::DiffusionConfig;

   using DiffusionBase = DiffusionBase<HConfig>;
   using Coordinates = DiffusionBase::Coordinates;
   using Fields = DiffusionBase::Fields;

   using DiffusionBase::_status;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Perpendicular diffusion coefficient (persistent)
   static constexpr double Dperp = Config::Dperp;

//! Parallel diffusion coefficient (persistent)
   static constexpr double Dpara = Config::Dpara;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(Component comp) override;

public:

//! Default constructor
   DiffusionFullConstant(void);

//! Copy constructor
   DiffusionFullConstant(const DiffusionFullConstant& other);

//! Destructor
   ~DiffusionFullConstant() override = default;

//! Compute derivative of diffusion coefficient in position or time
   double GetDirectionalDerivative(Component comp, int xyz, const DerivativeData& ddata) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(Component comp) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionFlowMomentumPowerLaw class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Diffusion as a power law of flow velocity magnitude and momentum
\author Juan G Alonso Guzman
\author Swati Sharma

Parameters: (DiffusionBase), double kap0, double u0, double power_law_U, double p0, double power_law_p, double kap_rat
*/
template <typename HConfig_>
class DiffusionFlowMomentumPowerLaw : public DiffusionBase<HConfig_> {

//! Readable name of the DiffusionFlowMomentumPowerLaw class
   static constexpr std::string_view diff_name = "DiffusionFlowMomentumPowerLaw";

public:

   using HConfig = HConfig_;
   using Config = HConfig::DiffusionConfig;

   using DiffusionBase = DiffusionBase<HConfig>;
   using Coordinates = DiffusionBase::Coordinates;
   using Fields = DiffusionBase::Fields;

   using DiffusionBase::_status;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Reference diffusion coefficient (persistent)
   static constexpr double kap0 = Config::kappa0;

//! Flow velocity normalization factor (persistent)
   static constexpr double U0 = Config::U0;

//! Power law slope for flow velocity (persistent)
   static constexpr double pow_law_U = Config::pow_law_U;
   
//! Momentum normalization factor (persistent)
   static constexpr double p0 = Config::p0;

//! Power law slope for momentum (persistent)
   static constexpr double pow_law_p = Config::pow_law_p;

//! Ratio of perpendicular to parallel diffusion (persistent)
   static constexpr double kap_rat = Config::kappa_ratio;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(Component comp) override;

public:

//! Default constructor
   DiffusionFlowMomentumPowerLaw(void);

//! Copy constructor
   DiffusionFlowMomentumPowerLaw(const DiffusionFlowMomentumPowerLaw& other);

//! Destructor
   ~DiffusionFlowMomentumPowerLaw() override = default;

//! Compute derivative of diffusion coefficient in position or time
   double GetDirectionalDerivative(Component comp, int xyz, const DerivativeData& ddata) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(Component comp) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionKineticEnergyRadialDistancePowerLaw class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Full (perpendicular + parallel) diffusion, kinetic energy and radial distance power law
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionBase), double kap0, double T0, double r0, double pow_law_T, double pow_law_r, double kap_rat
*/
template <typename HConfig_>
class DiffusionKineticEnergyRadialDistancePowerLaw : public DiffusionBase<HConfig_> {

//! Readable name of the DiffusionKineticEnergyRadialDistancePowerLaw class
   static constexpr std::string_view diff_name = "DiffusionKineticEnergyRadialDistancePowerLaw";

public:

   using HConfig = HConfig_;
   using Config = HConfig::DiffusionConfig;

   using DiffusionBase = DiffusionBase<HConfig>;
   using Coordinates = DiffusionBase::Coordinates;
   using Fields = DiffusionBase::Fields;

   using DiffusionBase::_status;
   using DiffusionBase::specie;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Diffusion coefficient normalization factor (persistent)
   static constexpr double kap0 = Config::kappa0;

//! Kinetic Energy normalization factor (persistent)
   static constexpr double T0 = Config::T0;

//! Radial distance normalization factor (persistent)
   static constexpr double r0 = Config::r0;

//! Power law slope for kinetic energy (persistent)
   static constexpr double pow_law_T = Config::pow_law_T;

//! Power law slope for radial distance (persistent)
   static constexpr double pow_law_r = Config::pow_law_r;

//! Ratio of perpendicular to parallel diffusion (persistent)
   static constexpr double kap_rat = Config::kappa_ratio;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(Component comp) override;

public:

//! Default constructor
   DiffusionKineticEnergyRadialDistancePowerLaw(void);

//! Copy constructor
   DiffusionKineticEnergyRadialDistancePowerLaw(const DiffusionKineticEnergyRadialDistancePowerLaw& other);

//! Destructor
   ~DiffusionKineticEnergyRadialDistancePowerLaw() override = default;

//! Compute derivative of diffusion coefficient in position or time
   double GetDirectionalDerivative(Component comp, int xyz, const DerivativeData& ddata) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(Component comp) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionRigidityMagneticFieldPowerLaw class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Full (perpendicular + parallel) diffusion, rigidity and magnetic field power law
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionBase), double lam0, double R0, double B0, double pow_law_R, double pow_law_B, double kap_rat
*/
template <typename HConfig_>
class DiffusionRigidityMagneticFieldPowerLaw : public DiffusionBase<HConfig_> {

//! Readable name of the DiffusionRigidityMagneticFieldPowerLaw class
   static constexpr std::string_view diff_name = "DiffusionRigidityMagneticFieldPowerLaw";

public:

   using HConfig = HConfig_;
   using Config = HConfig::DiffusionConfig;

   using DiffusionBase = DiffusionBase<HConfig>;
   using Coordinates = DiffusionBase::Coordinates;
   using Fields = DiffusionBase::Fields;

   using DiffusionBase::_status;
   using DiffusionBase::specie;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Parallel mean free path (persistent)
   static constexpr double lam0 = Config::lam0;

//! Rigidity normalization factor (persistent)
   static constexpr double R0 = Config::R0;

//! Magnetic field normalization factor (persistent)
   static constexpr double B0 = Config::B0;

//! Power law slope for rigidity (persistent)
   static constexpr double pow_law_R = Config::pow_law_R;

//! Power law slope for magnetic field (persistent)
   static constexpr double pow_law_B = Config::pow_law_B;

//! Ratio of perpendicular to parallel diffusion (persistent)
   static constexpr double kap_rat = Config::kappa_ratio;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(Component comp) override;

public:

//! Default constructor
   DiffusionRigidityMagneticFieldPowerLaw(void);

//! Copy constructor
   DiffusionRigidityMagneticFieldPowerLaw(const DiffusionRigidityMagneticFieldPowerLaw& other);

//! Destructor
   ~DiffusionRigidityMagneticFieldPowerLaw() override = default;

//! Compute derivative of diffusion coefficient in position or time
   double GetDirectionalDerivative(Component comp, int xyz, const DerivativeData& ddata) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(Component comp) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionStraussEtAl2013 class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Full (perpendicular + parallel) diffusion, rigidity and magnetic field power law according to Strauss et al 2013
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionBase), int LISM_idx, double lam_in, double lam_out, double R0, double B0, double kap_rat_in, double kap_rat_out
*/
template <typename HConfig_>
class DiffusionStraussEtAl2013 : public DiffusionBase<HConfig_> {

//! Readable name of the DiffusionStraussEtAl2013 class
   static constexpr std::string_view diff_name = "DiffusionStraussEtAl2013";

public:

   using HConfig = HConfig_;
   using Config = HConfig::DiffusionConfig;

   using DiffusionBase = DiffusionBase<HConfig>;
   using Coordinates = DiffusionBase::Coordinates;
   using Fields = DiffusionBase::Fields;

   using DiffusionBase::_status;
   using DiffusionBase::specie;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Index for LISM indicator variable (persistent)
   static constexpr int LISM_idx = Config::LISM_idx;

//! Parallel inner heliosphere mean free path (persistent)
   static constexpr double lam_in = Config::lam_inner;

//! Parallel outer heliosphere mean free path (persistent)
   static constexpr double lam_out = Config::lam_outer;

//! Rigidity normalization factor (persistent)
   static constexpr double R0 = Config::R0;

//! Magnetic field normalization factor for inner heliosphere (persistent)
   static constexpr double B0 = Config::B0;

//! Ratio of perpendicular to parallel diffusion inner heliosphere (persistent)
   static constexpr double kap_rat_in = Config::kappa_ratio_inner;

//! Ratio of perpendicular to parallel diffusion outer heliosphere (persistent)
   static constexpr double kap_rat_out = Config::kappa_ratio_outer;

//! LISM indicator variable: 0 means inside HP, 1 means outside HP (transient)
   double LISM_ind;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(Component comp) override;

public:

//! Default constructor
   DiffusionStraussEtAl2013(void);

//! Copy constructor
   DiffusionStraussEtAl2013(const DiffusionStraussEtAl2013& other);

//! Destructor
   ~DiffusionStraussEtAl2013() override = default;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(Component comp) override;
};

template <typename HConfig>
using DiffusionGuoEtAl2014 = DiffusionStraussEtAl2013<HConfig>;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionPotgieterEtAl2015 class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Full (perpendicular + parallel) diffusion, rigidity and magnetic field power law according to Potgieter et al 2015
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionBase), int LISM_idx, double kappa_in, double kappa_out, double R0, double B0, double kap_rat_in, double kap_rat_out
*/
template <typename HConfig_>
class DiffusionPotgieterEtAl2015 : public DiffusionBase<HConfig_> {

//! Readable name of the DiffusionPotgieterEtAl2015 class
   static constexpr std::string_view diff_name = "DiffusionPotgieterEtAl2015";

public:

   using HConfig = HConfig_;
   using Config = HConfig::DiffusionConfig;

   using DiffusionBase = DiffusionBase<HConfig>;
   using Coordinates = DiffusionBase::Coordinates;
   using Fields = DiffusionBase::Fields;

   using DiffusionBase::_status;
   using DiffusionBase::specie;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Index for LISM indicator variable (persistent)
   static constexpr int LISM_idx = Config::LISM_idx;

//! Parallel inner heliosphere diffusion coefficient (persistent)
   static constexpr double kappa_in = Config::kappa_inner;

//! Parallel outer heliosphere diffusion coefficient (persistent)
   static constexpr double kappa_out = Config::kappa_outer;

//! Rigidity normalization factor (persistent)
   static constexpr double R0 = Config::R0;

//! Magnetic field normalization factor for inner heliosphere (persistent)
   static constexpr double B0 = Config::B0;

//! Ratio of perpendicular to parallel diffusion inner heliosphere (persistent)
   static constexpr double kap_rat_in = Config::kappa_ratio_inner;

//! Ratio of perpendicular to parallel diffusion outer heliosphere (persistent)
   static constexpr double kap_rat_out = Config::kappa_ratio_outer;

//! LISM indicator variable: 0 means inside HP, 1 means outside HP (transient)
   double LISM_ind;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(Component comp) override;

public:

//! Default constructor
   DiffusionPotgieterEtAl2015(void);

//! Copy constructor
   DiffusionPotgieterEtAl2015(const DiffusionPotgieterEtAl2015& other);

//! Destructor
   ~DiffusionPotgieterEtAl2015() override = default;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(Component comp) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionEmpiricalSOQLTandUNLT class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Full (perpendicular + parallel) diffusion, parametrized empirical fit to SOQLT (parallel) and UNLT (perpendicular) using a bent power-law spectrum with change in perpendicular diffusion according to magnetic mixing and solar cycle indicator variables
\author Juan G Alonso Guzman
Parameters: (DiffusionBase), double lam_para, double lam_perp, double R0, double B0, int Bmix_idx, double kap_rat_red, double radial_limit_perp_red, int solar_cycle_idx, double solar_cycle_effect
*/
template <typename HConfig_>
class DiffusionEmpiricalSOQLTandUNLT : public DiffusionBase<HConfig_> {

//! Readable name of the DiffusionEmpiricalSOQLTandUNLT class
   static constexpr std::string_view diff_name = "DiffusionEmpiricalSOQLTandUNLT";

public:

   using HConfig = HConfig_;
   using Config = HConfig::DiffusionConfig;

   using DiffusionBase = DiffusionBase<HConfig>;
   using Coordinates = DiffusionBase::Coordinates;
   using Fields = DiffusionBase::Fields;

   using DiffusionBase::_status;
   using DiffusionBase::specie;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Parallel mean free path (persistent)
   static constexpr double lam_para = Config::lam_para;

//! Perpendicular mean free path (persistent)
   static constexpr double lam_perp = Config::lam_perp;

//! Rigidity normalization factor (persistent)
   static constexpr double R0 = Config::R0;

//! Magnetic field normalization factor for inner heliosphere (persistent)
   static constexpr double B0 = Config::B0;

//! Index for magnetic mixing indicator variable (persistent)
   static constexpr int Bmix_idx = Config::Bmix_idx;

//! Reduction factor for kappa in unipolar regions (persistent)
   static constexpr double kap_rat_red = Config::kappa_ratio_red;

//! Radial limit to apply unipolar reduction factor (persistent)
   static constexpr double radial_limit_perp_red = Config::radial_limit_perp_red;

//! Index for solar cycle indicator variable (persistent)
   static constexpr int solar_cycle_idx = Config::solar_cycle_idx;
   
//! Solar cycle effect constant (persistent)
   static constexpr double solar_cycle_effect = Config::solar_cycle_effect;

//! Magnetic mixing indicator variable: 0 means unipolar field, 1 means sectored field (transient)
   double Bmix_ind;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(Component comp) override;

public:

//! Default constructor
   DiffusionEmpiricalSOQLTandUNLT(void);

//! Copy constructor
   DiffusionEmpiricalSOQLTandUNLT(const DiffusionEmpiricalSOQLTandUNLT& other);

//! Destructor
   ~DiffusionEmpiricalSOQLTandUNLT() override = default;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(Component comp) override;
};

};

// Something like this is needed for templated classes
#include "diffusion_other.cc"

#endif
