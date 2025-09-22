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
// DiffusionIsotropicConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Isotropic scattering in pitch angle only, uniform in space
\author Vladimir Florinski

Parameters: (DiffusionBase), double D0
*/
template <typename Trajectory_>
class DiffusionIsotropicConstant : public DiffusionBase<Trajectory_> {

//! Readable name of the DiffusionIsotropicConstant class
   static constexpr std::string_view diff_name = "DiffusionIsotropicConstant";

public:

   using Trajectory = Trajectory_;
   using HConfig = Trajectory::HConfig;
   using DiffusionCoordinates = HConfig::DiffusionCoordinates;
   using DiffusionFields = HConfig::DiffusionFields;
   using DiffusionBase = DiffusionBase<Trajectory>;

   using DiffusionBase::_status;
   using DiffusionBase::container;
   using DiffusionBase::Kappa;
   using DiffusionBase::mu;
   using DiffusionBase::Stage;

   static_assert(!std::same_as<Trajectory, TrajectoryParker<HConfig>>, "DiffusionIsotropicConstant diffusion type cannot be applied to the Parker Trajectory type.");

protected:

//! Scattering frequency (persistent)
   double D0;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(int comp) override;

public:

//! Default constructor
   DiffusionIsotropicConstant(void);

//! Copy constructor
   DiffusionIsotropicConstant(const DiffusionIsotropicConstant& other);

//! Destructor
   ~DiffusionIsotropicConstant() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionIsotropicConstant);

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(int comp) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionQLTConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Scattering in pitch angle according to quasi-linear theory for slab turbulence, uniform in space
\author Vladimir Florinski

Parameters: (DiffusionBase), double A2A, double l_max, double ps_index
*/
template <typename Trajectory_>
class DiffusionQLTConstant : public DiffusionBase<Trajectory_> {

//! Readable name of the DiffusionQLTConstant class
   static constexpr std::string_view diff_name = "DiffusionQLTConstant";

public:

   using Trajectory = Trajectory_;
   using HConfig = Trajectory::HConfig;
   using DiffusionCoordinates = HConfig::DiffusionCoordinates;
   using DiffusionFields = HConfig::DiffusionFields;
   using DiffusionBase = DiffusionBase<Trajectory>;

//   using DiffusionBase::_status;
   using DiffusionBase::container;
   using DiffusionBase::Kappa;
   using DiffusionBase::mu;
   using DiffusionBase::Omega;
   using DiffusionBase::st2;
   using DiffusionBase::vmag;
   using DiffusionBase::Stage;

   static_assert(!std::same_as<Trajectory, TrajectoryParker<HConfig>>, "DiffusionQLTConstant diffusion type cannot be applied to the Parker Trajectory type.");

protected:

//! Alfven turbulence relative variance (persistent)
   double A2A;

//! Maximum turbulent lengthscale (persistent)
   double l_max;

//! Characteristic wavenumber (persistent)
   double k_min;

//! Power spectral index (persistent)
   double ps_index;

//! Power spectral index minus one (persistent)
   double ps_minus;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(int comp) override;

public:

//! Default constructor
   DiffusionQLTConstant(void);

//! Constructor with arguments (to speed up construction of derived classes)
   DiffusionQLTConstant(const std::string& name_in, uint16_t status_in);

//! Copy constructor
   DiffusionQLTConstant(const DiffusionQLTConstant& other);

//! Destructor
   ~DiffusionQLTConstant() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionQLTConstant);
};



//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionWNLTConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A weakly nonlinear model with pitch angle scattering and perpendicular diffusion with constant turbulent ratios.
\author Vladimir Florinski

Parameters: (DiffusionQLTConstant), double A2T, double A2L
*/
template <typename Trajectory_>
class DiffusionWNLTConstant : public DiffusionQLTConstant<Trajectory_> {

//! Readable name of the DiffusionWNLTConstant class
   static constexpr std::string_view diff_name = "DiffusionWNLTConstant";

public:

   using Trajectory = Trajectory_;
   using HConfig = Trajectory::HConfig;
   using DiffusionCoordinates = HConfig::DiffusionCoordinates;
   using DiffusionFields = HConfig::DiffusionFields;
   using HConfig::specie;
   using DiffusionBase = DiffusionBase<Trajectory>;
   using DiffusionQLTConstant = DiffusionQLTConstant<Trajectory>;

   //   using DiffusionBase::_status;
   using DiffusionBase::container;
   using DiffusionBase::Kappa;
   using DiffusionBase::mu;
   using DiffusionBase::Omega;
   using DiffusionBase::st2;
   using DiffusionBase::vmag;
   using DiffusionBase::Stage;

   using DiffusionQLTConstant::ps_index;
   using DiffusionQLTConstant::k_min;
   using DiffusionQLTConstant::ps_minus;

   static_assert(!std::same_as<Trajectory, TrajectoryParker<HConfig>>, "DiffusionWNLTConstant diffusion type cannot be applied to the Parker Trajectory type.");

protected:

//! Transverse turbulence relative variance (persistent)
   double A2T;

//! Longitudinal turbulence relative variance (persistent)
   double A2L;

//! Power spectral index plus one (persistent)
   double ps_plus;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(int comp) override;

public:

//! Default constructor
   DiffusionWNLTConstant(void);

//! Constructor with arguments (to speed up construction of derived classes)
   DiffusionWNLTConstant(const std::string& name_in, uint16_t status_in);

//! Copy constructor
   DiffusionWNLTConstant(const DiffusionWNLTConstant& other);

//! Destructor
   ~DiffusionWNLTConstant() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionWNLTConstant);
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
template <typename Trajectory_>
class DiffusionWNLTRampVLISM : public DiffusionWNLTConstant<Trajectory_> {

//! Readable name of the DiffusionWNLTRampVLISM class
   static constexpr std::string_view diff_name = "DiffusionWNLTRampVLISM";

public:

   using Trajectory = Trajectory_;
   using HConfig = Trajectory::HConfig;
   using DiffusionCoordinates = HConfig::DiffusionCoordinates;
   using DiffusionFields = HConfig::DiffusionFields;
   using DiffusionBase = DiffusionBase<Trajectory>;
   using DiffusionWNLTConstant = DiffusionWNLTConstant<Trajectory>;

   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::Stage;

   using DiffusionWNLTConstant::k_min;
   using DiffusionWNLTConstant::ps_minus;
   using DiffusionWNLTConstant::A2A;
   using DiffusionWNLTConstant::A2T;
   using DiffusionWNLTConstant::A2L;
   using DiffusionWNLTConstant::l_max;

   static_assert(!std::same_as<Trajectory, TrajectoryParker<HConfig>>, "DiffusionWNLTRampVLISM diffusion type cannot be applied to the Parker Trajectory type.");

protected:

//! Reference characteristic wavenumber (persistent) (k_min becomes transient)
   double k_min_ref;

//! Alfven turbulence relative variance (persistent) (A2A becomes transient)
   double A2A_ref;

//! Transverse turbulence relative variance (persistent) (A2T becomes transient)
   double A2T_ref;

//! Longitudinal turbulence relative variance (persistent) (A2L becomes transient)
   double A2L_ref;

//! Largest scale of turbulence at HP (persistent)
   double l_max_HP;

//! Difference between largest scales (persistent)
   double dl_max;

//! Extent of the HP in the nose direction (persistent)
   double z_nose;

//! Extent of the sheath in the nose direction (persistent)
   double z_sheath;

//! Difference between nose distance (persistent)
   double dz;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(int comp) override;

public:

//! Default constructor
   DiffusionWNLTRampVLISM(void);

//! Copy constructor
   DiffusionWNLTRampVLISM(const DiffusionWNLTRampVLISM& other);

//! Destructor
   ~DiffusionWNLTRampVLISM() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionWNLTRampVLISM);
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
template <typename Trajectory_>
class DiffusionParaConstant : public DiffusionBase<Trajectory_> {

//! Readable name of the DiffusionParaConstant class
   static constexpr std::string_view diff_name = "DiffusionParaConstant";

public:

   using Trajectory = Trajectory_;
   using HConfig = Trajectory::HConfig;
   using DiffusionCoordinates = HConfig::DiffusionCoordinates;
   using DiffusionFields = HConfig::DiffusionFields;
   using DiffusionBase = DiffusionBase<Trajectory>;

   using DiffusionBase::_status;
   using DiffusionBase::container;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Diffusion coefficient (persistent)
   double D0;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(int comp) override;

public:

//! Default constructor
   DiffusionParaConstant(void);

//! Copy constructor
   DiffusionParaConstant(const DiffusionParaConstant& other);

//! Destructor
   ~DiffusionParaConstant() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionParaConstant);

//! Compute derivative of diffusion coefficient in position or time
   double GetDirectionalDerivative(int comp, int xyz, const DerivativeData& ddata) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(int comp) override;
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
template <typename Trajectory_>
class DiffusionPerpConstant : public DiffusionBase<Trajectory_> {

//! Readable name of the DiffusionPerpConstant class
   static constexpr std::string_view diff_name = "DiffusionPerpConstant";

public:

   using Trajectory = Trajectory_;
   using HConfig = Trajectory::HConfig;
   using DiffusionCoordinates = HConfig::DiffusionCoordinates;
   using DiffusionFields = HConfig::DiffusionFields;
   using DiffusionBase = DiffusionBase<Trajectory>;

   using DiffusionBase::_status;
   using DiffusionBase::container;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Diffusion coefficient (persistent)
   double D0;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(int comp) override;

public:

//! Default constructor
   DiffusionPerpConstant(void);

//! Copy constructor
   DiffusionPerpConstant(const DiffusionPerpConstant& other);

//! Destructor
   ~DiffusionPerpConstant() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionPerpConstant);

//! Compute derivative of diffusion coefficient in position or time
   double GetDirectionalDerivative(int comp, int xyz, const DerivativeData& ddata) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(int comp) override;
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
template <typename Trajectory_>
class DiffusionFullConstant : public DiffusionBase<Trajectory_> {

//! Readable name of the DiffusionFullConstant class
   static constexpr std::string_view diff_name = "DiffusionFullConstant";

public:

   using Trajectory = Trajectory_;
   using HConfig = Trajectory::HConfig;
   using DiffusionCoordinates = HConfig::DiffusionCoordinates;
   using DiffusionFields = HConfig::DiffusionFields;
   using DiffusionBase = DiffusionBase<Trajectory>;

   using DiffusionBase::_status;
   using DiffusionBase::container;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Perpendicular diffusion coefficient (persistent)
   double Dperp;

//! Parallel diffusion coefficient (persistent)
   double Dpara;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(int comp) override;

public:

//! Default constructor
   DiffusionFullConstant(void);

//! Copy constructor
   DiffusionFullConstant(const DiffusionFullConstant& other);

//! Destructor
   ~DiffusionFullConstant() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionFullConstant);

//! Compute derivative of diffusion coefficient in position or time
   double GetDirectionalDerivative(int comp, int xyz, const DerivativeData& ddata) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(int comp) override;
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
template <typename Trajectory_>
class DiffusionFlowMomentumPowerLaw : public DiffusionBase<Trajectory_> {

//! Readable name of the DiffusionFlowMomentumPowerLaw class
   static constexpr std::string_view diff_name = "DiffusionFlowMomentumPowerLaw";


public:

   using Trajectory = Trajectory_;
   using HConfig = Trajectory::HConfig;
   using DiffusionCoordinates = HConfig::DiffusionCoordinates;
   using DiffusionFields = HConfig::DiffusionFields;
   using DiffusionBase = DiffusionBase<Trajectory>;

   using DiffusionBase::_status;
   using DiffusionBase::_mom;
   using DiffusionBase::_fields;
   using DiffusionBase::container;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Reference diffusion coefficient (persistent)
   double kap0;

//! Flow velocity normalization factor (persistent)
   double U0;

//! Power law slope for flow velocity (persistent)
   double pow_law_U;
   
   //! Momentum normalization factor (persistent)
   double p0;

//! Power law slope for momentum (persistent)
   double pow_law_p;

//! Ratio of perpendicular to parallel diffusion (persistent)
   double kap_rat;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(int comp) override;

public:

//! Default constructor
   DiffusionFlowMomentumPowerLaw(void);

//! Copy constructor
   DiffusionFlowMomentumPowerLaw(const DiffusionFlowMomentumPowerLaw& other);

//! Destructor
   ~DiffusionFlowMomentumPowerLaw() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionFlowMomentumPowerLaw);

//! Compute derivative of diffusion coefficient in position or time
   double GetDirectionalDerivative(int comp, int xyz, const DerivativeData& ddata) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(int comp) override;
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
template <typename Trajectory_>
class DiffusionKineticEnergyRadialDistancePowerLaw : public DiffusionBase<Trajectory_> {


//! Readable name of the DiffusionKineticEnergyRadialDistancePowerLaw class
   static constexpr std::string_view diff_name = "DiffusionKineticEnergyRadialDistancePowerLaw";

public:

   using Trajectory = Trajectory_;
   using HConfig = Trajectory::HConfig;
   using DiffusionCoordinates = HConfig::DiffusionCoordinates;
   using DiffusionFields = HConfig::DiffusionFields;
   using DiffusionBase = DiffusionBase<Trajectory>;

   using DiffusionBase::_status;
   using DiffusionBase::specie;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Diffusion coefficient normalization factor (persistent)
   double kap0;

//! Kinetic Energy normalization factor (persistent)
   double T0;

//! Radial distance normalization factor (persistent)
   double r0;

//! Power law slope for kinetic energy (persistent)
   double pow_law_T;

//! Power law slope for radial distance (persistent)
   double pow_law_r;

//! Ratio of perpendicular to parallel diffusion (persistent)
   double kap_rat;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(int comp) override;

public:

//! Default constructor
   DiffusionKineticEnergyRadialDistancePowerLaw(void);

//! Copy constructor
   DiffusionKineticEnergyRadialDistancePowerLaw(const DiffusionKineticEnergyRadialDistancePowerLaw& other);

//! Destructor
   ~DiffusionKineticEnergyRadialDistancePowerLaw() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionKineticEnergyRadialDistancePowerLaw);

//! Compute derivative of diffusion coefficient in position or time
   double GetDirectionalDerivative(int comp, int xyz, const DerivativeData& ddata) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(int comp) override;
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
template <typename Trajectory_>
class DiffusionRigidityMagneticFieldPowerLaw : public DiffusionBase<Trajectory_> {

//! Readable name of the DiffusionRigidityMagneticFieldPowerLaw class
   static constexpr std::string_view diff_name = "DiffusionRigidityMagneticFieldPowerLaw";

public:

   using Trajectory = Trajectory_;
   using HConfig = Trajectory::HConfig;
   using DiffusionCoordinates = HConfig::DiffusionCoordinates;
   using DiffusionFields = HConfig::DiffusionFields;
   using DiffusionBase = DiffusionBase<Trajectory>;

   using DiffusionBase::_status;
   using DiffusionBase::specie;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::vmag;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Parallel mean free path (persistent)
   double lam0;

//! Rigidity normalization factor (persistent)
   double R0;

//! Magnetic field normalization factor (persistent)
   double B0;

//! Power law slope for rigidity (persistent)
   double pow_law_R;

//! Power law slope for magnetic field (persistent)
   double pow_law_B;

//! Ratio of perpendicular to parallel diffusion (persistent)
   double kap_rat;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(int comp) override;

public:

//! Default constructor
   DiffusionRigidityMagneticFieldPowerLaw(void);

//! Copy constructor
   DiffusionRigidityMagneticFieldPowerLaw(const DiffusionRigidityMagneticFieldPowerLaw& other);

//! Destructor
   ~DiffusionRigidityMagneticFieldPowerLaw() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionRigidityMagneticFieldPowerLaw);

//! Compute derivative of diffusion coefficient in position or time
   double GetDirectionalDerivative(int comp, int xyz, const DerivativeData& ddata) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(int comp) override;
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
template <typename Trajectory_>
class DiffusionStraussEtAl2013 : public DiffusionBase<Trajectory_> {

//! Readable name of the DiffusionStraussEtAl2013 class
   static constexpr std::string_view diff_name = "DiffusionStraussEtAl2013";

public:

   using Trajectory = Trajectory_;
   using HConfig = Trajectory::HConfig;
   using DiffusionCoordinates = HConfig::DiffusionCoordinates;
   using DiffusionFields = HConfig::DiffusionFields;
   using DiffusionBase = DiffusionBase<Trajectory>;

   using DiffusionBase::_status;
   using DiffusionBase::specie;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::vmag;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Index for LISM indicator variable (persistent)
   int LISM_idx;

//! Parallel inner heliosphere mean free path (persistent)
   double lam_in;

//! Parallel outer heliosphere mean free path (persistent)
   double lam_out;

//! Rigidity normalization factor (persistent)
   double R0;

//! Magnetic field normalization factor for inner heliosphere (persistent)
   double B0;

//! Ratio of perpendicular to parallel diffusion inner heliosphere (persistent)
   double kap_rat_in;

//! Ratio of perpendicular to parallel diffusion outer heliosphere (persistent)
   double kap_rat_out;

//! LISM indicator variable: 0 means inside HP, 1 means outside HP (transient)
   double LISM_ind;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(int comp) override;

public:

//! Default constructor
   DiffusionStraussEtAl2013(void);

//! Copy constructor
   DiffusionStraussEtAl2013(const DiffusionStraussEtAl2013& other);

//! Destructor
   ~DiffusionStraussEtAl2013() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionStraussEtAl2013);

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(int comp) override;
};

template <typename Trajectory>
using DiffusionGuoEtAl2014 = DiffusionStraussEtAl2013<Trajectory>;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionPotgieterEtAl2015 class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Full (perpendicular + parallel) diffusion, rigidity and magnetic field power law according to Potgieter et al 2015
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionBase), int LISM_idx, double kappa_in, double kappa_out, double R0, double B0, double kap_rat_in, double kap_rat_out
*/
template <typename Trajectory_>
class DiffusionPotgieterEtAl2015 : public DiffusionBase<Trajectory_> {

//! Readable name of the DiffusionPotgieterEtAl2015 class
   static constexpr std::string_view diff_name = "DiffusionPotgieterEtAl2015";

public:

   using Trajectory = Trajectory_;
   using HConfig = Trajectory::HConfig;
   using DiffusionCoordinates = HConfig::DiffusionCoordinates;
   using DiffusionFields = HConfig::DiffusionFields;
   using DiffusionBase = DiffusionBase<Trajectory>;

   using DiffusionBase::_status;
   using DiffusionBase::specie;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::vmag;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Index for LISM indicator variable (persistent)
   int LISM_idx;

//! Parallel inner heliosphere diffusion coefficient (persistent)
   double kappa_in;

//! Parallel outer heliosphere diffusion coefficient (persistent)
   double kappa_out;

//! Rigidity normalization factor (persistent)
   double R0;

//! Magnetic field normalization factor for inner heliosphere (persistent)
   double B0;

//! Ratio of perpendicular to parallel diffusion inner heliosphere (persistent)
   double kap_rat_in;

//! Ratio of perpendicular to parallel diffusion outer heliosphere (persistent)
   double kap_rat_out;

//! LISM indicator variable: 0 means inside HP, 1 means outside HP (transient)
   double LISM_ind;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(int comp) override;

public:

//! Default constructor
   DiffusionPotgieterEtAl2015(void);

//! Copy constructor
   DiffusionPotgieterEtAl2015(const DiffusionPotgieterEtAl2015& other);

//! Destructor
   ~DiffusionPotgieterEtAl2015() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionPotgieterEtAl2015);

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(int comp) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionEmpiricalSOQLTandUNLT class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Full (perpendicular + parallel) diffusion, parametrized empirical fit to SOQLT (parallel) and UNLT (perpendicular) using a bent power-law spectrum with change in perpendicular diffusion according to magnetic mixing and solar cycle indicator variables
\author Juan G Alonso Guzman
Parameters: (DiffusionBase), double lam_para, double lam_perp, double R0, double B0, int Bmix_idx, double kap_rat_red, double radial_limit_perp_red, int solar_cycle_idx, double solar_cycle_effect
*/
template <typename Trajectory_>
class DiffusionEmpiricalSOQLTandUNLT : public DiffusionBase<Trajectory_> {

//! Readable name of the DiffusionEmpiricalSOQLTandUNLT class
   static constexpr std::string_view diff_name = "DiffusionEmpiricalSOQLTandUNLT";

public:

   using Trajectory = Trajectory_;
   using HConfig = Trajectory::HConfig;
   using DiffusionCoordinates = HConfig::DiffusionCoordinates;
   using DiffusionFields = HConfig::DiffusionFields;
   using DiffusionBase = DiffusionBase<Trajectory>;

   using DiffusionBase::_status;
   using DiffusionBase::specie;
   using DiffusionBase::container;
   using DiffusionBase::_coords;
   using DiffusionBase::_fields;
   using DiffusionBase::vmag;
   using DiffusionBase::Kappa;
   using DiffusionBase::Stage;

protected:

//! Parallel mean free path (persistent)
   double lam_para;

//! Perpendicular mean free path (persistent)
   double lam_perp;

//! Rigidity normalization factor (persistent)
   double R0;

//! Magnetic field normalization factor for inner heliosphere (persistent)
   double B0;

//! Index for magnetic mixing indicator variable (persistent)
   int Bmix_idx;

//! Reduction factor for kappa in unipolar regions (persistent)
   double kap_rat_red;

//! Radial limit to apply unipolar reduction factor (persistent)
   double radial_limit_perp_red;

//! Index for solar cycle indicator variable (persistent)
   int solar_cycle_idx;
   
//! Solar cycle effect constant (persistent)
   double solar_cycle_effect;

//! Magnetic mixing indicator variable: 0 means unipolar field, 1 means sectored field (transient)
   double Bmix_ind;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(int comp) override;

public:

//! Default constructor
   DiffusionEmpiricalSOQLTandUNLT(void);

//! Copy constructor
   DiffusionEmpiricalSOQLTandUNLT(const DiffusionEmpiricalSOQLTandUNLT& other);

//! Destructor
   ~DiffusionEmpiricalSOQLTandUNLT() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionEmpiricalSOQLTandUNLT);

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(int comp) override;
};

};

// Something like this is needed for templated classes
#include "diffusion_other.cc"

#endif
