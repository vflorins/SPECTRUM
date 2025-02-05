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

#if TRAJ_TYPE != TRAJ_PARKER

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionIsotropicConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DiffusionIsotropicConstant class
const std::string diff_name_isotropic_constant = "DiffusionIsotropicConstant";

/*!
\brief Isotropic scattering in pitch angle only, uniform in space
\author Vladimir Florinski

Parameters: (DiffusionBase), double D0
*/
class DiffusionIsotropicConstant : public DiffusionBase {

protected:

//! Scattering frequency (persistent)
   double D0;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(void) override;

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
   double GetMuDerivative(void) override;
};

#endif

#if TRAJ_TYPE != TRAJ_PARKER

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionQLTConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DiffusionQLTConstant class
const std::string diff_name_qlt_constant = "DiffusionQLTConstant";

/*!
\brief Scattering in pitch angle according to quasi-linear theory for slab turbulence, uniform in space
\author Vladimir Florinski

Parameters: (DiffusionBase), double A2A, double l_max, double ps_index
*/
class DiffusionQLTConstant : public DiffusionBase {

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
   void EvaluateDiffusion(void) override;

public:

//! Default constructor
   DiffusionQLTConstant(void);

//! Constructor with arguments (to speed up construction of derived classes)
   DiffusionQLTConstant(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor
   DiffusionQLTConstant(const DiffusionQLTConstant& other);

//! Destructor
   ~DiffusionQLTConstant() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionQLTConstant);
};

#endif

#if TRAJ_TYPE != TRAJ_PARKER

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionWNLTConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DiffusionWNLTConstant class
const std::string diff_name_wnlt_constant = "DiffusionWNLTConstant";

//! Flag to use QLT pitch angle scattering with WLNT perpendicular diffusion
#define USE_QLT_SCATT_WITH_WNLT_DIFF

/*!
\brief A weakly nonlinear model with pitch angle scattering and perpendicular diffusion with constant turbulent ratios.
\author Vladimir Florinski

Parameters: (DiffusionQLTConstant), double A2T, double A2L
*/
class DiffusionWNLTConstant : public DiffusionQLTConstant {

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
   void EvaluateDiffusion(void) override;

public:

//! Default constructor
   DiffusionWNLTConstant(void);

//! Constructor with arguments (to speed up construction of derived classes)
   DiffusionWNLTConstant(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

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

//! Readable name of the DiffusionWNLTRampVLISM class
const std::string diff_name_wnlt_ramp_vlism = "DiffusionWNLTRampVLISM";

/*!
\brief WNLT diffusion in VLISM where the largest scale to turbulence increases linearly (ramp) between the HP (Rankine half body) and the sheath surface (also a Rankine half body) and is constant elsewhere.
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionWNLTConstant)
*/
class DiffusionWNLTRampVLISM : public DiffusionWNLTConstant {

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
   void EvaluateDiffusion(void) override;

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

#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionParaConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DiffusionParaConstant class
const std::string diff_name_para_constant = "DiffusionParaConstant";

/*!
\brief Parallel diffusion uniform in space
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionBase), double D0
*/
class DiffusionParaConstant : public DiffusionBase {

protected:

//! Diffusion coefficient (persistent)
   double D0;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(void) override;

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
   double GetDirectionalDerivative(int xyz) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionPerpConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DiffusionPerpConstant class
const std::string diff_name_perp_constant = "DiffusionPerpConstant";

/*!
\brief Perpendicular diffusion uniform in space
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionBase), double D0
*/
class DiffusionPerpConstant : public DiffusionBase {

protected:

//! Diffusion coefficient (persistent)
   double D0;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(void) override;

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
   double GetDirectionalDerivative(int xyz) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionFullConstant class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DiffusionFullConstant class
const std::string diff_name_full_constant = "DiffusionFullConstant";

/*!
\brief Full (perpendicular + parallel) diffusion uniform in space
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionBase), double Dperp, double Dpara
*/
class DiffusionFullConstant : public DiffusionBase {

protected:

//! Perpendicular diffusion coefficient (persistent)
   double Dperp;

//! Parallel diffusion coefficient (persistent)
   double Dpara;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(void) override;

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
   double GetDirectionalDerivative(int xyz) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionFlowPowerLaw class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DiffusionFlowPowerLaw class
const std::string diff_name_flow_power_law = "DiffusionFlowPowerLaw";

/*!
\brief Diffusion as a power law of flow velocity magnitude
\author Juan G Alonso Guzman
\author Swati Sharma

Parameters: (DiffusionBase), double kappa0, double u0, double power_law_U, double kap_rat
*/
class DiffusionFlowPowerLaw : public DiffusionBase {

protected:

//! Reference diffusion coefficient (persistent)
   double kappa0;

//! Flow velocity normalization factor (persistent)
   double U0;

//! Power law slope for flow velocity (persistent)
   double pow_law_U;

//! Ratio of perpendicular to parallel diffusion (persistent)
   double kap_rat;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(void) override;

public:

//! Default constructor
   DiffusionFlowPowerLaw(void);

//! Copy constructor
   DiffusionFlowPowerLaw(const DiffusionFlowPowerLaw& other);

//! Destructor
   ~DiffusionFlowPowerLaw() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionFlowPowerLaw);

//! Compute derivative of diffusion coefficient in position or time
   double GetDirectionalDerivative(int xyz) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionMomentumPowerLaw class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DiffusionMomentumPowerLaw class
const std::string diff_name_momentum_power_law = "DiffusionMomentumPowerLaw";

/*!
\brief Full (perpendicular + parallel) diffusion, momentum (magnitude) power law
\author Juan G Alonso Guzman

Parameters: (DiffusionBase), kappa0, double p0, double pow_law_p, double kap_rat
*/
class DiffusionMomentumPowerLaw : public DiffusionBase {

protected:

//! Reference diffusion coefficient (persistent)
   double kappa0;

//! Momentum normalization factor (persistent)
   double p0;

//! Power law slope for momentum (persistent)
   double pow_law_p;

//! Ratio of perpendicular to parallel diffusion (persistent)
   double kap_rat;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(void) override;

public:

//! Default constructor
   DiffusionMomentumPowerLaw(void);

//! Copy constructor
   DiffusionMomentumPowerLaw(const DiffusionMomentumPowerLaw& other);

//! Destructor
   ~DiffusionMomentumPowerLaw() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionMomentumPowerLaw);

//! Compute derivative of diffusion coefficient in position or time
   double GetDirectionalDerivative(int xyz) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionKineticEnergyRadialDistancePowerLaw class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DiffusionKineticEnergyRadialDistancePowerLaw class
const std::string diff_name_kinetic_energy_radial_distance_power_law = "DiffusionKineticEnergyRadialDistancePowerLaw";

/*!
\brief Full (perpendicular + parallel) diffusion, kinetic energy and radial distance power law
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionBase), kap0, double T0, double r0, double pow_law_T, double pow_law_r, double kap_rat
*/
class DiffusionKineticEnergyRadialDistancePowerLaw : public DiffusionBase {

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
   void EvaluateDiffusion(void) override;

public:

//! Default constructor
   DiffusionKineticEnergyRadialDistancePowerLaw(void);

//! Copy constructor
   DiffusionKineticEnergyRadialDistancePowerLaw(const DiffusionKineticEnergyRadialDistancePowerLaw& other);

//! Destructor
   ~DiffusionKineticEnergyRadialDistancePowerLaw() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionKineticEnergyRadialDistancePowerLaw);

// TODO: Implement directional derivative (only spatial dependence is in radial distance)
//! Compute derivative of diffusion coefficient in position or time
   // double GetDirectionalDerivative(int xyz) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionRigidityMagneticFieldPowerLaw class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DiffusionRigidityMagneticFieldPowerLaw class
const std::string diff_name_rigidity_magnetic_field_power_law = "DiffusionRigidityMagneticFieldPowerLaw";

/*!
\brief Full (perpendicular + parallel) diffusion, rigidity and magnetic field power law
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionBase), lam0, double R0, double B0, double pow_law_R, double pow_law_B, double kap_rat
*/
class DiffusionRigidityMagneticFieldPowerLaw : public DiffusionBase {

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
   void EvaluateDiffusion(void) override;

public:

//! Default constructor
   DiffusionRigidityMagneticFieldPowerLaw(void);

//! Copy constructor
   DiffusionRigidityMagneticFieldPowerLaw(const DiffusionRigidityMagneticFieldPowerLaw& other);

//! Destructor
   ~DiffusionRigidityMagneticFieldPowerLaw() override = default;

//! Clone function
   CloneFunctionDiffusion(DiffusionRigidityMagneticFieldPowerLaw);

// TODO: Implement directional derivative (only spatial dependence is in Bmag)
//! Compute derivative of diffusion coefficient in position or time
   // double GetDirectionalDerivative(int xyz) override;

//! Compute derivative of diffusion coefficient in mu
   double GetMuDerivative(void) override;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionStraussEtAl2013 class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DiffusionStraussEtAl2013 class
const std::string diff_name_strauss_et_al_2013 = "DiffusionStraussEtAl2013";

/*!
\brief Full (perpendicular + parallel) diffusion, rigidity and magnetic field power law according to Strauss et al 2013, with change in perpendicular diffusion according to a magnetic mixing indicator variable
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionBase), lam_in, lam_out, double R0, double B0, double pow_law_R_low, double pow_law_R_high, double pow_law_B, double kap_rat_low, double kap_rat_high
*/
class DiffusionStraussEtAl2013 : public DiffusionBase {

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

//! Index for magnetic mixing indicator variable (persistent)
   int Bmix_idx;

//! Reduction factor for kappa in unipolar regions (persistent)
   double kap_rat_red;

//! LISM indicator variable: 0 means inside HP, 1 means outside HP (transient)
   double LISM_ind;

//! Magnetic mixing indicator variable: 0 means unipolar field, 1 means sectored field (transient)
   double Bmix_ind;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(void) override;

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
   double GetMuDerivative(void) override;
};

typedef DiffusionStraussEtAl2013 DiffusionGuoEtAl2014;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionPotgieterEtAl2015 class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the DiffusionPotgieterEtAl2015 class
const std::string diff_name_potgieter_et_al_2015 = "DiffusionPotgieterEtAl2015";

/*!
\brief Full (perpendicular + parallel) diffusion, rigidity and magnetic field power law according to Potgieter et al 2015, with change in perpendicular diffusion according to a magnetic mixing indicator variable
\author Juan G Alonso Guzman
\author Vladimir Florinski

Parameters: (DiffusionBase), lam_in, lam_out, double R0, double B0, double pow_law_R_low, double pow_law_R_high, double pow_law_B, double kap_rat_low, double kap_rat_high
*/
class DiffusionPotgieterEtAl2015 : public DiffusionBase {

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

//! Index for magnetic mixing indicator variable (persistent)
   int Bmix_idx;

//! Reduction factor for kappa in unipolar regions (persistent)
   double kap_rat_red;

//! LISM indicator variable: 0 means inside HP, 1 means outside HP (transient)
   double LISM_ind;

//! Magnetic mixing indicator variable: 0 means unipolar field, 1 means sectored field (transient)
   double Bmix_ind;

//! Set up the diffusion model based on "params"
   void SetupDiffusion(bool construct) override;

//! Compute the diffusion coefficients
   void EvaluateDiffusion(void) override;

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
   double GetMuDerivative(void) override;
};

};

#endif
