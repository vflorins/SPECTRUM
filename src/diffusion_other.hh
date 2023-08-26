/*!
\file diffusion_other.hh
\brief Declares several classes to compute difusion coefficients
\author Vladimir Florinski

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

Parameters: (DiffusionBase), double A2A, double k_min, double ps_index
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

Parameters: (DiffusionBase), double D0
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
};

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

};

#endif
