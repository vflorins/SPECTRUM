/*!
\file diffusion_other.cc
\brief Implements several classes to compute difusion coefficients
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "diffusion_other.hh"
#include <gsl/gsl_sf_hyperg.h>

namespace Spectrum {

#if TRAJ_TYPE != TRAJ_PARKER

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionIsotropicConstant methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/06/2022
*/
DiffusionIsotropicConstant::DiffusionIsotropicConstant(void)
                          : DiffusionBase(diff_name_isotropic_constant, 0, DIFF_NOBACKGROUND)
{
};

/*!
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
DiffusionIsotropicConstant::DiffusionIsotropicConstant(const DiffusionIsotropicConstant& other)
                          : DiffusionBase(other)
{
   RAISE_BITS(_status, DIFF_NOBACKGROUND);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Vladimir Florinski
\date 05/06/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void DiffusionIsotropicConstant::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(D0);
};

/*!
\author Vladimir Florinski
\date 05/06/2022
*/
void DiffusionIsotropicConstant::EvaluateDiffusion(void)
{
   if ((comp_eval == 0) || (comp_eval == 1)) return;
   Kappa[2] = D0 * (1.0 - Sqr(mu));
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinsk
\date 10/19/2022
\return double       Derivative in mu
*/
double DiffusionIsotropicConstant::GetMuDerivative(void)
{
   return -2.0 * D0 * mu;
};

#endif

#if TRAJ_TYPE != TRAJ_PARKER

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionQLTConstant methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/09/2022
*/
DiffusionQLTConstant::DiffusionQLTConstant(void)
                    : DiffusionBase(diff_name_qlt_constant, 0, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 05/27/2022
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
DiffusionQLTConstant::DiffusionQLTConstant(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                    : DiffusionBase(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
DiffusionQLTConstant::DiffusionQLTConstant(const DiffusionQLTConstant& other)
                    : DiffusionBase(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Vladimir Florinski
\date 05/09/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void DiffusionQLTConstant::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(A2A);
   container.Read(l_max);
   container.Read(ps_index);
   k_min = M_2PI / l_max;
   ps_minus = ps_index - 1.0;
};

/*!
\author Vladimir Florinski
\date 05/09/2022
*/
void DiffusionQLTConstant::EvaluateDiffusion(void)
{
   if ((comp_eval == 0) || (comp_eval == 1)) return;
   Kappa[2] = 0.25 * M_PI * ps_minus * fabs(Omega) * st2 * pow(vmag * k_min * fabs(mu / Omega), ps_minus) * A2A;
};

#endif

#if TRAJ_TYPE != TRAJ_PARKER

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionWNLTConstant methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/06/2022
*/
DiffusionWNLTConstant::DiffusionWNLTConstant(void)
                     : DiffusionQLTConstant(diff_name_wnlt_constant, 0, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 05/27/2022
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
DiffusionWNLTConstant::DiffusionWNLTConstant(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                     : DiffusionQLTConstant(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
DiffusionWNLTConstant::DiffusionWNLTConstant(const DiffusionWNLTConstant& other)
                     : DiffusionQLTConstant(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Vladimir Florinski
\date 05/09/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void DiffusionWNLTConstant::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionQLTConstant::SetupDiffusion(false);
   container.Read(A2T);
   container.Read(A2L);
   ps_plus = ps_index + 1.0;
};

/*!
\author Vladimir Florinski
\date 05/09/2022
*/
void DiffusionWNLTConstant::EvaluateDiffusion(void)
{
   if (comp_eval == 1) return;
   double CT, CL, xi1, xi2, F21, DT1, DT2;

// Kappa[0] is required for Kappa[2]
   CT = 0.5 * ps_minus / ps_plus * Sqr(mu / k_min) * A2T;
   CL = 0.125 * Sqr(vmag * st2 / Omega) * A2L;
   Kappa[0] = vmag * sqrt(CT + CL);

   if (comp_eval == 0) return;

// Hypergeometric function may crash if the last argument is close to 1
   DiffusionQLTConstant::EvaluateDiffusion();

#ifdef USE_QLT_SCATT_WITH_WNLT_DIFF
   return;
#endif

   if (A2T > sp_tiny) {
      xi1 = vmag * k_min * sqrt(st2) / M_SQRT2 / fabs(Omega);
      F21 = gsl_sf_hyperg_2F1(1.0, 1.0, (5.0 + ps_index) / 4.0, 1.0 / (1.0 + Sqr(xi1)));
      DT1 = 1.0 / (1.0 + Sqr(xi1)) * F21;
      xi2 = Kappa[0] * k_min * k_min / fabs(Omega);
      F21 = gsl_sf_hyperg_2F1(1.0, 1.0, (5.0 + ps_index) / 4.0, 1.0 / (1.0 + Sqr(xi2)));
      DT2 = 1.0 / (1.0 + Sqr(xi2)) * F21;
      Kappa[2] += (ps_minus / ps_plus) * Omega * st2 * A2T * DT1 * DT2 / sqrt(Sqr(DT1) + Sqr(DT2));
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionWNLTRampVLISM methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/06/2022
*/
DiffusionWNLTRampVLISM::DiffusionWNLTRampVLISM(void)
                      : DiffusionWNLTConstant(diff_name_wnlt_ramp_vlism, 0, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
DiffusionWNLTRampVLISM::DiffusionWNLTRampVLISM(const DiffusionWNLTRampVLISM& other)
                      : DiffusionWNLTConstant(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/09/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void DiffusionWNLTRampVLISM::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionWNLTConstant::SetupDiffusion(false);
   container.Read(l_max_HP);
   container.Read(z_nose);
   container.Read(z_sheath);
   k_min_ref = k_min;
   A2A_ref = A2A;
   A2T_ref = A2T;
   A2L_ref = A2L;
   dl_max = l_max - l_max_HP;
   dz = z_sheath - z_nose;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/09/2022
*/
void DiffusionWNLTRampVLISM::EvaluateDiffusion(void)
{
   if (comp_eval == 1) return;
   double z0, r0, k_ratio;

// Find distance to nose for Rankine half body which the particle is presently on
   r0 = _pos.Norm();
   z0 = sqrt(0.5 * r0 * (r0 + _pos[2]));

// Constant k_min beyond z_sheath
   if (z0 > z_sheath) {
      k_min = k_min_ref;
      A2A = A2A_ref;
      A2T = A2T_ref;
      A2L = A2L_ref;
   }
// Linearly interpolate l_max between z_nose and z_sheath
   else {
      k_min = M_2PI / (l_max_HP + dl_max * (z0 - z_nose) / dz);
      k_ratio = pow(k_min_ref / k_min, ps_minus);
      A2A = A2A_ref * k_ratio;
      A2T = A2T_ref * k_ratio;
      A2L = A2L_ref * k_ratio;
   };

// Evaluate WLNT diffusion
   DiffusionWNLTConstant::EvaluateDiffusion();
};

#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionParaConstant methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/07/2023
*/
DiffusionParaConstant::DiffusionParaConstant(void)
                     : DiffusionBase(diff_name_para_constant, 0, DIFF_NOBACKGROUND)
{
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/07/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
DiffusionParaConstant::DiffusionParaConstant(const DiffusionParaConstant& other)
                     : DiffusionBase(other)
{
   RAISE_BITS(_status, DIFF_NOBACKGROUND);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/07/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void DiffusionParaConstant::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(D0);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/07/2023
*/
void DiffusionParaConstant::EvaluateDiffusion(void)
{
   if ((comp_eval == 0) || (comp_eval == 2)) return;
   Kappa[1] = D0;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/07/2023
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\return double       Directional derivative
\note This is meant to be called after GetComponent() for the component for which the derivative is wanted
*/
double DiffusionParaConstant::GetDirectionalDerivative(int xyz)
{
   return 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
double DiffusionParaConstant::GetMuDerivative(void)
{
   return 0.0;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionPerpConstant methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/06/2022
*/
DiffusionPerpConstant::DiffusionPerpConstant(void)
                     : DiffusionBase(diff_name_perp_constant, 0, DIFF_NOBACKGROUND)
{
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
DiffusionPerpConstant::DiffusionPerpConstant(const DiffusionPerpConstant& other)
                     : DiffusionBase(other)
{
   RAISE_BITS(_status, DIFF_NOBACKGROUND);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/06/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void DiffusionPerpConstant::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(D0);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/06/2022
*/
void DiffusionPerpConstant::EvaluateDiffusion(void)
{
   if ((comp_eval == 1) || (comp_eval == 2)) return;
   Kappa[0] = D0;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/18/2022
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\return double       Directional derivative
\note This is meant to be called after GetComponent() for the component for which the derivative is wanted
*/
double DiffusionPerpConstant::GetDirectionalDerivative(int xyz)
{
   return 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
double DiffusionPerpConstant::GetMuDerivative(void)
{
   return 0.0;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionFullConstant methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/06/2022
*/
DiffusionFullConstant::DiffusionFullConstant(void)
                     : DiffusionBase(diff_name_full_constant, 0, DIFF_NOBACKGROUND)
{
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
DiffusionFullConstant::DiffusionFullConstant(const DiffusionFullConstant& other)
                     : DiffusionBase(other)
{
   RAISE_BITS(_status, DIFF_NOBACKGROUND);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/06/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void DiffusionFullConstant::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(Dperp);
   container.Read(Dpara);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/06/2022
*/
void DiffusionFullConstant::EvaluateDiffusion(void)
{
   if (comp_eval == 2) return;
   Kappa[0] = Dperp;
   Kappa[1] = Dpara;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/18/2022
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\return double       Directional derivative
\note This is meant to be called after GetComponent() for the component for which the derivative is wanted
*/
double DiffusionFullConstant::GetDirectionalDerivative(int xyz)
{
   return 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
double DiffusionFullConstant::GetMuDerivative(void)
{
   return 0.0;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionFlowPowerLaw methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Swati Sharma
\date 06/03/2024
*/
DiffusionFlowPowerLaw::DiffusionFlowPowerLaw(void)
                     : DiffusionBase(diff_name_flow_power_law, 0, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\author Swati Sharma
\date 06/03/2024
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
DiffusionFlowPowerLaw::DiffusionFlowPowerLaw(const DiffusionFlowPowerLaw& other)
                     : DiffusionBase(other)
{
   RAISE_BITS(_status, STATE_NONE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Juan G Alonso Guzman
\author Swati Sharma
\date 06/03/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void DiffusionFlowPowerLaw::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(kappa0);
   container.Read(U0);
   container.Read(pow_law_U);
   container.Read(kap_rat);
};

/*!
\author Juan G Alonso Guzman
\author Swati Sharma
\date 06/03/2024
*/
void DiffusionFlowPowerLaw::EvaluateDiffusion(void)
{
   if (comp_eval == 2) return;
   Kappa[1] = kappa0 * pow(_spdata.Uvec.Norm() / U0, pow_law_U);
   Kappa[0] = kap_rat * Kappa[1];
};

/*!
\author Juan G Alonso Guzman
\date 09/30/2024
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\return double       Directional derivative
\note This is meant to be called after GetComponent() for the component for which the derivative is wanted
*/
double DiffusionFlowPowerLaw::GetDirectionalDerivative(int xyz)
{
// Note that this doesn't work in regions were the flow is nearly zero.
   return Kappa[comp_eval] * pow_law_U * (_spdata.gradUvec.row[xyz] * _spdata.Uvec) / Sqr(_spdata.Uvec.Norm());
};

/*!
\author Juan G Alonso Guzman
\date 06/03/2024
\return double       Derivative in mu
*/
double DiffusionFlowPowerLaw::GetMuDerivative(void)
{
   return 0.0;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionMomentumPowerLaw methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 08/17/2023
*/
DiffusionMomentumPowerLaw::DiffusionMomentumPowerLaw(void)
                         : DiffusionBase(diff_name_momentum_power_law, 0, DIFF_NOBACKGROUND)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
DiffusionMomentumPowerLaw::DiffusionMomentumPowerLaw(const DiffusionMomentumPowerLaw& other)
                         : DiffusionBase(other)
{
   RAISE_BITS(_status, DIFF_NOBACKGROUND);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Juan G Alonso Guzman
\date 01/04/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void DiffusionMomentumPowerLaw::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(kappa0);
   container.Read(p0);
   container.Read(pow_law_p);
   container.Read(kap_rat);
};

/*!
\author Juan G Alonso Guzman
\date 01/04/2024
*/
void DiffusionMomentumPowerLaw::EvaluateDiffusion(void)
{
   if (comp_eval == 2) return;
   Kappa[1] = kappa0 * pow(_mom[0] / p0, pow_law_p);
   Kappa[0] = kap_rat * Kappa[1];
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\return double       Directional derivative
\note This is meant to be called after GetComponent() for the component for which the derivative is wanted
*/
double DiffusionMomentumPowerLaw::GetDirectionalDerivative(int xyz)
{
   return 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
double DiffusionMomentumPowerLaw::GetMuDerivative(void)
{
   return 0.0;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionKineticEnergyRadialDistancePowerLaw methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 08/18/2023
*/
DiffusionKineticEnergyRadialDistancePowerLaw::DiffusionKineticEnergyRadialDistancePowerLaw(void)
                                            : DiffusionBase(diff_name_kinetic_energy_radial_distance_power_law, 0, DIFF_NOBACKGROUND)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
DiffusionKineticEnergyRadialDistancePowerLaw::DiffusionKineticEnergyRadialDistancePowerLaw(const DiffusionKineticEnergyRadialDistancePowerLaw& other)
                                            : DiffusionBase(other)
{
   RAISE_BITS(_status, DIFF_NOBACKGROUND);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Juan G Alonso Guzman
\date 08/18/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void DiffusionKineticEnergyRadialDistancePowerLaw::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(kap0);
   container.Read(T0);
   container.Read(r0);
   container.Read(pow_law_T);
   container.Read(pow_law_r);
   container.Read(kap_rat);
};

/*!
\author Juan G Alonso Guzman
\date 08/18/2023
*/
void DiffusionKineticEnergyRadialDistancePowerLaw::EvaluateDiffusion(void)
{
   if (comp_eval == 2) return;
   Kappa[1] = kap0 * pow(EnrKin(_mom[0], specie) / T0, pow_law_T) * pow(_pos.Norm() / r0, pow_law_r);
   Kappa[0] = kap_rat * Kappa[1];
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
double DiffusionKineticEnergyRadialDistancePowerLaw::GetMuDerivative(void)
{
   return 0.0;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionRigidityMagneticFieldPowerLaw methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 08/17/2023
*/
DiffusionRigidityMagneticFieldPowerLaw::DiffusionRigidityMagneticFieldPowerLaw(void)
                                      : DiffusionBase(diff_name_rigidity_magnetic_field_power_law, 0, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
DiffusionRigidityMagneticFieldPowerLaw::DiffusionRigidityMagneticFieldPowerLaw(const DiffusionRigidityMagneticFieldPowerLaw& other)
                                      : DiffusionBase(other)
{
   RAISE_BITS(_status, STATE_NONE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Juan G Alonso Guzman
\date 01/04/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void DiffusionRigidityMagneticFieldPowerLaw::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(lam0);
   container.Read(R0);
   container.Read(B0);
   container.Read(pow_law_R);
   container.Read(pow_law_B);
   container.Read(kap_rat);
};

/*!
\author Juan G Alonso Guzman
\date 01/04/2024
*/
void DiffusionRigidityMagneticFieldPowerLaw::EvaluateDiffusion(void)
{
   if (comp_eval == 2) return;
   Kappa[1] = (lam0 * vmag / 3.0) * pow(Rigidity(_mom[0], specie) / R0, pow_law_R) * pow(_spdata.Bmag / B0, pow_law_B);
   Kappa[0] = kap_rat * Kappa[1];
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
double DiffusionRigidityMagneticFieldPowerLaw::GetMuDerivative(void)
{
   return 0.0;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionStraussEtAl2013 methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 12/06/2023
*/
DiffusionStraussEtAl2013::DiffusionStraussEtAl2013(void)
                        : DiffusionBase(diff_name_strauss_et_al_2013, 0, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/06/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
DiffusionStraussEtAl2013::DiffusionStraussEtAl2013(const DiffusionStraussEtAl2013& other)
                        : DiffusionBase(other)
{
   RAISE_BITS(_status, STATE_NONE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Juan G Alonso Guzman
\date 03/12/2024
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void DiffusionStraussEtAl2013::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(LISM_idx);
   container.Read(lam_in);
   container.Read(lam_out);
   container.Read(R0);
   container.Read(B0);
   container.Read(kap_rat_in);
   container.Read(kap_rat_out);
   container.Read(Bmix_idx);
   container.Read(kap_rat_red);
};

/*!
\author Juan G Alonso Guzman
\date 08/01/2024
*/
void DiffusionStraussEtAl2013::EvaluateDiffusion(void)
{
   if (comp_eval == 2) return;

// Find LISM indicator variable (convert -1:1 to 1:0) and interpolate inner/outer quantities. The "Cube" is to bias the indicator variable towards zero (inner heliosphere).
//   LISM_ind = Cube(fmin(fmax(0.0, -0.5 * _spdata.region[LISM_idx] + 0.5), 1.0));
   if (LISM_idx < 0) LISM_ind = 0.0;
   else LISM_ind = (_spdata.region[LISM_idx] > 0.0 ? 0.0 : 1.0);
   double lam_para = LISM_ind * lam_out + (1.0 - LISM_ind) * lam_in;
   double B0_eff = LISM_ind * _spdata.Bmag + (1.0 - LISM_ind) * B0;
   double rig = Rigidity(_mom[0], specie);
   Kappa[1] = (lam_para * vmag / 3.0) * (rig < R0 ? cbrt(rig / R0) : rig / R0) * (B0_eff / _spdata.Bmag);

// Find magnetic mixing indicator variable (convert -1:1 to 0:1) and interpolate perp-to-para diffusion ratio.
//   Bmix_ind = Cube(fmin(fmax(0.0, 0.5 * _spdata.region[Bmix_idx] + 0.5), 1.0));
   if (Bmix_idx < 0) Bmix_ind = 1.0;
   Bmix_ind = (_spdata.region[Bmix_idx] < 0.0 ? 0.0 : 1.0);
   double kap_rat = LISM_ind * kap_rat_out + (1.0 - LISM_ind) * kap_rat_in;
// Reduction factor based on lack of magnetic mixing (i.e. unipolar regions)
   kap_rat *= Bmix_ind + (1.0 - Bmix_ind) * kap_rat_red;
   Kappa[0] = kap_rat * Kappa[1];
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
double DiffusionStraussEtAl2013::GetMuDerivative(void)
{
   return 0.0;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionPotgieterEtAl2015 methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
*/
DiffusionPotgieterEtAl2015::DiffusionPotgieterEtAl2015(void)
                          : DiffusionBase(diff_name_potgieter_et_al_2015, 0, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
DiffusionPotgieterEtAl2015::DiffusionPotgieterEtAl2015(const DiffusionPotgieterEtAl2015& other)
                        : DiffusionBase(other)
{
   RAISE_BITS(_status, STATE_NONE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void DiffusionPotgieterEtAl2015::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(LISM_idx);
   container.Read(lam_in);
   container.Read(lam_out);
   container.Read(R0);
   container.Read(B0);
   container.Read(kap_rat_in);
   container.Read(kap_rat_out);
   container.Read(Bmix_idx);
   container.Read(kap_rat_red);
};

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
*/
void DiffusionPotgieterEtAl2015::EvaluateDiffusion(void)
{
   if (comp_eval == 2) return;

// Find LISM indicator variable (convert -1:1 to 1:0) and interpolate inner/outer quantities. The "Cube" is to bias the indicator variable towards zero (inner heliosphere).
   // LISM_ind = Cube(fmin(fmax(0.0, -0.5 * _spdata.region[LISM_idx] + 0.5), 1.0));
   if (LISM_idx < 0) LISM_ind = 0.0;
   else LISM_ind = (_spdata.region[LISM_idx] > 0.0 ? 0.0 : 1.0);
   double lam_para = LISM_ind * lam_out + (1.0 - LISM_ind) * lam_in;
   double B0_eff = LISM_ind * _spdata.Bmag + (1.0 - LISM_ind) * B0;
   double rig = Rigidity(_mom[0], specie);
   Kappa[1] = (lam_para * vmag / 3.0) * (rig < R0 ? 1.0 : sqrt(Cube(rig / R0))) * (B0_eff / _spdata.Bmag);

// Find magnetic mixing indicator variable (convert -1:1 to 0:1) and interpolate perp-to-para diffusion ratio.
   // Bmix_ind = Cube(fmin(fmax(0.0, 0.5 * _spdata.region[Bmix_idx] + 0.5), 1.0));
   if (Bmix_idx < 0) Bmix_ind = 1.0;
   Bmix_ind = (_spdata.region[Bmix_idx] < 0.0 ? 0.0 : 1.0);
   double kap_rat = LISM_ind * kap_rat_out + (1.0 - LISM_ind) * kap_rat_in;
// Reduction factor based on lack of magnetic mixing (i.e. unipolar regions)
   kap_rat *= Bmix_ind + (1.0 - Bmix_ind) * kap_rat_red;
   Kappa[0] = kap_rat * Kappa[1];
};

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
\return double       Derivative in mu
*/
double DiffusionPotgieterEtAl2015::GetMuDerivative(void)
{
   return 0.0;
};

};
