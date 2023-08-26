/*!
\file diffusion_other.cc
\brief Implements several classes to compute difusion coefficients
\author Vladimir Florinski

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
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
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
   if(!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(&D0);
};

/*!
\author Vladimir Florinski
\date 05/06/2022
*/
void DiffusionIsotropicConstant::EvaluateDiffusion(void)
{
   if((comp_eval == 0) || (comp_eval == 1)) return;
   Kappa[2] = D0 * (1.0 - Sqr(mu));
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinsk
\date 10/19/2022
\return double       Derivative in mu
\note This is meant to be called after GetComponent() for the componenent for which the derivative is wanted
*/
double DiffusionIsotropicConstant::GetMuDerivative(void)
{
   return -2.0 * D0 * mu;
}

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
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
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
   if(!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(&A2A);
   container.Read(&l_max);
   container.Read(&ps_index);
   k_min = twopi / l_max;
   ps_minus = ps_index - 1.0;
};

/*!
\author Vladimir Florinski
\date 05/09/2022
*/
void DiffusionQLTConstant::EvaluateDiffusion(void)
{
   if((comp_eval == 0) || (comp_eval == 1)) return;
   Kappa[2] = 0.25 * M_PI * ps_minus * fabs(Omega) * st2 * pow(vmag * k_min * fabs(mu / Omega), ps_minus) * A2A;
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
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
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
   if(!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(&D0);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/07/2023
*/
void DiffusionParaConstant::EvaluateDiffusion(void)
{
   if((comp_eval == 0) || (comp_eval == 2)) return;
   Kappa[1] = D0;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/07/2023
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\return double       Directional derivative
\note This is meant to be called after GetComponent() for the componenent for which the derivative is wanted
*/
double DiffusionParaConstant::GetDirectionalDerivative(int xyz)
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
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
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
   if(!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(&D0);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/06/2022
*/
void DiffusionPerpConstant::EvaluateDiffusion(void)
{
   if((comp_eval == 1) || (comp_eval == 2)) return;
   Kappa[0] = D0;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/18/2022
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\return double       Directional derivative
\note This is meant to be called after GetComponent() for the componenent for which the derivative is wanted
*/
double DiffusionPerpConstant::GetDirectionalDerivative(int xyz)
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
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
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
   if(!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(&Dperp);
   container.Read(&Dpara);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/06/2022
*/
void DiffusionFullConstant::EvaluateDiffusion(void)
{
   if((comp_eval == 2)) return;
   Kappa[0] = Dperp;
   Kappa[1] = Dpara;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/18/2022
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\return double       Directional derivative
\note This is meant to be called after GetComponent() for the componenent for which the derivative is wanted
*/
double DiffusionFullConstant::GetDirectionalDerivative(int xyz)
{
   return 0.0;
};

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
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
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
   if(!construct) DiffusionQLTConstant::SetupDiffusion(false);
   container.Read(&A2T);
   container.Read(&A2L);
   ps_plus = ps_index + 1.0;
};

/*!
\author Vladimir Florinski
\date 05/09/2022
*/
void DiffusionWNLTConstant::EvaluateDiffusion(void)
{
   if(comp_eval == 1) return;
   double CT, CL, xi1, xi2, F21, DT1, DT2;

// Kappa[0] is required for Kappa[2]
   CT = 0.5 * ps_minus / ps_plus * Sqr(mu / k_min) * A2T;
   CL = 0.125 * Sqr(vmag * st2 / Omega) * A2L;
   Kappa[0] = vmag * sqrt(CT + CL);

   if(comp_eval == 0) return;

// Hypergeometric function may crash if the last argument is close to 1
   DiffusionQLTConstant::EvaluateDiffusion();

#ifdef USE_QLT_SCATT_WITH_WNLT_DIFF
   return;
#endif

   if(A2T > tiny) {
      xi1 = vmag * k_min * sqrt(st2) / sqrttwo / fabs(Omega);
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
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
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
   if(!construct) DiffusionWNLTConstant::SetupDiffusion(false);
   container.Read(&l_max_HP);
   container.Read(&z_nose);
   container.Read(&z_sheath);
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
   if(comp_eval == 1) return;
   double z0, r0, k_ratio;

// Find distance to nose for Rankine half body which the particle is presently on
   r0 = _pos.Norm();
   z0 = sqrt(0.5 * r0 * (r0 + _pos[2]));
// Constant k_min beyond z_sheath
   if(z0 > z_sheath) {
      k_min = k_min_ref;
      A2A = A2A_ref;
      A2T = A2T_ref;
      A2L = A2L_ref;
   }
// Linearly interpolate l_max between z_nose and z_sheath
   else {
      k_min = twopi / (l_max_HP + dl_max * (z0 - z_nose) / dz);
      k_ratio = pow(k_min_ref / k_min, ps_minus);
      A2A = A2A_ref * k_ratio;
      A2T = A2T_ref * k_ratio;
      A2L = A2L_ref * k_ratio;
   };

// Evaluate WLNT diffusion
   DiffusionWNLTConstant::EvaluateDiffusion();
};

#endif

};
