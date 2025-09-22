/*!
\file diffusion_other.cc
\brief Implements several classes to compute difusion coefficients
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "diffusion_other.hh"
#ifdef USE_GSL
#include <gsl/gsl_sf_hyperg.h>
#endif

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionIsotropicConstant methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/06/2022
*/
template <typename Trajectory>
DiffusionIsotropicConstant<Trajectory>::DiffusionIsotropicConstant(void)
                          : DiffusionBase(diff_name, DIFF_NOBACKGROUND)
{
};

/*!
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename Trajectory>
DiffusionIsotropicConstant<Trajectory>::DiffusionIsotropicConstant(const DiffusionIsotropicConstant& other)
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
template <typename Trajectory>
void DiffusionIsotropicConstant<Trajectory>::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(D0);
};

/*!
\author Vladimir Florinski
\date 05/06/2022
*/
template <typename Trajectory>
void DiffusionIsotropicConstant<Trajectory>::EvaluateDiffusion(int comp)
{
   if ((comp == 0) || (comp == 1)) return;
   Kappa[2] = D0 * (1.0 - Sqr(mu));
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinsk
\date 10/19/2022
\return double       Derivative in mu
*/
template <typename Trajectory>
double DiffusionIsotropicConstant<Trajectory>::GetMuDerivative(int comp)
{
   return -2.0 * D0 * mu;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionQLTConstant methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/09/2022
*/
template <typename Trajectory>
DiffusionQLTConstant<Trajectory>::DiffusionQLTConstant(void)
                    : DiffusionBase(diff_name, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 05/27/2022
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename Trajectory>
DiffusionQLTConstant<Trajectory>::DiffusionQLTConstant(const std::string& name_in, uint16_t status_in)
                    : DiffusionBase(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename Trajectory>
DiffusionQLTConstant<Trajectory>::DiffusionQLTConstant(const DiffusionQLTConstant& other)
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
template <typename Trajectory>
void DiffusionQLTConstant<Trajectory>::SetupDiffusion(bool construct)
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
template <typename Trajectory>
void DiffusionQLTConstant<Trajectory>::EvaluateDiffusion(int comp)
{
   if ((comp == 0) || (comp == 1)) return;
   Kappa[2] = 0.25 * M_PI * ps_minus * fabs(Omega) * st2 * pow(vmag * k_min * fabs(mu / Omega), ps_minus) * A2A;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionWNLTConstant methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/06/2022
*/
template <typename Trajectory>
DiffusionWNLTConstant<Trajectory>::DiffusionWNLTConstant(void)
                     : DiffusionQLTConstant(diff_name, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 05/27/2022
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename Trajectory>
DiffusionWNLTConstant<Trajectory>::DiffusionWNLTConstant(const std::string& name_in, uint16_t status_in)
                     : DiffusionQLTConstant(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename Trajectory>
DiffusionWNLTConstant<Trajectory>::DiffusionWNLTConstant(const DiffusionWNLTConstant& other)
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
template <typename Trajectory>
void DiffusionWNLTConstant<Trajectory>::SetupDiffusion(bool construct)
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
template <typename Trajectory>
void DiffusionWNLTConstant<Trajectory>::EvaluateDiffusion(int comp)
{
   if (comp == 1) return;
   double CT, CL, xi1, xi2, F21, DT1, DT2;

// Kappa[0] is required for Kappa[2]
   CT = 0.5 * ps_minus / ps_plus * Sqr(mu / k_min) * A2T;
   CL = 0.125 * Sqr(vmag * st2 / Omega) * A2L;
   Kappa[0] = vmag * sqrt(CT + CL);

   if (comp == 0) return;

// Hypergeometric function may crash if the last argument is close to 1
   DiffusionQLTConstant::EvaluateDiffusion(comp);

   if constexpr (HConfig::use_qlt_scatt) {
      return;
   }
   else {
      if (A2T > sp_tiny) {
         xi1 = vmag * k_min * sqrt(st2) / M_SQRT2 / fabs(Omega);
         F21 = gsl_sf_hyperg_2F1(1.0, 1.0, (5.0 + ps_index) / 4.0, 1.0 / (1.0 + Sqr(xi1)));
         DT1 = 1.0 / (1.0 + Sqr(xi1)) * F21;
         xi2 = Kappa[0] * k_min * k_min / fabs(Omega);
         F21 = gsl_sf_hyperg_2F1(1.0, 1.0, (5.0 + ps_index) / 4.0, 1.0 / (1.0 + Sqr(xi2)));
         DT2 = 1.0 / (1.0 + Sqr(xi2)) * F21;
         Kappa[2] += (ps_minus / ps_plus) * Omega * st2 * A2T * DT1 * DT2 / sqrt(Sqr(DT1) + Sqr(DT2));
      };
   }
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionWNLTRampVLISM methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/06/2022
*/
template <typename Trajectory>
DiffusionWNLTRampVLISM<Trajectory>::DiffusionWNLTRampVLISM(void)
                      : DiffusionWNLTConstant(diff_name, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename Trajectory>
DiffusionWNLTRampVLISM<Trajectory>::DiffusionWNLTRampVLISM(const DiffusionWNLTRampVLISM& other)
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
template <typename Trajectory>
void DiffusionWNLTRampVLISM<Trajectory>::SetupDiffusion(bool construct)
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
template <typename Trajectory>
void DiffusionWNLTRampVLISM<Trajectory>::EvaluateDiffusion(int comp)
{
   if (comp == 1) return;
   double z0, r0, k_ratio;

// Find distance to nose for Rankine half body which the particle is presently on
   r0 = _coords.Pos().Norm();
   z0 = sqrt(0.5 * r0 * (r0 +_coords.Pos()[2]));

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
   DiffusionWNLTConstant::EvaluateDiffusion(comp);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionParaConstant methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/07/2023
*/
template <typename Trajectory>
DiffusionParaConstant<Trajectory>::DiffusionParaConstant(void)
                     : DiffusionBase(diff_name, DIFF_NOBACKGROUND)
{
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/07/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename Trajectory>
DiffusionParaConstant<Trajectory>::DiffusionParaConstant(const DiffusionParaConstant& other)
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
template <typename Trajectory>
void DiffusionParaConstant<Trajectory>::SetupDiffusion(bool construct)
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
template <typename Trajectory>
void DiffusionParaConstant<Trajectory>::EvaluateDiffusion(int comp)
{
   if ((comp == 0) || (comp == 2)) return;
   Kappa[1] = D0;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/07/2023
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\param[in] ddata_in Derivative data from computing background fields
\return double       Directional derivative
 \note This must be called after Stage() if the target coordinates have changed.
*/
template <typename Trajectory>
double DiffusionParaConstant<Trajectory>::GetDirectionalDerivative(int comp, int xyz, const DerivativeData& ddata)
{
   return 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
template <typename Trajectory>
double DiffusionParaConstant<Trajectory>::GetMuDerivative(int comp)
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
template <typename Trajectory>
DiffusionPerpConstant<Trajectory>::DiffusionPerpConstant(void)
                     : DiffusionBase(diff_name, DIFF_NOBACKGROUND)
{
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename Trajectory>
DiffusionPerpConstant<Trajectory>::DiffusionPerpConstant(const DiffusionPerpConstant& other)
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
template <typename Trajectory>
void DiffusionPerpConstant<Trajectory>::SetupDiffusion(bool construct)
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
template <typename Trajectory>
void DiffusionPerpConstant<Trajectory>::EvaluateDiffusion(int comp)
{
   if ((comp == 1) || (comp == 2)) return;
   Kappa[0] = D0;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/18/2022
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\param[in] ddata_in Derivative data from computing background fields
\return double       Directional derivative
 \note This must be called after Stage() if the target coordinates have changed.
*/
template <typename Trajectory>
double DiffusionPerpConstant<Trajectory>::GetDirectionalDerivative(int comp, int xyz, const DerivativeData& ddata)
{
   return 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
template <typename Trajectory>
double DiffusionPerpConstant<Trajectory>::GetMuDerivative(int comp)
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
template <typename Trajectory>
DiffusionFullConstant<Trajectory>::DiffusionFullConstant(void)
                     : DiffusionBase(diff_name, DIFF_NOBACKGROUND)
{
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename Trajectory>
DiffusionFullConstant<Trajectory>::DiffusionFullConstant(const DiffusionFullConstant& other)
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
template <typename Trajectory>
void DiffusionFullConstant<Trajectory>::SetupDiffusion(bool construct)
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
template <typename Trajectory>
void DiffusionFullConstant<Trajectory>::EvaluateDiffusion(int comp)
{
   if (comp == 2) return;
   Kappa[0] = Dperp;
   Kappa[1] = Dpara;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 10/18/2022
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\param[in] ddata_in Derivative data from computing background fields
\return double       Directional derivative
 \note This must be called after Stage() if the target coordinates have changed.
*/
template <typename Trajectory>
double DiffusionFullConstant<Trajectory>::GetDirectionalDerivative(int comp, int xyz, const DerivativeData& ddata)
{
   return 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
template <typename Trajectory>
double DiffusionFullConstant<Trajectory>::GetMuDerivative(int comp)
{
   return 0.0;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionFlowMomentumPowerLaw methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Swati Sharma
\date 01/03/2025
*/
template <typename Trajectory>
DiffusionFlowMomentumPowerLaw<Trajectory>::DiffusionFlowMomentumPowerLaw(void)
                             : DiffusionBase(diff_name, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\author Swati Sharma
\date 01/03/2025
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename Trajectory>
DiffusionFlowMomentumPowerLaw<Trajectory>::DiffusionFlowMomentumPowerLaw(const DiffusionFlowMomentumPowerLaw& other)
                             : DiffusionBase(other)
{
   RAISE_BITS(_status, STATE_NONE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Juan G Alonso Guzman
\author Swati Sharma
\date 01/03/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename Trajectory>
void DiffusionFlowMomentumPowerLaw<Trajectory>::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(kap0);
   container.Read(U0);
   container.Read(pow_law_U);
   container.Read(p0);
   container.Read(pow_law_p);
   container.Read(kap_rat);
};

/*!
\author Juan G Alonso Guzman
\author Swati Sharma
\date 01/03/2025
*/
template <typename Trajectory>
void DiffusionFlowMomentumPowerLaw<Trajectory>::EvaluateDiffusion(int comp)
{
   if (comp == 2) return;
   Kappa[1] = kap0 * pow(_fields.Vel().Norm() / U0, pow_law_U) * pow(_mom[0] / p0, pow_law_p);
   Kappa[0] = kap_rat * Kappa[1];
};

/*!
\author Juan G Alonso Guzman
\author Swati Sharma
\date 01/03/2025
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\param[in] ddata_in Derivative data from computing background fields
\return double       Directional derivative
 \note This must be called after Stage() if the target coordinates have changed.
*/
template <typename Trajectory>
double DiffusionFlowMomentumPowerLaw<Trajectory>::GetDirectionalDerivative(int comp, int xyz, const DerivativeData& ddata)
{
// Note that this doesn't work in regions were the flow is nearly zero.
   if ((0 <= xyz) && (xyz <= 2)) return Kappa[comp] * pow_law_U * (_fields.DelVel().row[xyz] * _fields.Vel()) / Sqr(_fields.Vel().Norm());
   else return Kappa[comp] * pow_law_U * (_fields.DotVel() * _fields.Vel()) / Sqr(_fields.Vel().Norm());
};

/*!
\author Juan G Alonso Guzman
\author Swati Sharma
\date 01/03/2025
\return double       Derivative in mu
*/
template <typename Trajectory>
double DiffusionFlowMomentumPowerLaw<Trajectory>::GetMuDerivative(int comp)
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
template <typename Trajectory>
DiffusionKineticEnergyRadialDistancePowerLaw<Trajectory>::DiffusionKineticEnergyRadialDistancePowerLaw(void)
                                            : DiffusionBase(diff_name, DIFF_NOBACKGROUND)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename Trajectory>
DiffusionKineticEnergyRadialDistancePowerLaw<Trajectory>::DiffusionKineticEnergyRadialDistancePowerLaw(const DiffusionKineticEnergyRadialDistancePowerLaw& other)
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
template <typename Trajectory>
void DiffusionKineticEnergyRadialDistancePowerLaw<Trajectory>::SetupDiffusion(bool construct)
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
template <typename Trajectory>
void DiffusionKineticEnergyRadialDistancePowerLaw<Trajectory>::EvaluateDiffusion(int comp)
{
   if (comp == 2) return;
   Kappa[1] = kap0 * pow(EnrKin(_coords.Mom()[0], specie) / T0, pow_law_T) * pow(_coords.Pos().Norm() / r0, pow_law_r);
   Kappa[0] = kap_rat * Kappa[1];
};

/*!
\author Juan G Alonso Guzman
\date 02/18/2025
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\param[in] ddata_in Derivative data from computing background fields
\return double       Directional derivative
 \note This must be called after Stage() if the target coordinates have changed.
*/
template <typename Trajectory>
double DiffusionKineticEnergyRadialDistancePowerLaw<Trajectory>::GetDirectionalDerivative(int comp, int xyz, const DerivativeData& ddata)
{
// Note that this doesn't work near the origin where the radial distance is close to zero.
   if ((0 <= xyz) && (xyz <= 2)) return Kappa[comp] * pow_law_r * _coords.Pos()[xyz] / Sqr(_coords.Pos().Norm());
   else return 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
template <typename Trajectory>
double DiffusionKineticEnergyRadialDistancePowerLaw<Trajectory>::GetMuDerivative(int comp)
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
template <typename Trajectory>
DiffusionRigidityMagneticFieldPowerLaw<Trajectory>::DiffusionRigidityMagneticFieldPowerLaw(void)
                                      : DiffusionBase(diff_name, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename Trajectory>
DiffusionRigidityMagneticFieldPowerLaw<Trajectory>::DiffusionRigidityMagneticFieldPowerLaw(const DiffusionRigidityMagneticFieldPowerLaw& other)
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
template <typename Trajectory>
void DiffusionRigidityMagneticFieldPowerLaw<Trajectory>::SetupDiffusion(bool construct)
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
template <typename Trajectory>
void DiffusionRigidityMagneticFieldPowerLaw<Trajectory>::EvaluateDiffusion(int comp)
{
   if (comp == 2) return;
   Kappa[1] = (lam0 * vmag / 3.0) * pow(Rigidity(_coords.Mom()[0], specie) / R0, pow_law_R) * pow(_fields.AbsMag() / B0, pow_law_B);
   Kappa[0] = kap_rat * Kappa[1];
};

/*!
\author Juan G Alonso Guzman
\date 02/18/2025
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\param[in] ddata_in Derivative data from computing background fields
\return double       Directional derivative
 \note This must be called after Stage() if the target coordinates have changed.
*/
template <typename Trajectory>
double DiffusionRigidityMagneticFieldPowerLaw<Trajectory>::GetDirectionalDerivative(int comp, int xyz, const DerivativeData& ddata)
{
// Note that this doesn't work in regions were the field is nearly zero, although in that case an error would be thrown elsewhere in the code.
   if ((0 <= xyz) && (xyz <= 2)) return Kappa[comp] * pow_law_B * _fields.DelAbsMag()[xyz] / _fields.AbsMag();
   else return Kappa[comp] * pow_law_B * _fields.DotAbsMag() / _fields.AbsMag();
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
template <typename Trajectory>
double DiffusionRigidityMagneticFieldPowerLaw<Trajectory>::GetMuDerivative(int comp)
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
template <typename Trajectory>
DiffusionStraussEtAl2013<Trajectory>::DiffusionStraussEtAl2013(void)
                        : DiffusionBase(diff_name, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/06/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename Trajectory>
DiffusionStraussEtAl2013<Trajectory>::DiffusionStraussEtAl2013(const DiffusionStraussEtAl2013& other)
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
template <typename Trajectory>
void DiffusionStraussEtAl2013<Trajectory>::SetupDiffusion(bool construct)
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
};

/*!
\author Juan G Alonso Guzman
\date 08/01/2024
*/
template <typename Trajectory>
void DiffusionStraussEtAl2013<Trajectory>::EvaluateDiffusion(int comp)
{
   if (comp == 2) return;

// Find LISM indicator variable (convert -1:1 to 1:0) and interpolate inner/outer quantities.
   if (LISM_idx < 0) LISM_ind = 0.0;
   else LISM_ind = (_fields.IvLISM() > 0.0 ? 0.0 : 1.0);
   double lam_para = LISM_ind * lam_out + (1.0 - LISM_ind) * lam_in;
   double B0_eff = LISM_ind * _fields.Mag() + (1.0 - LISM_ind) * B0;
   double rig = Rigidity(_coords.Mom()[0], specie);
   double kap_rat;

// Find diffusion coefficients
   Kappa[1] = (lam_para * vmag / 3.0) * (rig < R0 ? cbrt(rig / R0) : rig / R0) * (B0_eff / _fields.Mag());
   if (comp == 0) {
      kap_rat = LISM_ind * kap_rat_out + (1.0 - LISM_ind) * kap_rat_in;
      Kappa[0] = kap_rat * Kappa[1];
   };
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
template <typename Trajectory>
double DiffusionStraussEtAl2013<Trajectory>::GetMuDerivative(int comp)
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
template <typename Trajectory>
DiffusionPotgieterEtAl2015<Trajectory>::DiffusionPotgieterEtAl2015(void)
                          : DiffusionBase(diff_name, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename Trajectory>
DiffusionPotgieterEtAl2015<Trajectory>::DiffusionPotgieterEtAl2015(const DiffusionPotgieterEtAl2015& other)
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
template <typename Trajectory>
void DiffusionPotgieterEtAl2015<Trajectory>::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(LISM_idx);
   container.Read(kappa_in);
   container.Read(kappa_out);
   container.Read(R0);
   container.Read(B0);
   container.Read(kap_rat_in);
   container.Read(kap_rat_out);
};

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
*/
template <typename Trajectory>
void DiffusionPotgieterEtAl2015<Trajectory>::EvaluateDiffusion(int comp)
{
   if (comp == 2) return;

// Find LISM indicator variable (convert -1:1 to 1:0) and interpolate inner/outer quantities.
   if (LISM_idx < 0) LISM_ind = 0.0;
   else LISM_ind = (_fields.IvLISM() > 0.0 ? 0.0 : 1.0);
   double kappa_para = LISM_ind * kappa_out + (1.0 - LISM_ind) * kappa_in;
   double B0_eff = LISM_ind * _fields.Mag() + (1.0 - LISM_ind) * B0;
   double rig = Rigidity(_coords.Mom()[0], specie);
   double kap_rat;

// Find diffusion coefficients
   Kappa[1] = kappa_para * (vmag / c_code) * (rig < R0 ? 1.0 : sqrt(Cube(rig / R0))) * (B0_eff / _fields.Mag());
   if (comp == 0) {
      kap_rat = LISM_ind * kap_rat_out + (1.0 - LISM_ind) * kap_rat_in;
      Kappa[0] = kap_rat * Kappa[1];
   };
};

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
\return double       Derivative in mu
*/
template <typename Trajectory>
double DiffusionPotgieterEtAl2015<Trajectory>::GetMuDerivative(int comp)
{
   return 0.0;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionEmpiricalSOQLTandUNLT methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
*/
template <typename Trajectory>
DiffusionEmpiricalSOQLTandUNLT<Trajectory>::DiffusionEmpiricalSOQLTandUNLT(void)
                              : DiffusionBase(diff_name, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
\param[in] other Object to initialize from
A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename Trajectory>
DiffusionEmpiricalSOQLTandUNLT<Trajectory>::DiffusionEmpiricalSOQLTandUNLT(const DiffusionEmpiricalSOQLTandUNLT& other)
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
template <typename Trajectory>
void DiffusionEmpiricalSOQLTandUNLT<Trajectory>::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(lam_para);
   container.Read(lam_perp);
   container.Read(R0);
   container.Read(B0);
   container.Read(Bmix_idx);
   container.Read(kap_rat_red);
   container.Read(radial_limit_perp_red);
   container.Read(solar_cycle_idx);
   container.Read(solar_cycle_effect);
};

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
*/
template <typename Trajectory>
void DiffusionEmpiricalSOQLTandUNLT<Trajectory>::EvaluateDiffusion(int comp)
{
   if (comp == 2) return;
   double lam, rig_dep;
   double rig = Rigidity(_coords.Mom()[0], specie);

   if (comp == 1) {
// Compute mean free path and rigidity dependance with a bent power law
      rig_dep = cbrt((rig / R0) * (1.0 + Sqr(Sqr(rig / R0))));
      lam = lam_para;
   }
   else if (comp == 0) {
// Compute mean free path and rigidity dependance with a bent power law
      rig_dep = cbrt((rig / R0) * (1.0 + Sqr(rig / R0)));

// Find magnetic mixing indicator variable (convert -1:1 to 0:1) and interpolate perp-to-para diffusion ratio.
      if (Bmix_idx < 0) Bmix_ind = 1.0;
      Bmix_ind = (_fields.IvBmix() < 0.0 ? 0.0 : 1.0);
      if (_coords.Pos().Norm() < radial_limit_perp_red) lam = lam_perp * (Bmix_ind + (1.0 - Bmix_ind) * kap_rat_red);
      else lam = lam_perp;
   };
   Kappa[comp] = (lam * vmag / 3.0) * rig_dep * (B0 / _fields.Mag());
   Kappa[comp] /= 1.0 + solar_cycle_effect * Sqr(cos(0.5 * _fields.IvSolarCycle()));
};

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
\return double       Derivative in mu
*/
template <typename Trajectory>
double DiffusionEmpiricalSOQLTandUNLT<Trajectory>::GetMuDerivative(int comp)
{
   return 0.0;
};

};
