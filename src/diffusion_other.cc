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
template <typename HConfig>
DiffusionIsotropicConstant<HConfig>::DiffusionIsotropicConstant(void)
                          : DiffusionBase(diff_name, DIFF_NOBACKGROUND)
{
};

/*!
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename HConfig>
DiffusionIsotropicConstant<HConfig>::DiffusionIsotropicConstant(const DiffusionIsotropicConstant& other)
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
template <typename HConfig>
void DiffusionIsotropicConstant<HConfig>::SetupDiffusion(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DiffusionBase::SetupDiffusion(false);
   container.Read(D0);
};

/*!
\author Vladimir Florinski
\date 05/06/2022
*/
template <typename HConfig>
void DiffusionIsotropicConstant<HConfig>::EvaluateDiffusion(Component comp)
{
   if ((comp == Component::perp) || (comp == Component::para)) return;
   Kappa[Component::mu] = D0 * (1.0 - Sqr(_coords.MomMu()));
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinsk
\date 10/19/2022
\return double       Derivative in mu
*/
template <typename HConfig>
double DiffusionIsotropicConstant<HConfig>::GetMuDerivative(Component comp)
{
   return -2.0 * D0 * _coords.MomMu();
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionQLTConstant methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/09/2022
*/
template <typename HConfig>
DiffusionQLTConstant<HConfig>::DiffusionQLTConstant(void)
                    : DiffusionBase(diff_name, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 05/27/2022
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
DiffusionQLTConstant<HConfig>::DiffusionQLTConstant(const std::string& name_in, uint16_t status_in)
                    : DiffusionBase(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename HConfig>
DiffusionQLTConstant<HConfig>::DiffusionQLTConstant(const DiffusionQLTConstant& other)
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
template <typename HConfig>
void DiffusionQLTConstant<HConfig>::SetupDiffusion(bool construct)
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
template <typename HConfig>
void DiffusionQLTConstant<HConfig>::EvaluateDiffusion(Component comp)
{
   if ((comp == Component::perp) || (comp == Component::para)) return;
   Kappa[Component::mu] = 0.25 * M_PI * ps_minus * fabs(Omega) * st2 * pow(_coords.AbsVel() * k_min * fabs(_coords.MomMu() / Omega), ps_minus) * A2A;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionWNLTConstant methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/06/2022
*/
template <typename HConfig>
DiffusionWNLTConstant<HConfig>::DiffusionWNLTConstant(void)
                     : DiffusionQLTConstant(diff_name, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 05/27/2022
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
DiffusionWNLTConstant<HConfig>::DiffusionWNLTConstant(const std::string& name_in, uint16_t status_in)
                     : DiffusionQLTConstant(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename HConfig>
DiffusionWNLTConstant<HConfig>::DiffusionWNLTConstant(const DiffusionWNLTConstant& other)
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
template <typename HConfig>
void DiffusionWNLTConstant<HConfig>::SetupDiffusion(bool construct)
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
template <typename HConfig>
void DiffusionWNLTConstant<HConfig>::EvaluateDiffusion(Component comp)
{
   if (comp == Component::para) return;
   double CT, CL, xi1, xi2, F21, DT1, DT2;

// Kappa[Component::perp] is required for Kappa[Component::mu]
   CT = 0.5 * ps_minus / ps_plus * Sqr(_coords.MomMu() / k_min) * A2T;
   CL = 0.125 * Sqr(_coords.AbsVel() * st2 / Omega) * A2L;
   Kappa[Component::perp] = _coords.AbsVel() * sqrt(CT + CL);

   if (comp == Component::perp) return;

// Hypergeometric function may crash if the last argument is close to 1
   DiffusionQLTConstant::EvaluateDiffusion(comp);

   if constexpr (HConfig::use_qlt_scatt) {
      return;
   }
   else {
      if (A2T > sp_tiny) {
         xi1 = _coords.AbsVel() * k_min * sqrt(st2) / M_SQRT2 / fabs(Omega);
         F21 = gsl_sf_hyperg_2F1(1.0, 1.0, (5.0 + ps_index) / 4.0, 1.0 / (1.0 + Sqr(xi1)));
         DT1 = 1.0 / (1.0 + Sqr(xi1)) * F21;
         xi2 = Kappa[Component::perp] * k_min * k_min / fabs(Omega);
         F21 = gsl_sf_hyperg_2F1(1.0, 1.0, (5.0 + ps_index) / 4.0, 1.0 / (1.0 + Sqr(xi2)));
         DT2 = 1.0 / (1.0 + Sqr(xi2)) * F21;
         Kappa[Component::mu] += (ps_minus / ps_plus) * Omega * st2 * A2T * DT1 * DT2 / sqrt(Sqr(DT1) + Sqr(DT2));
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
template <typename HConfig>
DiffusionWNLTRampVLISM<HConfig>::DiffusionWNLTRampVLISM(void)
                      : DiffusionWNLTConstant(diff_name, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename HConfig>
DiffusionWNLTRampVLISM<HConfig>::DiffusionWNLTRampVLISM(const DiffusionWNLTRampVLISM& other)
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
template <typename HConfig>
void DiffusionWNLTRampVLISM<HConfig>::SetupDiffusion(bool construct)
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
template <typename HConfig>
void DiffusionWNLTRampVLISM<HConfig>::EvaluateDiffusion(Component comp)
{
   if (comp == Component::para) return;
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
template <typename HConfig>
DiffusionParaConstant<HConfig>::DiffusionParaConstant(void)
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
template <typename HConfig>
DiffusionParaConstant<HConfig>::DiffusionParaConstant(const DiffusionParaConstant& other)
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
template <typename HConfig>
void DiffusionParaConstant<HConfig>::SetupDiffusion(bool construct)
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
template <typename HConfig>
void DiffusionParaConstant<HConfig>::EvaluateDiffusion(Component comp)
{
   if ((comp == Component::perp) || (comp == Component::mu)) return;
   Kappa[Component::para] = D0;
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
template <typename HConfig>
double DiffusionParaConstant<HConfig>::GetDirectionalDerivative(Component comp, int xyz, const DerivativeData& ddata)
{
   return 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
template <typename HConfig>
double DiffusionParaConstant<HConfig>::GetMuDerivative(Component comp)
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
template <typename HConfig>
DiffusionPerpConstant<HConfig>::DiffusionPerpConstant(void)
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
template <typename HConfig>
DiffusionPerpConstant<HConfig>::DiffusionPerpConstant(const DiffusionPerpConstant& other)
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
template <typename HConfig>
void DiffusionPerpConstant<HConfig>::SetupDiffusion(bool construct)
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
template <typename HConfig>
void DiffusionPerpConstant<HConfig>::EvaluateDiffusion(Component comp)
{
   if ((comp == Component::para) || (comp == Component::mu)) return;
   Kappa[Component::perp] = D0;
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
template <typename HConfig>
double DiffusionPerpConstant<HConfig>::GetDirectionalDerivative(Component comp, int xyz, const DerivativeData& ddata)
{
   return 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
template <typename HConfig>
double DiffusionPerpConstant<HConfig>::GetMuDerivative(Component comp)
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
template <typename HConfig>
DiffusionFullConstant<HConfig>::DiffusionFullConstant(void)
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
template <typename HConfig>
DiffusionFullConstant<HConfig>::DiffusionFullConstant(const DiffusionFullConstant& other)
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
template <typename HConfig>
void DiffusionFullConstant<HConfig>::SetupDiffusion(bool construct)
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
template <typename HConfig>
void DiffusionFullConstant<HConfig>::EvaluateDiffusion(Component comp)
{
   if (comp == Component::mu) return;
   Kappa[Component::perp] = Dperp;
   Kappa[Component::para] = Dpara;
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
template <typename HConfig>
double DiffusionFullConstant<HConfig>::GetDirectionalDerivative(Component comp, int xyz, const DerivativeData& ddata)
{
   return 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
template <typename HConfig>
double DiffusionFullConstant<HConfig>::GetMuDerivative(Component comp)
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
template <typename HConfig>
DiffusionFlowMomentumPowerLaw<HConfig>::DiffusionFlowMomentumPowerLaw(void)
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
template <typename HConfig>
DiffusionFlowMomentumPowerLaw<HConfig>::DiffusionFlowMomentumPowerLaw(const DiffusionFlowMomentumPowerLaw& other)
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
template <typename HConfig>
void DiffusionFlowMomentumPowerLaw<HConfig>::SetupDiffusion(bool construct)
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
template <typename HConfig>
void DiffusionFlowMomentumPowerLaw<HConfig>::EvaluateDiffusion(Component comp)
{
   if (comp == Component::mu) return;
   // todo Fluv, AbsFluv, DotFluv - only for FlowMomentumPowerLaw!
   Kappa[Component::para] = kap0 * pow(_fields.AbsFluv() / U0, pow_law_U) * pow(_coords.AbsMom() / p0, pow_law_p);
   Kappa[Component::perp] = kap_rat * Kappa[Component::para];
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
template <typename HConfig>
double DiffusionFlowMomentumPowerLaw<HConfig>::GetDirectionalDerivative(Component comp, int xyz, const DerivativeData& ddata)
{
// Note that this doesn't work in regions were the flow is nearly zero.
   if ((0 <= xyz) && (xyz <= 2))
      return Kappa[comp] * pow_law_U * (_fields.DelFluv().row[xyz] * _fields.Fluv()) / Sqr(_fields.AbsFluv());
   else
      return Kappa[comp] * pow_law_U * (_fields.DotFluv() * _fields.Fluv()) / Sqr(_fields.AbsFluv());
};

/*!
\author Juan G Alonso Guzman
\author Swati Sharma
\date 01/03/2025
\return double       Derivative in mu
*/
template <typename HConfig>
double DiffusionFlowMomentumPowerLaw<HConfig>::GetMuDerivative(Component comp)
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
template <typename HConfig>
DiffusionKineticEnergyRadialDistancePowerLaw<HConfig>::DiffusionKineticEnergyRadialDistancePowerLaw(void)
                                            : DiffusionBase(diff_name, DIFF_NOBACKGROUND)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename HConfig>
DiffusionKineticEnergyRadialDistancePowerLaw<HConfig>::DiffusionKineticEnergyRadialDistancePowerLaw(const DiffusionKineticEnergyRadialDistancePowerLaw& other)
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
template <typename HConfig>
void DiffusionKineticEnergyRadialDistancePowerLaw<HConfig>::SetupDiffusion(bool construct)
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
template <typename HConfig>
void DiffusionKineticEnergyRadialDistancePowerLaw<HConfig>::EvaluateDiffusion(Component comp)
{
   if (comp == Component::mu) return;
   Kappa[Component::para] = kap0 * pow(EnrKin<specie>(_coords.AbsMom()) / T0, pow_law_T) * pow(_coords.Rad() / r0, pow_law_r);
   Kappa[Component::perp] = kap_rat * Kappa[Component::para];
};

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 09/28/2025
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\param[in] ddata_in Derivative data from computing background fields
\return double       Directional derivative
 \note This must be called after Stage() if the target coordinates have changed.
*/
template <typename HConfig>
double DiffusionKineticEnergyRadialDistancePowerLaw<HConfig>::GetDirectionalDerivative(Component comp, int xyz, const DerivativeData& ddata)
{
// Note that this doesn't work near the origin where the radial distance is close to zero.
   if ((0 <= xyz) && (xyz <= 2)) return Kappa[comp] * pow_law_r * _coords.Pos()[xyz] / Sqr(_coords.Rad());
   else return 0.0;
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
template <typename HConfig>
double DiffusionKineticEnergyRadialDistancePowerLaw<HConfig>::GetMuDerivative(Component comp)
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
template <typename HConfig>
DiffusionRigidityMagneticFieldPowerLaw<HConfig>::DiffusionRigidityMagneticFieldPowerLaw(void)
                                      : DiffusionBase(diff_name, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/09/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename HConfig>
DiffusionRigidityMagneticFieldPowerLaw<HConfig>::DiffusionRigidityMagneticFieldPowerLaw(const DiffusionRigidityMagneticFieldPowerLaw& other)
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
template <typename HConfig>
void DiffusionRigidityMagneticFieldPowerLaw<HConfig>::SetupDiffusion(bool construct)
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
template <typename HConfig>
void DiffusionRigidityMagneticFieldPowerLaw<HConfig>::EvaluateDiffusion(Component comp)
{
   if (comp == Component::mu) return;
   Kappa[Component::para] = (lam0 * _coords.AbsVel() / 3.0) * pow(Rigidity<specie>(_coords.AbsMom()) / R0, pow_law_R) * pow(_fields.AbsMag() / B0, pow_law_B);
   Kappa[Component::perp] = kap_rat * Kappa[Component::para];
};

/*!
\author Juan G Alonso Guzman
\date 02/18/2025
\param[in] xyz       Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\param[in] ddata_in Derivative data from computing background fields
\return double       Directional derivative
 \note This must be called after Stage() if the target coordinates have changed.
*/
template <typename HConfig>
double DiffusionRigidityMagneticFieldPowerLaw<HConfig>::GetDirectionalDerivative(Component comp, int xyz, const DerivativeData& ddata)
{
// Note that this doesn't work in regions were the field is nearly zero, although in that case an error would be thrown elsewhere in the code.
   if ((0 <= xyz) && (xyz <= 2))
      return Kappa[comp] * pow_law_B * _fields.DelAbsMag()[xyz] / _fields.AbsMag();
   else
      return Kappa[comp] * pow_law_B * _fields.DotAbsMag() / _fields.AbsMag();
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
template <typename HConfig>
double DiffusionRigidityMagneticFieldPowerLaw<HConfig>::GetMuDerivative(Component comp)
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
template <typename HConfig>
DiffusionStraussEtAl2013<HConfig>::DiffusionStraussEtAl2013(void)
                        : DiffusionBase(diff_name, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/06/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename HConfig>
DiffusionStraussEtAl2013<HConfig>::DiffusionStraussEtAl2013(const DiffusionStraussEtAl2013& other)
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
template <typename HConfig>
void DiffusionStraussEtAl2013<HConfig>::SetupDiffusion(bool construct)
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
template <typename HConfig>
void DiffusionStraussEtAl2013<HConfig>::EvaluateDiffusion(Component comp)
{
   if (comp == Component::mu) return;

// Find LISM indicator variable (convert -1:1 to 1:0) and interpolate inner/outer quantities.
   if (LISM_idx < 0) LISM_ind = 0.0;
   else LISM_ind = (_fields.IvLISM() > 0.0 ? 0.0 : 1.0);
   double lam_para = LISM_ind * lam_out + (1.0 - LISM_ind) * lam_in;
   double B0_eff = LISM_ind * _fields.Mag() + (1.0 - LISM_ind) * B0;
   double rig = Rigidity<specie>(_coords.AbsMom());
   double kap_rat;

// Find diffusion coefficients
   Kappa[Component::para] = (lam_para * _coords.AbsVel() / 3.0) * (rig < R0 ? cbrt(rig / R0) : rig / R0) * (B0_eff / _fields.Mag());
   if (comp == Component::perp) {
      kap_rat = LISM_ind * kap_rat_out + (1.0 - LISM_ind) * kap_rat_in;
      Kappa[Component::perp] = kap_rat * Kappa[Component::para];
   };
};

/*!
\author Juan G Alonso Guzman
\date 05/13/2024
\return double       Derivative in mu
*/
template <typename HConfig>
double DiffusionStraussEtAl2013<HConfig>::GetMuDerivative(Component comp)
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
template <typename HConfig>
DiffusionPotgieterEtAl2015<HConfig>::DiffusionPotgieterEtAl2015(void)
                          : DiffusionBase(diff_name, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename HConfig>
DiffusionPotgieterEtAl2015<HConfig>::DiffusionPotgieterEtAl2015(const DiffusionPotgieterEtAl2015& other)
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
template <typename HConfig>
void DiffusionPotgieterEtAl2015<HConfig>::SetupDiffusion(bool construct)
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
template <typename HConfig>
void DiffusionPotgieterEtAl2015<HConfig>::EvaluateDiffusion(Component comp)
{
   if (comp == Component::mu) return;

// Find LISM indicator variable (convert -1:1 to 1:0) and interpolate inner/outer quantities.
   if (LISM_idx < 0) LISM_ind = 0.0;
   else LISM_ind = (_fields.IvLISM() > 0.0 ? 0.0 : 1.0);
   double kappa_para = LISM_ind * kappa_out + (1.0 - LISM_ind) * kappa_in;
   double B0_eff = LISM_ind * _fields.Mag() + (1.0 - LISM_ind) * B0;
   double rig = Rigidity<specie>(_coords.AbsMom());
   double kap_rat;

// Find diffusion coefficients
   Kappa[Component::para] = kappa_para * (_coords.AbsVel() / c_code) * (rig < R0 ? 1.0 : sqrt(Cube(rig / R0))) * (B0_eff / _fields.Mag());
   if (comp == Component::perp) {
      kap_rat = LISM_ind * kap_rat_out + (1.0 - LISM_ind) * kap_rat_in;
      Kappa[Component::perp] = kap_rat * Kappa[Component::para];
   };
};

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
\return double       Derivative in mu
*/
template <typename HConfig>
double DiffusionPotgieterEtAl2015<HConfig>::GetMuDerivative(Component comp)
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
template <typename HConfig>
DiffusionEmpiricalSOQLTandUNLT<HConfig>::DiffusionEmpiricalSOQLTandUNLT(void)
                              : DiffusionBase(diff_name, STATE_NONE)
{
};

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
\param[in] other Object to initialize from
A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename HConfig>
DiffusionEmpiricalSOQLTandUNLT<HConfig>::DiffusionEmpiricalSOQLTandUNLT(const DiffusionEmpiricalSOQLTandUNLT& other)
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
template <typename HConfig>
void DiffusionEmpiricalSOQLTandUNLT<HConfig>::SetupDiffusion(bool construct)
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
template <typename HConfig>
void DiffusionEmpiricalSOQLTandUNLT<HConfig>::EvaluateDiffusion(Component comp)
{
   if (comp == Component::mu) return;
   double lam, rig_dep;
   double rig = Rigidity<specie>(_coords.AbsMom());

   if (comp == Component::para) {
// Compute mean free path and rigidity dependance with a bent power law
      rig_dep = cbrt((rig / R0) * (1.0 + Sqr(Sqr(rig / R0))));
      lam = lam_para;
   }
   else if (comp == Component::perp) {
// Compute mean free path and rigidity dependance with a bent power law
      rig_dep = cbrt((rig / R0) * (1.0 + Sqr(rig / R0)));

// Find magnetic mixing indicator variable (convert -1:1 to 0:1) and interpolate perp-to-para diffusion ratio.
      if (Bmix_idx < 0) Bmix_ind = 1.0;
      Bmix_ind = (_fields.IvBmix() < 0.0 ? 0.0 : 1.0);
      if (_coords.Rad() < radial_limit_perp_red) lam = lam_perp * (Bmix_ind + (1.0 - Bmix_ind) * kap_rat_red);
      else lam = lam_perp;
   };
   Kappa[comp] = (lam * _coords.AbsVel() / 3.0) * rig_dep * (B0 / _fields.Mag());
   Kappa[comp] /= 1.0 + solar_cycle_effect * Sqr(cos(0.5 * _fields.IvSolarCycle()));
};

/*!
\author Juan G Alonso Guzman
\date 01/09/2025
\return double       Derivative in mu
*/
template <typename HConfig>
double DiffusionEmpiricalSOQLTandUNLT<HConfig>::GetMuDerivative(Component comp)
{
   return 0.0;
};

};
