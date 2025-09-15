/*!
\file diffusion_base.cc
\brief Implements a base class to compute difusion coefficients
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "diffusion_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/06/2022
*/
template <typename Trajectory>
DiffusionBase<Trajectory>::DiffusionBase(void)
             : Params("", STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 05/06/2022
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
template <typename Trajectory>
DiffusionBase<Trajectory>::DiffusionBase(const std::string& name_in,  uint16_t status_in)
             : Params(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 05/06/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename Trajectory>
DiffusionBase<Trajectory>::DiffusionBase(const DiffusionBase& other)
             : Params(other)
{
// Params' constructor resets all flags
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Vladimir Florinski
\date 05/06/2022
\param[in] cont_in Container with parameters

This is the default method to set up an object. It should only be defined in the base class (XXXXBase). Derived classes should _not_ modify it! This version always calls the correct virtual "SetupDiffusion()" method.
*/
template <typename Trajectory>
void DiffusionBase<Trajectory>::SetupObject(const DataContainer& cont_in)
{
   Params::SetContainer(cont_in);
   SetupDiffusion(false);
};

/*!
\author Vladimir Florinski
\date 05/06/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename Trajectory>
void DiffusionBase<Trajectory>::SetupDiffusion(bool construct)
{
// Only needed in the parent version
   container.Reset();
   RAISE_BITS(_status, STATE_SETUP_COMPLETE);
   LOWER_BITS(_status, STATE_INVALID);

// Reset Kappa
   Kappa = gv_zeros;
};

/*!
\author Vladimir Florinski
\date 05/09/2022
\note This is only a stub.
*/
template <typename Trajectory>
void DiffusionBase<Trajectory>::EvaluateDiffusion(void)
{
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/10/2025
\param[in] comp      Which component to evaluate
\param[in] t_in      Time
\param[in] pos_in    Position
\param[in] mom_in    Momentum (p,mu,phi) coordinates
\param[in] fields_in Spatial fields at the required location
\return One diffusion component
\note This is a common routine that the derived classes should not change.
*/
template <typename Trajectory>
double DiffusionBase<Trajectory>::GetComponent(int comp, double t_in, const GeoVector& pos_in, const GeoVector& mom_in, const Fields& fields_in)
{
   SetState(t_in, pos_in, mom_in);
   vmag = Vel(_mom[0], Trajectory::specie);
   _fields = fields_in;
   Omega = CyclotronFrequency(vmag, _fields.AbsMag(), Trajectory::specie);

   if constexpr (!std::same_as<Trajectory, TrajectoryParker<Fields>>) {
      mu = _mom[1];
      st2 = 1.0 - Sqr(mu);
   }

   comp_eval = comp;
   EvaluateDiffusion();
   return Kappa[comp];
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 07/12/2024
\param[in] xyz Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\param[in] ddata_in Derivative data from computing background fields
\return Directional derivative
\note This is meant to be called after GetComponent() for the componenent for which the derivative is wanted
*/
template <typename Trajectory>
double DiffusionBase<Trajectory>::GetDirectionalDerivative(int xyz, DerivativeData& ddata_in)
{
   double _t_saved, Bmag_saved, derivative, _dr, _dt;
   GeoVector _pos_saved, Bvec_saved, Kappa_saved, Kappa_forw, Kappa_back;

// Save diffusion and field values at "current" position
   Kappa_saved = Kappa;
   Bvec_saved = _fields.Mag();
   Bmag_saved = _fields.AbsMag();

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Spatial derivatives
   if ((0 <= xyz) && (xyz <= 2)) {

// Save position, compute increment
      _pos_saved = _pos;

// This computation of "Bvec" at a displaced position is exact if numerical derivatives are used, and a good estimate if "_dr" is small enough.
      _dr = 0.5 * _ddata._dr[xyz];
      if (_ddata._dr_forw_fail[xyz]) Kappa_forw[comp_eval] = Kappa_saved[comp_eval];
      else {
         _pos[xyz] += _dr;
         _fields.Mag() += _fields.DelMag().row[xyz] * _dr;
         _fields.AbsMag() += _fields.DelAbsMag()[xyz] * _dr;
         EvaluateDiffusion();
         Kappa_forw[comp_eval] = Kappa[comp_eval];
      };

      _dr *= 2.0;
      if (_ddata._dr_back_fail[xyz]) Kappa_back[comp_eval] = Kappa_saved[comp_eval];
      else {
         _pos[xyz] -= _dr;
         _fields.Mag() -= _fields.DelMag().row[xyz] * _dr;
         _fields.AbsMag() -= _fields.DelAbsMag()[xyz] * _dr;
         EvaluateDiffusion();
         Kappa_back[comp_eval] = Kappa[comp_eval];
      };

      if (_ddata._dr_forw_fail[xyz] || _ddata._dr_back_fail[xyz]) _dr *= 0.5;
      derivative = (Kappa_forw[comp_eval] - Kappa_back[comp_eval]) / _dr;

// Restore position
      _pos = _pos_saved;
   }

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time derivatives
   else {

// Save time, compute increment
      _t_saved = _t;

//A similar comment as the one in the spatial derivatives applies here for "_dt".
      _dt = 0.5 * _ddata._dt;
      if (_ddata._dt_forw_fail) {
         _t += _dt;
         _fields.Mag() += _fields.DotMag() * _dt;
         _fields.AbsMag() += _fields.DotAbsMag() * _dt;
         EvaluateDiffusion();
         Kappa_forw[comp_eval] = Kappa[comp_eval];
      }
      else Kappa_forw[comp_eval] = Kappa_saved[comp_eval];

      _t += 2.0;
      if (_ddata._dt_back_fail) {
         _t -= _dt;
         _fields.Mag() -= _fields.DotMag() * _dt;
         _fields.AbsMag() -= _fields.DotAbsMag() * _dt;
         EvaluateDiffusion();
         Kappa_back[comp_eval] = Kappa[comp_eval];
      }
      else Kappa_back[comp_eval] = Kappa_saved[comp_eval];

      if (_ddata._dt_forw_fail || _ddata._dt_back_fail) _dt *= 0.5;
      derivative = (Kappa_forw[comp_eval] - Kappa_back[comp_eval]) / _dt;

// Restore time
      _t = _t_saved;
   };

// Restore diffusion and field values at "current" position
   Kappa = Kappa_saved;
   _fields.Mag() = Bvec_saved;
   _fields.AbsMag() = Bmag_saved;

   return derivative;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/10/2025
\return Derivative in mu
\note This is meant to be called after GetComponent() for the componenent for which the derivative is wanted
*/
template <typename Trajectory>
double DiffusionBase<Trajectory>::GetMuDerivative(void)
{
   double mu_saved, dmu, derivative;
   GeoVector Kappa_saved;

// Save diffusion and field values at "current" position
   Kappa_saved = Kappa;
   mu_saved = mu;

// Mu derivative (momentum is in (p,mu,phi) coordinates)
   dmu = sp_small * (mu + sp_small < 1.0 ? 1.0 : -1.0);

   if constexpr (!std::same_as<Trajectory, TrajectoryParker<Fields>>) {
      mu += dmu;
      st2 = 1.0 - Sqr(mu);
   }

   EvaluateDiffusion();
   derivative = (Kappa[comp_eval] - Kappa_saved[comp_eval]) / dmu;

// Restore diffusion and field values at "current" position
   Kappa = Kappa_saved;

   if constexpr (!std::same_as<Trajectory, TrajectoryParker<Fields>>) {
      mu = mu_saved;
      st2 = 1.0 - Sqr(mu);
   }

   return derivative;
};

};
