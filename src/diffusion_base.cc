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
DiffusionBase::DiffusionBase(void)
             : Params("", 0, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 05/06/2022
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
DiffusionBase::DiffusionBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
             : Params(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 05/06/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
DiffusionBase::DiffusionBase(const DiffusionBase& other)
             : Params(other)
{
// Params' constructor resets all flags
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDiffusion(true);
};

/*!
\author Vladimir Florinski
\date 05/06/2022
\param[in] cont_in Container with parameters

This is the default method to set up an object. It should only be defined in the base class (XXXXBase). Derived classes should _not_ modify it! This version always calls the correct virtual "SetupDiffusion()" method.
*/
void DiffusionBase::SetupObject(const DataContainer& cont_in)
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
void DiffusionBase::SetupDiffusion(bool construct)
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
void DiffusionBase::EvaluateDiffusion(void)
{
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\date 07/06/2022
\param[in] comp      Which component to evaluate
\param[in] t_in      Time
\param[in] pos_in    Position
\param[in] mom_in    Momentum (p,mu,phi) coordinates
\param[in] spdata_in Spatial data at the required location
\return One diffusion component
\note This is a common routine that the derived classes should not change.
*/
double DiffusionBase::GetComponent(int comp, double t_in, const GeoVector& pos_in, const GeoVector& mom_in, const SpatialData& spdata_in)
{
   SetState(t_in, pos_in, mom_in);
   vmag = Vel(_mom[0], specie);
   _spdata._mask = spdata_in._mask;
   _spdata = spdata_in;
   Omega = CyclotronFrequency(vmag, _spdata.Bmag, specie);

#if TRAJ_TYPE != TRAJ_PARKER
   mu = _mom[1];
   st2 = 1.0 - Sqr(mu);
#endif

   comp_eval = comp;
   EvaluateDiffusion();
   return Kappa[comp];
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 07/12/2024
\param[in] xyz Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\return Directional derivative
\note This is meant to be called after GetComponent() for the componenent for which the derivative is wanted
*/
double DiffusionBase::GetDirectionalDerivative(int xyz)
{
   double _t_saved, Bmag_saved, derivative, _dr, _dt;
   GeoVector _pos_saved, Bvec_saved, Kappa_saved, Kappa_forw, Kappa_back;

// Save diffusion and field values at "current" position
   Kappa_saved = Kappa;
   Bvec_saved = _spdata.Bvec;
   Bmag_saved = _spdata.Bmag;

// Spatial derivatives
   if(0 <= xyz && xyz <= 2) {
// Save position, compute increment
      _pos_saved = _pos;

//This computation of "Bvec" at a displaced position is exact if numerical derivatives are used, and a good estimate if "_dr" is small enough.
      _dr = 0.5 * _spdata._dr[xyz];
      if(_spdata._dr_forw_fail[xyz]) Kappa_forw[comp_eval] = Kappa_saved[comp_eval];
      else {
         _pos[xyz] += _dr;
         _spdata.Bvec += _spdata.gradBvec.row[xyz] * _dr;
         _spdata.Bmag += _spdata.gradBmag[xyz] * _dr;
         EvaluateDiffusion();
         Kappa_forw[comp_eval] = Kappa[comp_eval];
      };
      _dr *= 2.0;
      if(_spdata._dr_back_fail[xyz]) Kappa_back[comp_eval] = Kappa_saved[comp_eval];
      else {
         _pos[xyz] -= _dr;
         _spdata.Bvec -= _spdata.gradBvec.row[xyz] * _dr;
         _spdata.Bmag -= _spdata.gradBmag[xyz] * _dr;
         EvaluateDiffusion();
         Kappa_back[comp_eval] = Kappa[comp_eval];
      };
      if(_spdata._dr_forw_fail[xyz] || _spdata._dr_back_fail[xyz]) _dr *= 0.5;

      derivative = (Kappa_forw[comp_eval] - Kappa_back[comp_eval]) / _dr;

// Restore position
      _pos = _pos_saved;
   }
// Time derivatives
   else {
// Save time, compute increment
      _t_saved = _t;

//A similar comment as the one in the spatial derivatives applies here for "_dt".
      _dt = 0.5 * _spdata._dt;
      if(_spdata._dt_forw_fail) {
         _t += _dt;
         _spdata.Bvec += _spdata.dBvecdt * _dt;
         _spdata.Bmag += _spdata.dBmagdt * _dt;
         EvaluateDiffusion();
         Kappa_forw[comp_eval] = Kappa[comp_eval];
      }
      else Kappa_forw[comp_eval] = Kappa_saved[comp_eval];
      _t += 2.0;
      if(_spdata._dt_back_fail) {
         _t -= _dt;
         _spdata.Bvec -= _spdata.dBvecdt * _dt;
         _spdata.Bmag -= _spdata.dBmagdt * _dt;
         EvaluateDiffusion();
         Kappa_back[comp_eval] = Kappa[comp_eval];
      }
      else Kappa_back[comp_eval] = Kappa_saved[comp_eval];
      if(_spdata._dt_forw_fail || _spdata._dt_back_fail) _dt *= 0.5;

      derivative = (Kappa_forw[comp_eval] - Kappa_back[comp_eval]) / _dt;

// Restore position
      _t = _t_saved;
   };

// Restore diffusion and field values at "current" position
   Kappa = Kappa_saved;
   _spdata.Bvec = Bvec_saved;
   _spdata.Bmag = Bmag_saved;

   return derivative;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinsk
\date 10/19/2022
\return Derivative in mu
\note This is meant to be called after GetComponent() for the componenent for which the derivative is wanted
*/
double DiffusionBase::GetMuDerivative(void)
{
   double mu_saved, dmu, derivative;
   GeoVector Kappa_saved;

// Save diffusion and field values at "current" position
   Kappa_saved = Kappa;
   mu_saved = mu;

// Mu derivative (momentum is in (p,mu,phi) coordinates)
   dmu = sp_small * (mu + sp_small < 1.0 ? 1.0 : -1.0);
#if TRAJ_TYPE != TRAJ_PARKER
   mu += dmu;
   st2 = 1.0 - Sqr(mu);
#endif
   EvaluateDiffusion();
   derivative = (Kappa[comp_eval] - Kappa_saved[comp_eval]) / dmu;

// Restore diffusion and field values at "current" position
   Kappa = Kappa_saved;
#if TRAJ_TYPE != TRAJ_PARKER
   mu = mu_saved;
   st2 = 1.0 - Sqr(mu);
#endif

   return derivative;
};

};
