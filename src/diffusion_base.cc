/*!
\file diffusion_base.cc
\brief Implements a base class to compute difusion coefficients
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "diffusion_base.hh"

namespace Spectrum {

using namespace DiffusionOptions;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DiffusionBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/06/2022
*/
template <typename HConfig>
DiffusionBase<HConfig>::DiffusionBase(void)
             : Params("", STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 05/06/2022
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
DiffusionBase<HConfig>::DiffusionBase(const std::string& name_in,  uint16_t status_in)
             : Params(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 05/06/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupDiffusion()" with the argument of "true".
*/
template <typename HConfig>
DiffusionBase<HConfig>::DiffusionBase(const DiffusionBase& other)
             : Params(other)
{
// Params' constructor resets all flags
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) DiffusionBase::SetupDiffusion(true);
};

/*!
\author Vladimir Florinski
\date 05/06/2022
\param[in] cont_in Container with parameters

This is the default method to set up an object. It should only be defined in the base class (XXXXBase). Derived classes should _not_ modify it! This version always calls the correct virtual "SetupDiffusion()" method.
*/
template <typename HConfig>
void DiffusionBase<HConfig>::SetupObject(const DataContainer& cont_in)
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
template <typename HConfig>
void DiffusionBase<HConfig>::SetupDiffusion(bool construct)
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
template <typename HConfig>
void DiffusionBase<HConfig>::EvaluateDiffusion(Component comp)
{
   LOWER_BITS(_status, STATE_INVALID);
};


/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/21/2025
\param[in] coords     Coordinates managed by trajectory type
\param[in] fields       Coordinate-dependent fields managed by trajectory type
\note This is a common routine that the derived classes should not change.
This method must be called prior to all diffusion computations on 'new' coordinates.
After that, methods Get, GetDirectionalDerivative, and GetMuDerivative
can be called any number of times. To use a particular diffusion type,
the trajectory fields and coordinates must be convertible to that diffusion
type's fields and coordinates. Conversion is performed by the caller to Stage()
because situations exist where the coordinates are most easily adjusted after conversion, not before.
Note that fields passed by the caller to Stage() are *inputs* to diffusion class,
that can only be evaluated once coordinates are known, so here, too, the
caller should make use of the Diffusion class's fields type and populate it with the values
required by the diffusion coefficient computation. This is normally done by calling the background,
but it is not necessarily the case that the fields need to be computed at the point of staging,
so the caller is again given the responsibility to prepare the argument.
*/
template <typename HConfig>
void DiffusionBase<HConfig>::Stage(const DiffusionCoordinates& coords, const DiffusionFields& fields)
{
   _coords = coords;
   _fields = fields;
   // todo this is in Convert
   // todo MomMu, VelMu, MomSt2, VelSt2
//   vmag = Vel<specie>(coords.AbsMom());
   Omega = CyclotronFrequency<specie>(_coords.AbsVel(), _fields.AbsMag());
   // todo review
   st2 = 1.0 - Sqr(_coords.MomMu());
};



/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/21/2025
\param[in] comp      Which component to evaluate
\return One diffusion component
\note This is a common routine that the derived classes should not change.
\note This must be called after Stage() if the target coordinates have changed.
*/
template <typename HConfig>
double DiffusionBase<HConfig>::Get(Component comp)
{
   EvaluateDiffusion(comp);
   return Kappa[comp];
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/12/2025
\param[in] comp Which component to evaluate
\param[in] xyz Index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\param[in] ddata_in Derivative data from computing background fields (read only)
\return Directional derivative
\note This must be called after Stage() if the target coordinates have changed.
*/
template <typename HConfig>
double DiffusionBase<HConfig>::GetDirectionalDerivative(Component comp, int xyz, const DerivativeData& ddata)
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
      _pos_saved = _coords.Pos();

// This computation of "Bvec" at a displaced position is exact if numerical derivatives are used, and a good estimate if "_dr" is small enough.
      _dr = 0.5 * ddata._dr[xyz];
      if (ddata._dr_forw_fail[xyz]) Kappa_forw[comp] = Kappa_saved[comp];
      else {
         _coords.Pos()[xyz] += _dr;
         _fields.Mag() += _fields.DelMag().row[xyz] * _dr;
         _fields.AbsMag() += _fields.DelAbsMag()[xyz] * _dr;
         EvaluateDiffusion(comp);
         Kappa_forw[comp] = Kappa[comp];
      };

      _dr *= 2.0;
      if (ddata._dr_back_fail[xyz]) Kappa_back[comp] = Kappa_saved[comp];
      else {
         _coords.Pos()[xyz] -= _dr;
         _fields.Mag() -= _fields.DelMag().row[xyz] * _dr;
         _fields.AbsMag() -= _fields.DelAbsMag()[xyz] * _dr;
         EvaluateDiffusion(comp);
         Kappa_back[comp] = Kappa[comp];
      };

      if (ddata._dr_forw_fail[xyz] || ddata._dr_back_fail[xyz]) _dr *= 0.5;
      derivative = (Kappa_forw[comp] - Kappa_back[comp]) / _dr;

// Restore position
      _coords.Pos() = _pos_saved;
   }

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time derivatives
   else {

// Save time, compute increment
      _t_saved = _coords.Time();

//A similar comment as the one in the spatial derivatives applies here for "_dt".
      _dt = 0.5 * ddata._dt;
      if (ddata._dt_forw_fail) {
         _coords.Time() += _dt;
         _fields.Mag() += _fields.DotMag() * _dt;
         _fields.AbsMag() += _fields.DotAbsMag() * _dt;
         EvaluateDiffusion(comp);
         Kappa_forw[comp] = Kappa[comp];
      }
      else Kappa_forw[comp] = Kappa_saved[comp];

      _coords.Time() += 2.0;
      if (ddata._dt_back_fail) {
         _coords.Time() -= _dt;
         _fields.Mag() -= _fields.DotMag() * _dt;
         _fields.AbsMag() -= _fields.DotAbsMag() * _dt;
         EvaluateDiffusion(comp);
         Kappa_back[comp] = Kappa[comp];
      }
      else Kappa_back[comp] = Kappa_saved[comp];

      if (ddata._dt_forw_fail || ddata._dt_back_fail) _dt *= 0.5;
      derivative = (Kappa_forw[comp] - Kappa_back[comp]) / _dt;

// Restore time
      _coords.Time() = _t_saved;
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
\note This must be called after Stage() if the target coordinates have changed.
*/
template <typename HConfig>
double DiffusionBase<HConfig>::GetMuDerivative(Component comp)
{
// Save diffusion and field values at "current" position
   double mu_saved = _coords.MomMu();
   GeoVector Kappa_saved = Kappa;
   double derivative;
   double dmu = sp_small * (mu_saved + sp_small < 1.0 ? 1.0 : -1.0);

   _coords.MomMu() += dmu;
   st2 = 1.0 - Sqr(_coords.MomMu());

   EvaluateDiffusion(comp);
   derivative = (Kappa[comp] - Kappa_saved[comp]) / dmu;

// Restore diffusion and field values at "current" position
   Kappa = Kappa_saved;

   _coords.MomMu() = mu_saved;
   st2 = 1.0 - Sqr(mu_saved);

   return derivative;
};

};


