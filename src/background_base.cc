/*!
\file background_base.cc
\brief Implements a base class to compute the plasma background
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_base.hh"
#include <stdexcept>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 11/24/2020
*/
template <typename Fields>
BackgroundBase<Fields>::BackgroundBase(void)
              : Params("", 0, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 12/17/2020
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
template <typename Fields>
BackgroundBase<Fields>::BackgroundBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
              : Params(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 09/26/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename Fields>
BackgroundBase<Fields>::BackgroundBase(const BackgroundBase& other)
              : Params(other)
{
// Params' constructor resets all flags
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\date 12/14/2020
\return Maximum distance based on the grid or other properties
*/
template <typename Fields>
double BackgroundBase<Fields>::GetDmax(void) const
{
   return _ddata.dmax;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/05/2025
\param[in] int index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\param[in] Fields temporary fields for case of status failure.
\note This is a common routine that the derived classes should not change.
*/
template <typename Fields>
void BackgroundBase<Fields>::DirectionalDerivative(int xyz, Fields& fields_tmp)
{
   double _t_saved;
   GeoVector _pos_saved;
// temporaries for forward- and backward-stepped evaluation
   Fields fields_forw, fields_back;

   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) {
      RAISE_BITS(_status, STATE_INVALID);
      throw ExUninitialized();
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Spatial derivatives
   if ((0 <= xyz) && (xyz <= 2)) {
      _ddata._dr_forw_fail[xyz] = false;
      _ddata._dr_back_fail[xyz] = false;

// Save position, compute increment      
      _pos_saved = _pos;

// Increment forward, then back.
// In either case, if the attempted state is invalid,
// revert to first order differences, otherwise use second order differences.

// Forward increment
      _ddata._dr[xyz] = r_g;
      _pos += _ddata._dr[xyz] * fa_basis.row[xyz];
      EvaluateBackground();

      if (BITS_LOWERED(_status, STATE_INVALID)) {
         fields_forw = _fields;
      }
      else {
         fields_forw = fields_tmp;
         _ddata._dr_forw_fail[xyz] = true;
      };

// Backward increment
      _ddata._dr[xyz] *= 2.0;
      _pos -= _ddata._dr[xyz] * fa_basis.row[xyz];
      EvaluateBackground();

      if (BITS_LOWERED(_status, STATE_INVALID)) {
         fields_back = _fields;
      }
      else {
         fields_back = fields_tmp;
         _ddata._dr_back_fail[xyz] = true;
      };

// If at least one increment failed, half the increment. If both increments failed, throw an error.
      if (_ddata._dr_forw_fail[xyz] || _ddata._dr_back_fail[xyz]) _ddata._dr[xyz] *= 0.5;
      if (_ddata._dr_forw_fail[xyz] && _ddata._dr_back_fail[xyz]) throw ExFieldError();

// Restore position
      _pos = _pos_saved;

// Compute spatial derivatives. This calculation gives gradV[i][j] = dV_j / ds^i which is the transpose of the Jacobian.
      if constexpr (Fields::DelVel_found())
         _fields.DelVel()[xyz] = (fields_forw.Vel() - fields_back.Vel()) / _ddata._dr[xyz];
      if constexpr (Fields::DelMag_found())
         _fields.DelMag()[xyz] = (fields_forw.Mag() - fields_back.Mag()) / _ddata._dr[xyz];
      if constexpr (Fields::DelAbsMag_found())
         _fields.DelAbsMag()[xyz] = (fields_forw.AbsMag() - fields_back.AbsMag()) / _ddata._dr[xyz];
      if constexpr (Fields::DelElc_found())
         _fields.DelElc()[xyz] = (fields_forw.Elc() - fields_back.Elc()) / _ddata._dr[xyz];
   }

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time derivatives
   else {
      _ddata._dt_forw_fail = false;
      _ddata._dt_back_fail = false;

// Save time, compute increment      
      _t_saved = _t;

// Increment forward, then back.
// In either case, if the attempted state is invalid,
// revert to first order differences, otherwise use second order differences.

// Forward increment
      _ddata._dt = 1.0 / w_g;
      _t += _ddata._dt;
      EvaluateBackground();

      if (BITS_LOWERED(_status, STATE_INVALID)) {
         fields_forw = _fields;
      }
      else {
         fields_forw = fields_tmp;
         _ddata._dt_forw_fail = true;
      };

// Backward increment
      _ddata._dt *= 2.0;
      _t -= _ddata._dt;
      EvaluateBackground();

      if (BITS_LOWERED(_status, STATE_INVALID)) {
         fields_back = _fields;
      }
      else {
         fields_back = fields_tmp;
         _ddata._dt_back_fail = true;
      };

// If at least one increment failed, half the increment. If both increments failed, throw an error.
      if (_ddata._dt_forw_fail || _ddata._dt_back_fail) _ddata._dt *= 0.5;
      if (_ddata._dt_forw_fail && _ddata._dt_back_fail) throw ExFieldError();

// Restore time
      _t = _t_saved;

// Compute time derivatives
      if constexpr (Fields::DdtVel_found())
         _fields.DdtVel() = (fields_forw.Vel() - fields_back.Vel()) / _ddata._dt;
      if constexpr (Fields::DdtMag_found())
         _fields.DdtMag() = (fields_forw.Mag() - fields_back.Mag()) / _ddata._dt;
      if constexpr (Fields::DdtAbsMag_found())
         _fields.DdtAbsMag() = (fields_forw.AbsMag() - fields_back.AbsMag()) / _ddata._dt;
      if constexpr (Fields::DdtElc_found())
         _fields.DdtElc() = (fields_forw.Elc() - fields_back.Elc()) / _ddata._dt;
   };
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/05/2025
*/
template <typename Fields>
void BackgroundBase<Fields>::NumericalDerivatives(void)
{

// copy _fields into temporary in case of averaging.
   Fields fields_tmp = _fields;
   double AbsMag;
   GeoVector HatMag;
   constexpr bool gradients = Fields::DelVel_found() || Fields::DelMag_found() || Fields::DelAbsMag_found() || Fields::DelElc_found();
   constexpr bool time_derivatives =  Fields::DdtVel_found() || Fields::DdtMag_found() || Fields::DdtAbsMag_found() || Fields::DdtElc_found();

// This method is only valid in the presence of a magnetic field.
// Normally magnetic field magnitude+direction is tracked but we can make do with magnetic field only.
   if constexpr (gradients || time_derivatives) {
      if constexpr (!Fields::AbsMag_found())
         AbsMag = _fields.Mag().Norm();
      else
         AbsMag = _fields.AbsMag();
      if constexpr (!Fields::HatMag_found())
         HatMag = _fields.Mag()/AbsMag;
      else
         HatMag = _fields.HatMag();
   }

   if constexpr (gradients) {

// Derivatives are only needed for trajectory types whose transport assumes the background changes on scales larger than the gyro-radius.
      r_g = fmin(LarmorRadius(_mom[0], AbsMag, specie), _ddata.dmax);


// Get field aligned basis in (transpose) of rotation matrix
      fa_basis.row[2] = HatMag;
      fa_basis.row[0] = GetSecondUnitVec(HatMag);
      fa_basis.row[1] = fa_basis.row[2] ^ fa_basis.row[0];

// Compute derivatives in field-aligned basis
      for (auto xyz = 0; xyz < 3; xyz++) DirectionalDerivative(xyz, fields_tmp);

// Transform basis back to global cartesian frame
      rot_mat.Transpose(fa_basis);
      if constexpr (Fields::DelVel_found())
         fields_tmp.DelVel() = rot_mat * _fields.DelVel();
      if constexpr (Fields::DelMag_found())
         fields_tmp.DelMag() = rot_mat * _fields.DelMag();
      if constexpr (Fields::DelAbsMag_found())
         fields_tmp.DelAbsMag() = rot_mat * _fields.DelAbsMag();
      if constexpr (Fields::DelElc_found())
         fields_tmp.DelElc() = rot_mat * _fields.DelElc();

#if BACKGROUND_NUM_GRAD_EVALS > 1

// Repeat for any additional rotations
      for (auto n_rot = 1; n_rot < BACKGROUND_NUM_GRAD_EVALS; n_rot++) {
         fa_basis[0].Rotate(fa_basis.row[2], sin_lra, cos_lra);
         fa_basis[1].Rotate(fa_basis.row[2], sin_lra, cos_lra);

         for (auto xyz = 0; xyz < 3; xyz++) DirectionalDerivative(xyz, fields_tmp);

         rot_mat.Transpose(fa_basis);
         if constexpr (Fields::DelVel_found())
            fields_tmp.DelVel() += rot_mat * _fields.DelVel();
         if constexpr (Fields::DelMag_found())
            fields_tmp.DelMag() += rot_mat * _fields.DelMag();
         if constexpr (Fields::DelAbsMag_found())
            fields_tmp.DelAbsMag() += rot_mat * _fields.DelAbsMag();
         if constexpr (Fields::DelElc_found())
            fields_tmp.DelElc() += rot_mat * _fields.DelElc();
      };

// Average results
      if constexpr (Fields::DelVel_found())
         fields_tmp.DelVel() /= BACKGROUND_NUM_GRAD_EVALS;
      if constexpr (Fields::DelMag_found())
         fields_tmp.DelMag() /= BACKGROUND_NUM_GRAD_EVALS;
      if constexpr (Fields::DelAbsMag_found())
         fields_tmp.DelAbsMag() /= BACKGROUND_NUM_GRAD_EVALS;
      if constexpr (Fields::DelElc_found())
         fields_tmp.DelElc() /= BACKGROUND_NUM_GRAD_EVALS;

#endif

   };

// Time derivatives.
   if constexpr (time_derivatives) {
// Derivatives are only needed for trajectory types whose transport assumes the background changes on scales longer than the gyro-frequency.
      w_g = fmin(CyclotronFrequency(Vel(_mom[0]), AbsMag, specie), Vel(_mom[0]) / _ddata.dmax);
      DirectionalDerivative(3, fields_tmp);
   };

// Copy data from temporary fields into _fields
   _fields = fields_tmp;
};

/*!
\author Vladimir Florinski
\date 11/25/2020
\param[in] cont_in Container with parameters

This is the default method to set up an object. It should only be defined in the base class (XXXXBase). Derived classes should _not_ modify it! This version always calls the correct virtual "SetupBackground()" method.
*/
template <typename Fields>
void BackgroundBase<Fields>::SetupObject(const DataContainer& cont_in)
{
   Params::SetContainer(cont_in);
   SetupBackground(false);
};

/*!
\author Vladimir Florinski
\date 09/26/2021
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename Fields>
void BackgroundBase<Fields>::SetupBackground(bool construct)
{
// Only needed in the parent version
   container.Reset();
   RAISE_BITS(_status, STATE_SETUP_COMPLETE);
   LOWER_BITS(_status, STATE_INVALID);

// Unpack parameters
   container.Read(t0);
   container.Read(r0);
   container.Read(u0);
   container.Read(B0);
   container.Read(dmax0);

// Initialize "safe" box for derivatives
   for (auto xyz = 0; xyz < 3; xyz++) _ddata._dr[xyz] = incr_dmax_ratio * dmax0;
   _ddata._dt = incr_dmax_ratio * dmax0 / c_code;
   _ddata.dmax = dmax0;
};

/*!
\author Vladimir Florinski
\date 12/08/2021
\note This is only a stub
*/
template <typename Fields>
void BackgroundBase<Fields>::EvaluateBackground(void)
{
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\date 10/13/2022
\note This is only a stub
*/
template <typename Fields>
void BackgroundBase<Fields>::EvaluateBackgroundDerivatives(void)
{
};

/*!
\author Vladimir Florinski
\date 09/15/2021
\note The default method should be good enough for all grid-free backgrounds
*/
template <typename Fields>
void BackgroundBase<Fields>::EvaluateDmax(void)
{
   _ddata.dmax = dmax0;
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/08/2025
\note The default method should be good enough for all grid-free backgrounds
*/
template <typename Fields>
void BackgroundBase<Fields>::EvaluateBmag(void)
{
   _fields.AbsMag() = _fields.Mag().Norm();
};

/*!
\author Vladimir Florinski
\date 02/17/2023
\note The default method should be good enough for all grid-free backgrounds
*/
template <typename Fields>
void BackgroundBase<Fields>::StopServerFront(void)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/19/2022
\param[in] dir Direction
\return Safe increment in some direction (potential negative) to stay inside domain
*/
template <typename Fields>
double BackgroundBase<Fields>::GetSafeIncr(const GeoVector& dir)
{
//FIXME: This is incomplete.
   return incr_dmax_ratio * _ddata.dmax;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/05/2025
\param[in]  t_in   Time
\param[in]  pos_in Position
\param[in]  mom_in Momentum (p,mu,phi) coordinates
\param[out] fields All fields data
\note This is a common routine that the derived classes should not change.
*/
template <typename Fields>
void BackgroundBase<Fields>::GetFields(double t_in, const GeoVector& pos_in, const GeoVector& mom_in, Fields& fields)
{
// When the call is finished, the `fields` member will be updated for the caller's state, and the `fields` argument will be equal to that member.

// Check that state setup is complete
   if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) {
      RAISE_BITS(_status, STATE_INVALID);
      throw ExUninitialized();
   };

// If "EvaluateDmax()" fails, the state will be set to "STATE_INVALID" and background will not be evaluated
   SetState(t_in, pos_in, mom_in);
   EvaluateDmax();
   if (BITS_RAISED(_status, STATE_INVALID)) throw ExCoordinates();

// Compute u, B, E
   EvaluateBackground();
   if (BITS_RAISED(_status, STATE_INVALID)) throw ExFieldError();
// Compute Bmag, bhat (always the same)
   if constexpr (Fields::Mag_found()) {
      double AbsMag;
      if constexpr (Fields::AbsMag_found()) {
         EvaluateBmag();
         AbsMag = _fields.AbsMag();
      } else {
         AbsMag = _fields.Mag().Norm();
      }
      if constexpr (Fields::HatMag_found()) {
         if (AbsMag < sp_tiny) RAISE_BITS(_status, STATE_INVALID);
         else _fields.AbsMag() = _fields.Mag() / AbsMag;
      };
   }
   if (BITS_RAISED(_status, STATE_INVALID)) throw ExFieldError();
// Compute derivatives of u, B, E
   EvaluateBackgroundDerivatives();

// Copy the internal fields into arguments
   fields = _fields;
};

};



