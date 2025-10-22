/*!
\file background_base.cc
\brief Implements a base class to compute the plasma background
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_base.hh"
#include <stdexcept>

namespace Spectrum {

using namespace BackgroundOptions;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 11/24/2020
*/
template <typename HConfig>
BackgroundBase<HConfig>::BackgroundBase(void)
              : Params("", STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 12/17/2020
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
BackgroundBase<HConfig>::BackgroundBase(const std::string& name_in, uint16_t status_in)
              : Params(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 09/26/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
template <typename HConfig>
BackgroundBase<HConfig>::BackgroundBase(const BackgroundBase& other)
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
template <typename HConfig>
double BackgroundBase<HConfig>::GetDmax(void) const
{
   return _ddata.dmax;
};

/*!
\author Lucius Schoenbaum
\date 08/15/2025
\return DerivativeData, a small struct containing information about the most recent derivative computed, and dmax
\note This information from the background is currently only needed by Diffusion classes
when numerical directional derivatives are computed.
 */
template <typename HConfig>
DerivativeData BackgroundBase<HConfig>::GetDerivativeData(void) const
{
   return _ddata;
}

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/05/2025
\param[in] int index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\param[in] coords coordinates requesting the computation (copy consumed by computation)
\param[out] fields requested fields (writeout argument)
\note This is a common routine that the derived classes should not change.
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
void BackgroundBase<HConfig>::DirectionalDerivative(const int xyz, Coordinates coords, Fields& fields)
{
// temporaries for forward- and backward-stepped evaluation
   Fields fields_forw, fields_back;

   if constexpr (HConfig::debug) {
      if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) {
         RAISE_BITS(_status, STATE_INVALID);
         throw ExUninitialized();
      };
   };

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Spatial derivatives
   if ((0 <= xyz) && (xyz <= 2)) {

      _ddata._dr_forw_fail[xyz] = false;
      _ddata._dr_back_fail[xyz] = false;

// Increment forward, then back.
// In either case, if the attempted state is invalid,
// revert to first order differences, otherwise use second order differences.

// Forward increment
      _ddata._dr[xyz] = r_g;
      coords.Pos() += _ddata._dr[xyz] * fa_basis.row[xyz];
      EvaluateBackground<Coordinates, Fields, RequestedFields>(coords, fields_forw);

      if (BITS_RAISED(_status, STATE_INVALID)) {
         fields_forw = fields;
         _ddata._dr_forw_fail[xyz] = true;
      };

// Backward increment
      _ddata._dr[xyz] *= 2.0;
      coords.Pos() -= _ddata._dr[xyz] * fa_basis.row[xyz];
      EvaluateBackground<Coordinates, Fields, RequestedFields>(coords, fields_back);

      if (BITS_RAISED(_status, STATE_INVALID)) {
         fields_back = fields;
         _ddata._dr_back_fail[xyz] = true;
      };

// If at least one increment failed, half the increment. If both increments failed, throw an error.
      if (_ddata._dr_forw_fail[xyz] || _ddata._dr_back_fail[xyz]) _ddata._dr[xyz] *= 0.5;
      if (_ddata._dr_forw_fail[xyz] && _ddata._dr_back_fail[xyz]) throw ExFieldError();

// Compute spatial derivatives. This calculation gives gradV[i][j] = dV_j / ds^i which is the transpose of the Jacobian.
      if constexpr (RequestedFields::DelFluv_found())
         fields.DelFluv()[xyz] = (fields_forw.Fluv() - fields_back.Fluv()) / _ddata._dr[xyz];
      if constexpr (RequestedFields::DelMag_found())
         fields.DelMag()[xyz] = (fields_forw.Mag() - fields_back.Mag()) / _ddata._dr[xyz];
      if constexpr (RequestedFields::DelAbsMag_found())
         fields.DelAbsMag()[xyz] = (fields_forw.AbsMag() - fields_back.AbsMag()) / _ddata._dr[xyz];
      if constexpr (RequestedFields::DelElc_found())
         fields.DelElc()[xyz] = (fields_forw.Elc() - fields_back.Elc()) / _ddata._dr[xyz];
   }

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time derivatives
   else {
      _ddata._dt_forw_fail = false;
      _ddata._dt_back_fail = false;

// Increment forward, then back.
// In either case, if the attempted state is invalid,
// revert to first order differences, otherwise use second order differences.

// Forward increment
      _ddata._dt = 1.0 / w_g;
      coords.Time() += _ddata._dt;
      EvaluateBackground<Coordinates, Fields, RequestedFields>(coords, fields_forw);

      if (BITS_RAISED(_status, STATE_INVALID)) {
         fields_forw = fields;
         _ddata._dt_forw_fail = true;
      };

// Backward increment
      _ddata._dt *= 2.0;
      coords.Time() -= _ddata._dt;
      EvaluateBackground<Coordinates, Fields, RequestedFields>(coords, fields_back);

      if (BITS_LOWERED(_status, STATE_INVALID)) {
         fields_back = fields;
         _ddata._dt_back_fail = true;
      };

// If at least one increment failed, half the increment. If both increments failed, throw an error.
      if (_ddata._dt_forw_fail || _ddata._dt_back_fail) _ddata._dt *= 0.5;
      if (_ddata._dt_forw_fail && _ddata._dt_back_fail) throw ExFieldError();

// Compute time derivatives
      if constexpr (RequestedFields::DotFluv_found())
         fields.DotFluv() = (fields_forw.Fluv() - fields_back.Fluv()) / _ddata._dt;
      if constexpr (RequestedFields::DotMag_found())
         fields.DotMag() = (fields_forw.Mag() - fields_back.Mag()) / _ddata._dt;
      if constexpr (RequestedFields::DotAbsMag_found())
         fields.DotAbsMag() = (fields_forw.AbsMag() - fields_back.AbsMag()) / _ddata._dt;
      if constexpr (RequestedFields::DotElc_found())
         fields.DotElc() = (fields_forw.Elc() - fields_back.Elc()) / _ddata._dt;
   };
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 09/08/2025
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
void BackgroundBase<HConfig>::NumericalDerivatives(Coordinates& coords, Fields& fields)
{
   double AbsMom = fields.AbsMom();
   double AbsMag = fields.AbsMag();
   GeoVector HatMag = fields.HatMag();
   constexpr Specie specie = HConfig::specie;
   constexpr bool gradients = RequestedFields::DelFluv_found() || RequestedFields::DelMag_found() || RequestedFields::DelAbsMag_found() || RequestedFields::DelElc_found();
   constexpr bool time_derivatives =  RequestedFields::DotFluv_found() || RequestedFields::DotMag_found() || RequestedFields::DotAbsMag_found() || RequestedFields::DotElc_found();

// This method is only valid in the presence of a magnetic field.
// Normally magnetic field magnitude+direction is tracked but we can make do with magnetic field only.

   if constexpr (gradients) {

// Derivatives are only needed for trajectory types whose transport assumes the background changes on scales larger than the gyro-radius.
      r_g = fmin(LarmorRadius<specie>(AbsMom, AbsMag), _ddata.dmax);

// Get field aligned basis in (transpose) of rotation matrix
      fa_basis.row[2] = HatMag;
      fa_basis.row[0] = GetSecondUnitVec(HatMag);
      fa_basis.row[1] = fa_basis.row[2] ^ fa_basis.row[0];

// Compute derivatives in field-aligned basis
      for (auto xyz = 0; xyz < 3; xyz++)
         DirectionalDerivative<Coordinates, Fields, RequestedFields>(xyz, coords, fields);

// Transform basis back to global cartesian frame
      rot_mat.Transpose(fa_basis);
      if constexpr (RequestedFields::DelFluv_found())
         fields.DelFluv() = rot_mat * fields.DelFluv();
      if constexpr (RequestedFields::DelMag_found())
         fields.DelMag() = rot_mat * fields.DelMag();
      if constexpr (RequestedFields::DelAbsMag_found())
         fields.DelAbsMag() = rot_mat * fields.DelAbsMag();
      if constexpr (RequestedFields::DelElc_found())
         fields.DelElc() = rot_mat * fields.DelElc();

      if constexpr (num_numeric_grad_evals > 1) {
         constexpr int num_evals = num_numeric_grad_evals;
         // copy fields into temporary
         Fields fields_tmp = fields;

// Repeat for each additional rotation
         for (auto n_rot = 1; n_rot < num_evals; ++n_rot) {
            fa_basis[0].Rotate(fa_basis.row[2], sin_lra, cos_lra);
            fa_basis[1].Rotate(fa_basis.row[2], sin_lra, cos_lra);

            for (auto xyz = 0; xyz < 3; xyz++) DirectionalDerivative<Coordinates, Fields, RequestedFields>(xyz, coords, fields_tmp);

            rot_mat.Transpose(fa_basis);
            if constexpr (RequestedFields::DelFluv_found())
               fields.DelFluv() += rot_mat * fields_tmp.DelFluv();
            if constexpr (RequestedFields::DelMag_found())
               fields.DelMag() += rot_mat * fields_tmp.DelMag();
            if constexpr (RequestedFields::DelAbsMag_found())
               fields.DelAbsMag() += rot_mat * fields_tmp.DelAbsMag();
            if constexpr (RequestedFields::DelElc_found())
               fields.DelElc() += rot_mat * fields_tmp.DelElc();

         };

// Average results
         if constexpr (RequestedFields::DelFluv_found())
            fields.DelFluv() /= num_evals;
         if constexpr (RequestedFields::DelMag_found())
            fields.DelMag() /= num_evals;
         if constexpr (RequestedFields::DelAbsMag_found())
            fields.DelAbsMag() /= num_evals;
         if constexpr (RequestedFields::DelElc_found())
            fields.DelElc() /= num_evals;
      };
   };

// Time derivatives.
   if constexpr (time_derivatives) {
// Derivatives are only needed for trajectory types whose transport assumes the background changes on scales longer than the gyro-frequency.
      auto vel = Vel<specie>(AbsMom);
      w_g = fmin(CyclotronFrequency(vel, AbsMag), vel / _ddata.dmax);
      DirectionalDerivative<Coordinates, Fields, RequestedFields>(3, coords, fields);
   };

};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/05/2025
\param[in] cont_in Container with parameters

This is the default method to set up an object. It should only be defined in the base class (XXXXBase). Derived classes should _not_ modify it! This version always calls the correct virtual "SetupBackground()" method.
*/
template <typename HConfig>
void BackgroundBase<HConfig>::SetupObject(const DataContainer& cont_in)
{
   Params::SetContainer(cont_in);
   SetupBackground(false);
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/05/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void BackgroundBase<HConfig>::SetupBackground(bool construct)
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
\author Lucius Schoenbaum
\date 08/05/2025
\note This is only a stub
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
void BackgroundBase<HConfig>::EvaluateBackground(Coordinates& coords, Fields& fields)
{
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/30/2025
\note This is a stub to be filled in by particular backgrounds.
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
void BackgroundBase<HConfig>::EvaluateBackgroundDerivatives(Coordinates& coords, Fields& fields)
{
   NumericalDerivatives<Coordinates, Fields, RequestedFields>(coords, fields);
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/09/2025
\note The default method should be good enough for all grid-free backgrounds
*/
template <typename HConfig>
template <typename Coordinates>
void BackgroundBase<HConfig>::EvaluateDmax(Coordinates& coords)
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
template <typename HConfig>
template <typename Fields>
void BackgroundBase<HConfig>::EvaluateAbsMag(Fields& fields)
{
   return fields.AbsMag();
};

/*!
\author Vladimir Florinski
\date 02/17/2023
\note The default method should be good enough for all grid-free backgrounds
*/
template <typename HConfig>
void BackgroundBase<HConfig>::StopServerFront(void)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/19/2022
\param[in] dir Direction
\return Safe increment in some direction (potential negative) to stay inside domain
*/
template <typename HConfig>
double BackgroundBase<HConfig>::GetSafeIncr(const GeoVector& dir)
{
//FIXME: This is incomplete.
   return incr_dmax_ratio * _ddata.dmax;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 09/08/2025
\param[in] coords coordinates, any coordinate system providing (time, position, p*) with p* the magnitude of momentum (access via AbsMom), and position in cartesian system.
\param[out] fields All fields data requested by caller. Optional type RequestedFields specifies which fields to evaluate, if only a subset is needed.
\note This is a common routine that the derived classes should not change.
\note This public method is valid for any Coordinate type that includes
position and time (Pos_t, Time_t) and magnitude of momentum (AbsMom_t or at least Mom_t).
The fields type must always contain magnetic field (Mag_t).
Magnetic field magnitude and/or direction can also be tracked but magnetic field is sufficient.
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
void BackgroundBase<HConfig>::GetFields(Coordinates& coords, Fields& fields)
{

// Check that state setup is complete
   if (HConfig::build_mode == BuildMode::debug) {
      if (BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) {
         RAISE_BITS(_status, STATE_INVALID);
         throw ExUninitialized();
      };
   }

// If "EvaluateDmax()" fails, the state will be set to "STATE_INVALID" and background will not be evaluated
   EvaluateDmax(coords);
   if (BITS_RAISED(_status, STATE_INVALID)) throw ExCoordinates();

// Compute physical fields
   EvaluateBackground<Coordinates, Fields, RequestedFields>(coords, fields);
   if (BITS_RAISED(_status, STATE_INVALID)) throw ExFieldError();

   fields.template MakeConsistent<RequestedFields>();

   if (fields.AbsMag() < sp_tiny) RAISE_BITS(_status, STATE_INVALID);
   if (BITS_RAISED(_status, STATE_INVALID)) throw ExFieldError();

// Compute derivatives of fields
   if (RequestedFields::Derived_found())
      EvaluateBackgroundDerivatives<Coordinates, Fields, RequestedFields>(coords, fields);

};

};



