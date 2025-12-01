/*!
\file utils_numerical_derivatives.cc
\brief Implements numerical evaluation of gradients based on the plasma background
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "utils_numerical_derivatives.hh"
#include "common/status.hh"

namespace Spectrum {

using namespace BackgroundOptions;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// NumericalDerivatives methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/05/2025
\param[in] int index for which derivative to take (0 = x, 1 = y, 2 = z, else = t)
\param[in] coords coordinates requesting the computation (copy consumed by computation)
\param[out] fields requested fields (writeout argument)
\param[in] scale_factor a factor establishing the scale at which derivatives are significant.
\note This is a common routine that the derived classes should not change.
*/
template <typename Background>
template <typename Coordinates, typename Fields, typename RequestedFields>
void NumericalDerivatives<Background, true>::DirectionalDerivative(const int xyz, Coordinates coords, Fields& fields, double scale_factor, Background& background)
{
// temporaries for forward- and backward-stepped evaluation
   Fields fields_forw, fields_back;
   status_t status;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Spatial derivatives
   if ((0 <= xyz) && (xyz <= 2)) {

      _ddata._dr_forw_fail[xyz] = false;
      _ddata._dr_back_fail[xyz] = false;

// Increment forward, then back.
// In either case, if the attempted state is invalid,
// revert to first order differences, otherwise use second order differences.

// Forward increment
      _ddata._dr[xyz] = scale_factor;
      coords.Pos('w') += _ddata._dr[xyz] * fa_basis.row[xyz];
      status = background.template EvaluateBackground<Coordinates, Fields, RequestedFields>(coords, fields_forw);

      if (BITS_RAISED(status, STATE_INVALID)) {
         fields_forw = fields;
         _ddata._dr_forw_fail[xyz] = true;
      };

// Backward increment
      _ddata._dr[xyz] *= 2.0;
      coords.Pos('w') -= _ddata._dr[xyz] * fa_basis.row[xyz];
      status = background.template EvaluateBackground<Coordinates, Fields, RequestedFields>(coords, fields_back);

      if (BITS_RAISED(status, STATE_INVALID)) {
         fields_back = fields;
         _ddata._dr_back_fail[xyz] = true;
      };

// If at least one increment failed, half the increment. If both increments failed, throw an error.
      if (_ddata._dr_forw_fail[xyz] || _ddata._dr_back_fail[xyz]) _ddata._dr[xyz] *= 0.5;
      if (_ddata._dr_forw_fail[xyz] && _ddata._dr_back_fail[xyz]) throw ExFieldError();

// Compute spatial derivatives. This calculation gives gradV[i][j] = dV_j / ds^i which is the transpose of the Jacobian.
      if constexpr (RequestedFields::DelFluv_found())
         fields.DelFluv('w')[xyz] = (fields_forw.Fluv() - fields_back.Fluv()) / _ddata._dr[xyz];
      if constexpr (RequestedFields::DelMag_found())
         fields.DelMag('w')[xyz] = (fields_forw.Mag() - fields_back.Mag()) / _ddata._dr[xyz];
      if constexpr (RequestedFields::DelAbsMag_found())
         fields.DelAbsMag('w')[xyz] = (fields_forw.AbsMag() - fields_back.AbsMag()) / _ddata._dr[xyz];
      if constexpr (RequestedFields::DelElc_found())
         fields.DelElc('w')[xyz] = (fields_forw.Elc() - fields_back.Elc()) / _ddata._dr[xyz];
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
      _ddata._dt = 1.0 / scale_factor;
      coords.Time('w') += _ddata._dt;
      status = background.template EvaluateBackground<Coordinates, Fields, RequestedFields>(coords, fields_forw);

      if (BITS_RAISED(status, STATE_INVALID)) {
         fields_forw = fields;
         _ddata._dt_forw_fail = true;
      };

// Backward increment
      _ddata._dt *= 2.0;
      coords.Time('w') -= _ddata._dt;
      status = background.template EvaluateBackground<Coordinates, Fields, RequestedFields>(coords, fields_back);

      if (BITS_LOWERED(status, STATE_INVALID)) {
         fields_back = fields;
         _ddata._dt_back_fail = true;
      };

// If at least one increment failed, half the increment. If both increments failed, throw an error.
      if (_ddata._dt_forw_fail || _ddata._dt_back_fail) _ddata._dt *= 0.5;
      if (_ddata._dt_forw_fail && _ddata._dt_back_fail) throw ExFieldError();

// Compute time derivatives
      if constexpr (RequestedFields::DotFluv_found())
         fields.DotFluv('w') = (fields_forw.Fluv() - fields_back.Fluv()) / _ddata._dt;
      if constexpr (RequestedFields::DotMag_found())
         fields.DotMag('w') = (fields_forw.Mag() - fields_back.Mag()) / _ddata._dt;
      if constexpr (RequestedFields::DotAbsMag_found())
         fields.DotAbsMag('w') = (fields_forw.AbsMag() - fields_back.AbsMag()) / _ddata._dt;
      if constexpr (RequestedFields::DotElc_found())
         fields.DotElc('w') = (fields_forw.Elc() - fields_back.Elc()) / _ddata._dt;
   };
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 11/28/2025
*/
template <typename Background>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t NumericalDerivatives<Background, true>::EvaluateBackgroundDerivatives(Coordinates& coords, Fields& fields, Background& background)
{
   double AbsMom = coords.AbsMom();
   double AbsMag = fields.AbsMag();
   GeoVector HatMag = fields.HatMag();
   constexpr bool gradients = RequestedFields::DelFluv_found() || RequestedFields::DelMag_found() || RequestedFields::DelAbsMag_found() || RequestedFields::DelElc_found();
   constexpr bool time_derivatives =  RequestedFields::DotFluv_found() || RequestedFields::DotMag_found() || RequestedFields::DotAbsMag_found() || RequestedFields::DotElc_found();

// This method is only valid in the presence of a magnetic field.
// Normally magnetic field magnitude+direction is tracked but we can make do with magnetic field only.

   if constexpr (gradients) {

// Derivatives are only needed for trajectory types whose transport assumes the background changes on scales larger than the gyro-radius.
      auto r_g = fmin(LarmorRadius<Config::specie>(AbsMom, AbsMag), _ddata.dmax);

// Get field aligned basis in (transpose) of rotation matrix
      fa_basis.row[2] = HatMag;
      fa_basis.row[0] = GetSecondUnitVec(HatMag);
      fa_basis.row[1] = fa_basis.row[2] ^ fa_basis.row[0];

// Compute derivatives in field-aligned basis
      for (auto xyz = 0; xyz < 3; xyz++)
         DirectionalDerivative<Background, Coordinates, Fields, RequestedFields>(xyz, coords, fields, r_g, background);

// Transform basis back to global cartesian frame
      rot_mat.Transpose(fa_basis);
      if constexpr (RequestedFields::DelFluv_found())
         fields.DelFluv('w') = rot_mat * fields.DelFluv();
      if constexpr (RequestedFields::DelMag_found())
         fields.DelMag('w') = rot_mat * fields.DelMag();
      if constexpr (RequestedFields::DelAbsMag_found())
         fields.DelAbsMag('w') = rot_mat * fields.DelAbsMag();
      if constexpr (RequestedFields::DelElc_found())
         fields.DelElc('w') = rot_mat * fields.DelElc();

      if constexpr (Config::num_numeric_grad_evals > 1) {
         constexpr int num_evals = Config::num_numeric_grad_evals;
         // copy fields into temporary
         Fields fields_tmp = fields;

// Repeat for each additional rotation
         for (auto n_rot = 1; n_rot < num_evals; ++n_rot) {
            fa_basis[0].Rotate(fa_basis.row[2], sin_lra, cos_lra);
            fa_basis[1].Rotate(fa_basis.row[2], sin_lra, cos_lra);

            for (auto xyz = 0; xyz < 3; xyz++)
               DirectionalDerivative<Background, Coordinates, Fields, RequestedFields>(xyz, coords, fields_tmp, r_g, background);

            rot_mat.Transpose(fa_basis);
            if constexpr (RequestedFields::DelFluv_found())
               fields.DelFluv('w') += rot_mat * fields_tmp.DelFluv();
            if constexpr (RequestedFields::DelMag_found())
               fields.DelMag('w') += rot_mat * fields_tmp.DelMag();
            if constexpr (RequestedFields::DelAbsMag_found())
               fields.DelAbsMag('w') += rot_mat * fields_tmp.DelAbsMag();
            if constexpr (RequestedFields::DelElc_found())
               fields.DelElc('w') += rot_mat * fields_tmp.DelElc();

         };

// Average results
         if constexpr (RequestedFields::DelFluv_found())
            fields.DelFluv('w') /= num_evals;
         if constexpr (RequestedFields::DelMag_found())
            fields.DelMag('w') /= num_evals;
         if constexpr (RequestedFields::DelAbsMag_found())
            fields.DelAbsMag('w') /= num_evals;
         if constexpr (RequestedFields::DelElc_found())
            fields.DelElc('w') /= num_evals;
      };
   };

// Time derivatives.
   if constexpr (time_derivatives) {
// Derivatives are only needed for trajectory types whose transport assumes the background changes on scales longer than the gyro-frequency.
      auto vel = Vel<Config::specie>(AbsMom);
      auto w_g = fmin(CyclotronFrequency<Config::specie>(vel, AbsMag), vel / _ddata.dmax);
      DirectionalDerivative<Background, Coordinates, Fields, RequestedFields>(3, coords, fields, 1.0/w_g, background);
   };

   return 0;
};



};



