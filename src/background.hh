//
// Created by Lucius Schoenbaum on 9/8/25.
//

#ifndef SPECTRUM_BACKGROUND_HH
#define SPECTRUM_BACKGROUND_HH

#include "common/compiletime_lists.hh"
#include "common/fields.hh"

#include "background_cylindrical_obstacle.hh"
#include "background_dipole.hh"
#include "background_discontinuity.hh"
#include "background_magnetized_cylinder.hh"
#include "background_magnetized_sphere.hh"
#include "background_data_batl.hh"
#include "background_data_cartesian.hh"
#include "background_shock.hh"
#include "background_smooth_discontinuity.hh"
#include "background_smooth_shock.hh"
#include "background_solarwind.hh"
#include "background_solarwind_termshock.hh"
#include "background_spherical_obstacle.hh"
#include "background_uniform.hh"
#include "background_vlism_bochum.hh"
#include "background_waves.hh"

namespace Spectrum {


template<typename HConfig>
using BackgroundList = Fields<
FConfig<>,
BackgroundCylindricalObstacle<HConfig>,
BackgroundDipole<HConfig>,
BackgroundDiscontinuity<HConfig>,
BackgroundMagnetizedCylinder<HConfig>,
BackgroundMagnetizedSphere<HConfig>,
BackgroundShock<HConfig>,
BackgroundSmoothDiscontinuity<HConfig>,
BackgroundSmoothShock<HConfig>,
BackgroundSolarWind<HConfig>,
BackgroundSolarWindTermShock<HConfig>,
BackgroundSphericalObstacle<HConfig>,
BackgroundUniform<HConfig>,
BackgroundVLISMBochum<HConfig>,
BackgroundWaves<HConfig>
>;


template<typename HConfig>
using Background = FieldOps::Nth<BackgroundList<HConfig>, static_cast<int>(HConfig::background)>;


}



/*
 *
 * Documentation (Work in Progress)
 *
 * The Background methods Evaluate, EvaluateDerivatives, and EvaluateDmax
 * are valid for any Coordinate type that includes position and time (Pos_t, Time_t) and
 * magnitude of momentum (AbsMom_t or at least Mom_t).
 * The fields type (Fields) must always contain magnetic field (Mag_t).
 *
 */

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


#endif