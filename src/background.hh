//
// Created by Lucius Schoenbaum on 9/8/25.
//

#ifndef SPECTRUM_BACKGROUND_HH
#define SPECTRUM_BACKGROUND_HH

#include "common/compiletime_lists.hh"

#include "background_cylindrical_obstacle.hh"
#include "background_dipole.hh"
#include "background_discontinuity.hh"
#include "background_magnetized_cylinder.hh"
#include "background_magnetized_sphere.hh"
#include "background_server.hh"
#include "background_server_batl.hh"
#include "background_server_cartesian.hh"
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
// todo fix - ServerFront
BackgroundServer<HConfig, nullptr_t>,
BackgroundServerBATL<HConfig, nullptr_t>,
BackgroundServerCartesian<HConfig, nullptr_t>,
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
using Background = FieldOps::Nth<BackgroundList<HConfig>, reinterpret_cast<int>(HConfig::BackgroundConfig::backgroundid)>;


}



#endif
