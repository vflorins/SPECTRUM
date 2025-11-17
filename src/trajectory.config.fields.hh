/*!
\file trajectory.config.fields.hh
\brief Default field and coordinate types for a SPECTRUM trajectory class
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_TRAJECTORY_CONFIG_FIELDS_HH
#define SPECTRUM_TRAJECTORY_CONFIG_FIELDS_HH

#include "common/compiletime_lists.hh"
#include "common/fields.hh"


/*
 * These types are only default choices.
 * A custom/unique choice can be defined via the TrajectoryConfig type.
 *
 */


namespace Spectrum {

template<Config::Trajectory trajectoryid, SpecieId specieid>
class TrajectoryCoordinates;

template<Config::Trajectory trajectoryid, SpecieId specieid>
class TrajectoryFields;


template<SpecieId specieid>
class TrajectoryCoordinates<Config::Trajectory::Fieldline, specieid> {
   using type = Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::anisotropic>, Pos_t, Time_t, Mom_t, Vel_t>;
};

template<SpecieId specieid>
class TrajectoryFields<Config::Trajectory::Fieldline, specieid> {
   /* Unused */
   using type = Fields<FConfig<specieid>>;
};


template<SpecieId specieid>
class TrajectoryCoordinates<Config::Trajectory::Lorentz, specieid> {
   using type = Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::cartesian>, Pos_t, Time_t, Mom_t, Vel_t>;
};

template<SpecieId specieid>
class TrajectoryFields<Config::Trajectory::Lorentz, specieid> {
   using type = Fields<FConfig<specieid>, Mag_t, Elc_t>;
};


template<SpecieId specieid>
class TrajectoryCoordinates<Config::Trajectory::Focused, specieid> {
   using type = Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Mom_t, Vel_t>;
};

template<SpecieId specieid>
class TrajectoryFields<Config::Trajectory::Focused, specieid> {
   using type = Fields<FConfig<specieid>, Fluv_t, Mag_t, AbsMag_t, HatMag_t, DelFluv_t, DelAbsMag_t, DelMag_t, DotFluv_t>;
};


template<SpecieId specieid>
class TrajectoryCoordinates<Config::Trajectory::Guiding, specieid> {
   using type = Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::anisotropic>, Pos_t, Time_t, Mom_t, Vel_t>;
};

template<SpecieId specieid>
class TrajectoryFields<Config::Trajectory::Guiding, specieid> {
   using type = Fields<FConfig<specieid>, Fluv_t, Mag_t, Elc_t, AbsMag_t, HatMag_t, DelMag_t, DelAbsMag_t, DotAbsMag_t>;
};


template<SpecieId specieid>
class TrajectoryCoordinates<Config::Trajectory::Parker, specieid> {
   using type = Fields<FConfig<specieid, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Mom_t, Vel_t>;
};

template<SpecieId specieid>
class TrajectoryFields<Config::Trajectory::Parker, specieid> {
   using type = Fields<FConfig<specieid>, Fluv_t, Mag_t, AbsMag_t, HatMag_t, DelMag_t, DelAbsMag_t>;
};


}


#endif
