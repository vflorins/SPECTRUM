/*!
\file compiletime_lists.hh
\brief Defines global compile-time constant options that are available, and some associated constant lists.
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_COMPILETIME_LISTS_HH
#define SPECTRUM_COMPILETIME_LISTS_HH

#ifdef USE_GSL
#include <gsl/gsl_const_cgsm.h>
#endif

namespace Spectrum {

//! The largest number of species (distinct particle mass and charge).
#define MAX_PARTICLE_SPECIES 16

//! Elementary charge (esu)
#define SPC_CONST_CGSM_ELECTRON_CHARGE 4.8032044E-10

//! Electron mass (g)
#ifdef USE_GSL
#define SPC_CONST_CGSM_MASS_ELECTRON GSL_CONST_CGSM_MASS_ELECTRON
#else
#define SPC_CONST_CGSM_MASS_ELECTRON 9.11E-26
#endif

//! Proton mass (g)
#ifdef USE_GSL
#define SPC_CONST_CGSM_MASS_PROTON GSL_CONST_CGSM_MASS_PROTON
#else
#define SPC_CONST_CGSM_MASS_PROTON 1.67E-24
#endif

//! Neutron mass (g)
#ifdef USE_GSL
#define SPC_CONST_CGSM_MASS_NEUTRON GSL_CONST_CGSM_MASS_NEUTRON
#else
#define SPC_CONST_CGSM_MASS_NEUTRON SPC_CONST_CGSM_MASS_PROTON
#endif

//! Electron volt (cgs)
#ifdef USE_GSL
#define SPC_CONST_CGSM_ELECTRON_VOLT GSL_CONST_CGSM_ELECTRON_VOLT
#else
#define SPC_CONST_CGSM_ELECTRON_VOLT 1.6E-12
#endif

//! Various particle energy units (cgs)
#define SPC_CONST_CGSM_KILO_ELECTRON_VOLT (1.0E3 * SPC_CONST_CGSM_ELECTRON_VOLT)
#define SPC_CONST_CGSM_MEGA_ELECTRON_VOLT (1.0E6 * SPC_CONST_CGSM_ELECTRON_VOLT)
#define SPC_CONST_CGSM_GIGA_ELECTRON_VOLT (1.0E9 * SPC_CONST_CGSM_ELECTRON_VOLT)


//! Electric charge: typically equal to elementary charge
#define unit_charge_particle SPC_CONST_CGSM_ELECTRON_CHARGE

//! Energy: 1 eV or the proton rest energy are reasonable values
#define unit_energy_particle SPC_CONST_CGSM_ELECTRON_VOLT

//! Momentum: derived, based on the assumption that the particles use the fluid unit for velocity
#define unit_momentum_particle (unit_energy_particle / unit_velocity_fluid)

//! Mass: derived, based on the assumption that the particles use the fluid unit for velocity
#define unit_mass_particle (unit_energy_particle / unit_velocity_fluid / unit_velocity_fluid)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Particle units used in the code - user configurable. Length and time units are the same as for the fluid.
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Electric charge: typically equal to elementary charge
#define unit_charge_particle SPC_CONST_CGSM_ELECTRON_CHARGE

//! Energy: 1 eV or the proton rest energy are reasonable values
#define unit_energy_particle SPC_CONST_CGSM_ELECTRON_VOLT

//! Momentum: derived, based on the assumption that the particles use the fluid unit for velocity
#define unit_momentum_particle (unit_energy_particle / unit_velocity_fluid)

//! Mass: derived, based on the assumption that the particles use the fluid unit for velocity
#define unit_mass_particle (unit_energy_particle / unit_velocity_fluid / unit_velocity_fluid)

//! Ratio of charge/mass has a conversion factor required to calculate particle's cyclotron frequency and radius.
#define charge_mass_particle ((unit_charge_particle * unit_magnetic_fluid / unit_mass_particle / unit_velocity_fluid) / unit_frequency_fluid)

//! Rigidity: derived, 1 V = 1 eV / 1 e
#define unit_rigidity_particle (unit_energy_particle / unit_charge_particle)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Fluid units used in the code - user configurable
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Length: 1 au is a good value for the heliosphere
#ifdef USE_GSL
#define unit_length_fluid GSL_CONST_CGSM_ASTRONOMICAL_UNIT
#else
#define unit_length_fluid 1.5E13
#endif

//! Velocity: 100 km/s is a sensible value; do not use the speed of light because some of the tests will return infinity.
#define unit_velocity_fluid 1.0E7

//! Time: derived, should be hours or days
#define unit_time_fluid (unit_length_fluid / unit_velocity_fluid)

//! Frequency: derived
#define unit_frequency_fluid (1.0 / unit_time_fluid)

//! Number density: a good number is 1 per cm^3. We can choose a number density unit independently of the length unit because certain conservation laws (namely MHD) are invariant under the transformation rho->a*rho, p->a*p, B^2->a*B^2. However, this is not the case for other physical laws such as those governing collisions, where corrections must are applied.
#define unit_number_density_fluid 1.0

//! Density: derived, should be 1 _particle_ mass unit times 1 number density unit
#define unit_density_fluid (unit_mass_particle * unit_number_density_fluid)

//! Momentum: derived
#define unit_momentum_fluid (unit_density_fluid * unit_velocity_fluid)

//! Magnetic field: derived, should be uG-mG
#define unit_magnetic_fluid (unit_velocity_fluid * sqrt(unit_density_fluid))

//! Electric field: derived, should be the same as magnetic field because of cgs-Gaussian unit system
#define unit_electric_fluid unit_magnetic_fluid

//! Pressure: derived, should be ~10^-10
#define unit_pressure_fluid (unit_magnetic_fluid * unit_magnetic_fluid)

//! Diffusion coefficient: derived, should be 10^20-10^23 (good numbers for GCRs)
#define unit_diffusion_fluid (unit_velocity_fluid * unit_length_fluid)

//! Temperature: define unit for a typical energy of 1 eV (10^4 K)
#ifdef USE_GSL
#define unit_temperature_fluid (GSL_CONST_CGSM_ELECTRON_VOLT / GSL_CONST_CGSM_BOLTZMANN)
#else
#define unit_temperature_fluid 1.6E-12 / 1.36E-16
#endif

//! Speed of light
#ifdef USE_GSL
#define c_code (GSL_CONST_CGSM_SPEED_OF_LIGHT / unit_velocity_fluid)
#else
#define c_code (3.0E10 / unit_velocity_fluid)
#endif

//! Speed of light squared
#define c2_code (c_code * c_code)

//! Boltzmann constant. Note that for fluids the ideal gas EOS is p=n*kb*T.
#ifdef USE_GSL
#define kb_code (GSL_CONST_CGSM_BOLTZMANN / unit_pressure_fluid * unit_temperature_fluid * unit_number_density_fluid)
#else
#define kb_code (1.38E-16 / unit_pressure_fluid * unit_temperature_fluid * unit_number_density_fluid)
#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// General options
//----------------------------------------------------------------------------------------------------------------------------------------------------


//! Build mode (debug build or production build)
enum class BuildMode {
   debug,
   release
};

//! Specie (alternative impl)
enum class Specie {
   electron_core,
   electron_halo,
   electron_beam,
   //
   proton_core,
   proton_halo,
   proton_beam,
   //
   proton_pickup,
   alpha_core,
   alpha_halo,
   //
   heliumII_core,
   heliumII_pickup,
   hydrogenII_core,
   //
   hydrogenI_core,
   hydrogenI_halo,
   hydrogenI_beam,
   heliumI_core
};

//! Names of fluid species
constexpr std::array<std::string_view, MAX_PARTICLE_SPECIES> SpeciesNames
      = {
            "Electron (core)", "Electron (halo)", "Electron (beam)",
            "Proton (core)", "Proton (halo)", "Proton (beam)",
            "Proton (pickup)", "Alpha (core)", "Alpha (halo)",
            "Helium II (core)", "Helium II (pickup)", "Hydrogen II (core)",
            "Hydrogen I (core)", "Hydrogen I (halo)", "Hydrogen I (beam)","Helium I (core)"};


//! Masses of particles comprising the fluids
constexpr std::array<double, MAX_PARTICLE_SPECIES> SpeciesMasses = {   SPC_CONST_CGSM_MASS_ELECTRON  / unit_mass_particle,
                                                                       SPC_CONST_CGSM_MASS_ELECTRON  / unit_mass_particle,
                                                                       SPC_CONST_CGSM_MASS_ELECTRON  / unit_mass_particle,
                                                                       SPC_CONST_CGSM_MASS_PROTON    / unit_mass_particle,
                                                                       SPC_CONST_CGSM_MASS_PROTON    / unit_mass_particle,
                                                                       SPC_CONST_CGSM_MASS_PROTON    / unit_mass_particle,
                                                                       SPC_CONST_CGSM_MASS_PROTON    / unit_mass_particle,
                                                                       2.0 * (SPC_CONST_CGSM_MASS_PROTON   + SPC_CONST_CGSM_MASS_NEUTRON)  / unit_mass_particle,
                                                                       2.0 * (SPC_CONST_CGSM_MASS_PROTON   + SPC_CONST_CGSM_MASS_NEUTRON)  / unit_mass_particle,
                                                                       (2.0 * (SPC_CONST_CGSM_MASS_PROTON + SPC_CONST_CGSM_MASS_NEUTRON) + SPC_CONST_CGSM_MASS_ELECTRON) / unit_mass_particle,
                                                                       (2.0 * (SPC_CONST_CGSM_MASS_PROTON + SPC_CONST_CGSM_MASS_NEUTRON) + SPC_CONST_CGSM_MASS_ELECTRON) / unit_mass_particle,
                                                                       (SPC_CONST_CGSM_MASS_PROTON   + SPC_CONST_CGSM_MASS_ELECTRON) / unit_mass_particle,
                                                                       (SPC_CONST_CGSM_MASS_PROTON   + SPC_CONST_CGSM_MASS_ELECTRON) / unit_mass_particle,
                                                                       (SPC_CONST_CGSM_MASS_PROTON   + SPC_CONST_CGSM_MASS_ELECTRON) / unit_mass_particle,
                                                                       (SPC_CONST_CGSM_MASS_PROTON   + SPC_CONST_CGSM_MASS_ELECTRON) / unit_mass_particle,
                                                                       2.0 * (SPC_CONST_CGSM_MASS_PROTON + SPC_CONST_CGSM_MASS_NEUTRON  + SPC_CONST_CGSM_MASS_ELECTRON) / unit_mass_particle};

//! Charges of particles comprising the fluids
constexpr std::array<double, MAX_PARTICLE_SPECIES> SpeciesCharges = {-SPC_CONST_CGSM_ELECTRON_CHARGE / unit_charge_particle,
                                                                     -SPC_CONST_CGSM_ELECTRON_CHARGE / unit_charge_particle,
                                                                     -SPC_CONST_CGSM_ELECTRON_CHARGE / unit_charge_particle,
                                                                     SPC_CONST_CGSM_ELECTRON_CHARGE / unit_charge_particle,
                                                                     SPC_CONST_CGSM_ELECTRON_CHARGE / unit_charge_particle,
                                                                     SPC_CONST_CGSM_ELECTRON_CHARGE / unit_charge_particle,
                                                                     SPC_CONST_CGSM_ELECTRON_CHARGE / unit_charge_particle,
                                                                     2.0 * SPC_CONST_CGSM_ELECTRON_CHARGE / unit_charge_particle,
                                                                     2.0 * SPC_CONST_CGSM_ELECTRON_CHARGE / unit_charge_particle,
                                                                     SPC_CONST_CGSM_ELECTRON_CHARGE / unit_charge_particle,
                                                                     SPC_CONST_CGSM_ELECTRON_CHARGE / unit_charge_particle,
                                                                     0.0, 0.0, 0.0, 0.0, 0.0};

//! Polytropic indices of fluids
constexpr std::array<double, MAX_PARTICLE_SPECIES> SpeciesPolytropicIndices
      = {5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0,
         5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0};


//! Derivative method (generally either analytic or numeric)
enum class DerivativeMethod {
   analytic,
   numeric
};

//! Direction of trajectory integration (time flow direction)
enum class TimeFlow {
   forward,
   backward
};


//----------------------------------------------------------------------------------------------------------------------------------------------------
// Background options
//----------------------------------------------------------------------------------------------------------------------------------------------------


namespace SolarWindOptions {

//! Heliospheric current sheet (0: disabled, 1: flat, 2: wavy (Jokipii-Thomas 1981) and static, 3: wavy and time-dependent).
enum class CurrentSheet {
   disabled,
   flat,
   wavy_static,
   wavy_time_dependent
};

//! Magnetic topology region (0: nowhere, 1: same as HCS)
enum class SectoredRegion {
   nowhere,
   HCS
};

//! Correction to Parker Spiral, mainly for polar regions (0: none, 1: Smith-Bieber 1991, 2: Zurbuchen et al. 1997, 3: Schwadron-McComas 2003)
enum class PolarCorrection {
   none,
   Smith_Bieber,
   Zurbuchen_etal,
   Schwadron_McComas
};

//! Latitudinal profile for bulk speed (0: constant, 1: linear step, 2: smooth step)
enum class SpeedLatitudeProfile {
   constant,
   linear_step,
   smooth_step
};

//! Integer exponent of decrease of solar wind speed beyond the termination shock
enum class TermShockSpeedExponent {
   zero,
   one,
   square,
   cube
};

}

namespace VLISMBochumOptions {

//! What function to use within 'get_ampfactor' (0 = none, 1 = zero, 2 = constant, 3 = scaled)
enum class ModType {
   none,
   zero,
   constant,
   scaled
};

//! Whether to scale relative to s=0 (0) or s=+inf (1)
enum class ModRPos {
   scale_rel_zero,
   scale_rel_inf
};

}

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Trajectory options
//----------------------------------------------------------------------------------------------------------------------------------------------------



namespace TrajectoryOptions {

//! How many time steps to allow before recording a mirror event
//const int mirror_threshold;


enum class SafetyLevel {
   low,
   medium,
   high
};

// Switch controlling how to calculate mu. "0" means computing it at the end of the step from magnetic moment conservation. "1" means advancing it in time according to the scheme (does not guarantee conservation of MM, but can be used with non-adiabatic terms).
enum class PPerpMethod {
   mag_moment_conservation,
   scheme
};

//! Flag to use gradient and curvature drifts in drift velocity calculation
enum class UseBDrifts {
   none,
   gradient_curvature
};

//! Which stochastic method to use for PA scattering, 0 = Euler, 1 = Milstein, 2 = RK2
enum class StochasticMethod {
   Euler,
   Milstein,
   RK2
};

//! Whether to use constant dmumax or constant dthetamax, 0 = constant dthetamax, 1 = constant dmumax
enum class ConstDmumax {
   constant_dtheta_max,
   constant_dmu_max
};

//! Which method of computation to use for divK: 0 = using direct central FD, 1 = using _spdata.grad quantities
enum class DivkMethod {
   direct,
   gradients
};



}




};

#endif
