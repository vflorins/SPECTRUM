/*!
\file specie.hh
\brief Declares some elementary particle properties
\author Lucius Schoenbaum
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SPECIE_HH
#define SPECTRUM_SPECIE_HH

#include <array>
#include <string_view>

#include "config.h"
#include "common/gpu_config.hh"
#include "common/physical_constants.hh"

namespace Spectrum {

//! The largest number of species (distinct particle mass and charge).
#define MAX_PARTICLE_SPECIES 17

namespace Particle {

// TODO: The following four constants should be moved to a user-configurable location

//! Length unit: primary (choose based on on Larmor radius, such as 1000 km)
SPECTRUM_CONSTEXPR double unit_length = 1.0E10;
//SPECTRUM_CONSTEXPR double unit_length = 1.5E13;

//! Velocity unit: primary (~c)
SPECTRUM_CONSTEXPR double unit_velocity = 0.1 * SPC_CONST_CGSM_SPEED_OF_LIGHT;
//SPECTRUM_CONSTEXPR double unit_velocity = 1.0E7;

//! Mass unit: primary (m_p, m_e, etc.)
SPECTRUM_CONSTEXPR double unit_mass = SPC_CONST_CGSM_MASS_PROTON;

//! Charge unit: arbitrary (set for convenience to e and not usable without conversion)
SPECTRUM_CONSTEXPR double unit_charge = SPC_CONST_CGSM_ELECTRON_CHARGE;

//! Temperature: primary (use for PUIs and other low-energy populations)
SPECTRUM_CONSTEXPR double unit_temperature = SPC_CONST_CGSM_KILO_ELECTRON_VOLT / SPC_CONST_CGSM_BOLTZMANN;

//! Number density: arbitrary (not usable without conversion)
SPECTRUM_CONSTEXPR double unit_number_density = 1.0;

//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Momentum unit: derived from mass and velocity
SPECTRUM_CONSTEXPR double unit_momentum = unit_mass * unit_velocity;

//! Energy unit: derived from mass and velocity
SPECTRUM_CONSTEXPR double unit_energy = unit_mass * unit_velocity * unit_velocity;

//! Time unit: derived from velocity and length (should be ~ gyro-period)
SPECTRUM_CONSTEXPR double unit_time = unit_length / unit_velocity;

//! Field unit: derived from momentum, time, and charge (should be ~1 nT)
SPECTRUM_CONSTEXPR double unit_magnetic = unit_momentum / (unit_charge * unit_time);

//! Pressure and energy density: derived from magnetic field
SPECTRUM_CONSTEXPR double unit_pressure = unit_magnetic * unit_magnetic;

//! Force unit: derived from momentum and time
SPECTRUM_CONSTEXPR double unit_force = unit_momentum / unit_time;

//! Rigidity unit: derived from energy and charge
SPECTRUM_CONSTEXPR double unit_rigidity = unit_energy / unit_charge;

//! Diffusion coefficient: derived, should be 10^20-10^23 (good numbers for GCRs)
SPECTRUM_CONSTEXPR double unit_diffusion = unit_velocity * unit_length;

//! Magnetic moment: derived from charge and length
SPECTRUM_CONSTEXPR double unit_magnetic_moment = unit_charge * unit_length;

//! Conversion factor for charge
//SPECTRUM_CONSTEXPR double charge_conv = unit_charge * unit_magnetic * unit_time / unit_momentum;

//! Speed of light in particle units
SPECTRUM_CONSTEXPR double c_code = SPC_CONST_CGSM_SPEED_OF_LIGHT / unit_velocity;

//! Boltzmann constant
SPECTRUM_CONSTEXPR double kb_code = SPC_CONST_CGSM_BOLTZMANN * unit_temperature / unit_energy;

};

//! Specie identifiers for lookup
enum class SpecieId
{
//! Electrons - core
   electron_core,

//! Electrons - energetic
   electron_halo,

//! Electrons - energetic/anisotropic (e.g., strahl)
   electron_beam,

//! Muons - core
   muon_core,

//! Protons - core
   proton_core,

//! Protons - energetic (e.g., suprathermal)
   proton_halo,

//! Protons - energetic/anisotropic
   proton_beam,

//! Protons - pickup (as proton_halo, but different origin)
   proton_pickup,

//! Alpha particles - core
   alpha_core,

//! Alpha particles - energetic
   alpha_halo,

//! Singly ionized helium
   heliumII_core,

//! Singly ionized helium - pickup
   heliumII_pickup,

//! Proton-electron plasma - core
   hydrogenII_core,

//! Hydrogen atoms - core
   hydrogenI_core,

//! Hydrogen atoms - energetic (e.g., from heliosheath)
   hydrogenI_halo,

//! Hydrogen atoms - energetic (e.g., from solar wind)
   hydrogenI_beam,

//! Helium atoms - core
   heliumI_core
};

//! Default specie
SPECTRUM_CONSTEXPR SpecieId default_specie = SpecieId::proton_core;

namespace Specie_impl {

#ifndef __CUDACC__

//! Names of fluid species
constexpr std::array<std::string_view, MAX_PARTICLE_SPECIES> SpeciesNames
= {"Electron (core)", "Electron (halo)", "Electron (beam)",
   "Muon (core)",
   "Proton (core)", "Proton (halo)", "Proton (beam)", "Proton (pickup)",
   "Alpha (core)", "Alpha (halo)",
   "Helium II (core)", "Helium II (pickup)",
   "Hydrogen II (core)",
   "Hydrogen I (core)", "Hydrogen I (halo)", "Hydrogen I (beam)",
   "Helium I (core)"};

#endif

//! Masses of particles comprising the fluids
SPECTRUM_CONSTEXPR double SpeciesMasses[MAX_PARTICLE_SPECIES] =       {SPC_CONST_CGSM_MASS_ELECTRON  / Particle::unit_mass,
                                                                       SPC_CONST_CGSM_MASS_ELECTRON  / Particle::unit_mass,
                                                                       SPC_CONST_CGSM_MASS_ELECTRON  / Particle::unit_mass,
                                                                       SPC_CONST_CGSM_MASS_MUON      / Particle::unit_mass,
                                                                       SPC_CONST_CGSM_MASS_PROTON    / Particle::unit_mass,
                                                                       SPC_CONST_CGSM_MASS_PROTON    / Particle::unit_mass,
                                                                       SPC_CONST_CGSM_MASS_PROTON    / Particle::unit_mass,
                                                                       SPC_CONST_CGSM_MASS_PROTON    / Particle::unit_mass,
                                 2.0 * (SPC_CONST_CGSM_MASS_PROTON   + SPC_CONST_CGSM_MASS_NEUTRON)  / Particle::unit_mass,
                                 2.0 * (SPC_CONST_CGSM_MASS_PROTON   + SPC_CONST_CGSM_MASS_NEUTRON)  / Particle::unit_mass,
   (2.0 * (SPC_CONST_CGSM_MASS_PROTON + SPC_CONST_CGSM_MASS_NEUTRON) + SPC_CONST_CGSM_MASS_ELECTRON) / Particle::unit_mass,
   (2.0 * (SPC_CONST_CGSM_MASS_PROTON + SPC_CONST_CGSM_MASS_NEUTRON) + SPC_CONST_CGSM_MASS_ELECTRON) / Particle::unit_mass,
                                       (SPC_CONST_CGSM_MASS_PROTON   + SPC_CONST_CGSM_MASS_ELECTRON) / Particle::unit_mass,
                                       (SPC_CONST_CGSM_MASS_PROTON   + SPC_CONST_CGSM_MASS_ELECTRON) / Particle::unit_mass,
                                       (SPC_CONST_CGSM_MASS_PROTON   + SPC_CONST_CGSM_MASS_ELECTRON) / Particle::unit_mass,
                                       (SPC_CONST_CGSM_MASS_PROTON   + SPC_CONST_CGSM_MASS_ELECTRON) / Particle::unit_mass,
    2.0 * (SPC_CONST_CGSM_MASS_PROTON + SPC_CONST_CGSM_MASS_NEUTRON  + SPC_CONST_CGSM_MASS_ELECTRON) / Particle::unit_mass};

//! Charges of particles comprising the fluids
SPECTRUM_CONSTEXPR double SpeciesCharges[MAX_PARTICLE_SPECIES] = {-SPC_CONST_CGSM_ELECTRON_CHARGE / Particle::unit_charge,
                                                                  -SPC_CONST_CGSM_ELECTRON_CHARGE / Particle::unit_charge,
                                                                  -SPC_CONST_CGSM_ELECTRON_CHARGE / Particle::unit_charge,
                                                                  -SPC_CONST_CGSM_ELECTRON_CHARGE / Particle::unit_charge,
                                                                   SPC_CONST_CGSM_ELECTRON_CHARGE / Particle::unit_charge,
                                                                   SPC_CONST_CGSM_ELECTRON_CHARGE / Particle::unit_charge,
                                                                   SPC_CONST_CGSM_ELECTRON_CHARGE / Particle::unit_charge,
                                                                   SPC_CONST_CGSM_ELECTRON_CHARGE / Particle::unit_charge,
                                                             2.0 * SPC_CONST_CGSM_ELECTRON_CHARGE / Particle::unit_charge,
                                                             2.0 * SPC_CONST_CGSM_ELECTRON_CHARGE / Particle::unit_charge,
                                                                   SPC_CONST_CGSM_ELECTRON_CHARGE / Particle::unit_charge,
                                                                   SPC_CONST_CGSM_ELECTRON_CHARGE / Particle::unit_charge,
                                                                   0.0, 0.0, 0.0, 0.0, 0.0};

//! Polytropic indices (degrees of freedom)
SPECTRUM_CONSTEXPR double SpeciesPolytropicIndices[MAX_PARTICLE_SPECIES] = {5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0,
                                                                            5.0 / 3.0,
                                                                            5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0,
                                                                            5.0 / 3.0, 5.0 / 3.0,
                                                                            5.0 / 3.0, 5.0 / 3.0,
                                                                            5.0 / 3.0,
                                                                            5.0 / 3.0, 5.0 / 3.0, 5.0 / 3.0,
                                                                            5.0 / 3.0};

};

//! All specie properties
template <SpecieId specieid>
struct Specie
{
//! Index of the specie from "SpecieId"
   static constexpr SpecieId id = specieid;

#ifndef __CUDACC__
//! Readable name of the specie
   static constexpr std::string_view name = Specie_impl::SpeciesNames[static_cast<size_t>(specieid)];
#endif

//! Mass of the particle
   static constexpr double mass = Specie_impl::SpeciesMasses[static_cast<size_t>(specieid)];

//! Signed charge of the particle
   static constexpr double charge = Specie_impl::SpeciesCharges[static_cast<size_t>(specieid)];

//! Degrees of freedom
   static constexpr double pt_idx = Specie_impl::SpeciesPolytropicIndices[static_cast<size_t>(specieid)];
};

};

#endif
