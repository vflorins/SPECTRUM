/*!
\file physics.hh
\brief Declares physical constants and conversion routines for various physical quantities
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_PHYSICS_HH
#define SPECTRUM_PHYSICS_HH

#include <array>

#include <config.h>

#ifdef USE_GSL
#include <gsl/gsl_const_cgsm.h>
#endif

#include <common/vectors.hh>

namespace Spectrum {

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

//! The largest nmumber of species (distinct particle mass and charge).
#define MAX_PARTICLE_SPECIES 16

//! Electrons - core
#define SPECIES_ELECTRON_CORE 0

//! Electrons - energetic
#define SPECIES_ELECTRON_HALO 1

//! Electrons - energetic/anisotropic (e.g., strahl)
#define SPECIES_ELECTRON_BEAM 2

//! Protons - core
#define SPECIES_PROTON_CORE 3

//! Protons - energetic (e.g., suprathermal)
#define SPECIES_PROTON_HALO 4

//! Protons - energetic/anisotropic
#define SPECIES_PROTON_BEAM 5

//! Protons - pickup (as SPECIES_PROTON_HALO, but different origin)
#define SPECIES_PROTON_PICKUP 6

//! Alpha particles - core
#define SPECIES_ALPHA_CORE 7

//! Alpha particles - energetic
#define SPECIES_ALPHA_HALO 8

//! Singly ionized helium
#define SPECIES_HELIUM_SINGLE_CORE 9

//! Singly ionized helium - pickup
#define SPECIES_HELIUM_SINGLE_PICKUP 10

//! Proton-electron plasma - core
#define SPECIES_HYDROGEN_PLASMA_CORE 11

//! Hydrogen atoms - core
#define SPECIES_HYDROGEN_CORE 12

//! Hydrogen atoms - energetic (e.g., from heliosheath)
#define SPECIES_HYDROGEN_HALO 13

//! Hydrogen atoms - energetic (e.g., from solar wind)
#define SPECIES_HYDROGEN_BEAM 14

//! Helium atoms - core
#define SPECIES_HELIUM_CORE 15

//! Names of fluids
constexpr std::array<std::string_view, MAX_PARTICLE_SPECIES> SpeciesNames
= {"Electrons (core)", "Electrons (halo)", "Electrons (beam)", "Proton (core)", "Proton (halo)", "Proton (beam)", "Proton (pickup)",
   "Alpha (core)", "Alpha (halo)", "Helium II (core)", "Helium II (pickup)", "Hydrogen II (core)", "Hydrogen I (core)", "Hydrogen I (halp)",
   "Hydrogen I (beam)", "Helium I (core)"};

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

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Handy time conversion
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Seconds -> hours
#define HOURS(sec) ((sec) / 3600.0)

//! Seconds -> days
#define DAYS(sec) ((sec) / 86400.0)

//! Seconds -> years
#define YEARS(sec) ((sec) / 31536000.0)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Fluid flow routines
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Compute sound speed
\author Vladimir Florinski
\date 07/31/2019
\param[in] den Mass density
\param[in] pre Pressure
\param[in] ifl Index of the fluid
\return Sound speed
*/
SPECTRUM_DEVICE_FUNC inline double SoundSpeed(double den, double pre, unsigned int ifl = SPECIES_HYDROGEN_PLASMA_CORE)
{
   return sqrt(SpeciesPolytropicIndices[ifl] * pre / den);
};

/*!
\brief Compute square of the sound speed
\author Vladimir Florinski
\date 07/31/2019
\param[in] den Mass density
\param[in] pre Pressure
\param[in] ifl Index of the fluid
\return Sound speed squared
*/
SPECTRUM_DEVICE_FUNC inline double Sound2(double den, double pre, unsigned int ifl = SPECIES_HYDROGEN_PLASMA_CORE)
{
   return SpeciesPolytropicIndices[ifl] * pre / den;
};

/*!
\brief Compute Alfven speed
\author Vladimir Florinski
\date 07/31/2019
\param[in] den Mass density
\param[in] B2  Square of the magnetic field
\return Alfven speed
*/
SPECTRUM_DEVICE_FUNC inline double AlfvenSpeed(double den, double B2)
{
   return sqrt(B2 / den / M_4PI);
};

/*!
\brief Compute square of the Alfven speed
\author Vladimir Florinski
\date 07/31/2019
\param[in] den Mass density
\param[in] B2  Square of the magnetic field
\return Alfven speed squared
*/
SPECTRUM_DEVICE_FUNC inline double Alfven2(double den, double B2)
{
   return B2 / den / M_4PI;
};

/*!
\brief Compute fast magnetosonic speed
\author Vladimir Florinski
\date 07/31/2019
\param[in] Va2  Square of the Alfven speed
\param[in] Vax2 Square of the component of Alfven speed along direction of interest
\param[in] Cs2  Square of the sound speed
\return Fast magnetosonic speed
*/
SPECTRUM_DEVICE_FUNC inline double FastMagnetosonicSpeed(double Va2, double Vax2, double Cs2)
{
   return sqrt(0.5 * (Va2 + Cs2 + sqrt(Sqr(Va2 + Cs2) - 4.0 * Vax2 * Cs2)));
};

/*!
\brief Compute slow magnetosonic speed
\author Vladimir Florinski
\date 07/31/2019
\param[in] Va2  Square of the Alfven speed
\param[in] Vax2 Square of the component of Alfven speed along direction of interest
\param[in] Cs2  Square of the sound speed
\return Slow magnetosonic speed
*/
SPECTRUM_DEVICE_FUNC inline double SlowMagnetosonicSpeed(double Va2, double Vax2, double Cs2)
{
   return sqrt(0.5 * (Va2 + Cs2 - sqrt(Sqr(Va2 + Cs2) - 4.0 * Vax2 * Cs2)));
};

/*!
\brief Compute squares of both fast and slow magnetosonic speeds
\author Vladimir Florinski
\date 04/24/2024
\param[in]  Va2  Square of the Alfven speed
\param[in]  Vax2 Square of the component of Alfven speed along direction of interest
\param[in]  Cs2  Square of the sound speed
\param[out] Vf   Fast magnetosonic speed
\param[out] Vs   Slow magnetosonic speed
*/
SPECTRUM_DEVICE_FUNC inline void Magnetosonic2(double Va2, double Vax2, double Cs2, double& Vf2, double& Vs2)
{
   double Vd2 = sqrt(Sqr(Va2 + Cs2) - 4.0 * Vax2 * Cs2);
   Vf2 = 0.5 * (Va2 + Cs2 + Vd2);
   Vs2 = Vf2 - Vd2;
};

/*!
\brief Calculate the MHD energy density
\author Vladimir Florinski
\date 07/31/2019
\param[in] den Mass density
\param[in] u2  Square of velocity
\param[in] B2  Square of the magnetic field
\param[in] pre Gas pressure
\param[in] ifl Index of the fluid
\return Energy per unit volume
*/
SPECTRUM_DEVICE_FUNC inline double Energy(double den, double u2, double B2, double pre, unsigned int ifl = SPECIES_HYDROGEN_PLASMA_CORE)
{
   return den * u2 / 2.0 + pre / (SpeciesPolytropicIndices[ifl] - 1.0) + B2 / M_8PI;
};

/*!
\brief Calculate gas pressure from temperature
\author Vladimir Florinski
\date 07/15/2020
\param[in] den Mass density
\param[in] T   Temperature
\return Sound speed
*/
SPECTRUM_DEVICE_FUNC inline double Pressure(double den, double T)
{
   return den * kb_code * T;
};

/*!
\brief Calculate gas pressure from MHD energy density
\author Vladimir Florinski
\date 07/31/2019
\param[in] den Mass density
\param[in] u2  Square of velocity
\param[in] B2  Square of the magnetic field
\param[in] enr MHD energy density
\param[in] ifl Index of the fluid
\return gas pressure
*/
SPECTRUM_DEVICE_FUNC inline double Pressure(double den, double u2, double B2, double enr, unsigned int ifl = SPECIES_HYDROGEN_PLASMA_CORE)
{
   return (enr - den * u2 / 2.0 - B2 / M_8PI) * (SpeciesPolytropicIndices[ifl] - 1.0);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Kinetic/particle routines
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Relativistic factor
\author Vladimir Florinski
\date 07/15/2020
\param[in] vel Velocity
\return Lorentz factor
*/
SPECTRUM_DEVICE_FUNC inline double RelFactor(double vel)
{
   return 1.0 / sqrt(1.0 - Sqr(vel / c_code));
};

/*!
\brief Relativistic factor based on momentum
\author Vladimir Florinski
\date 12/29/2020
\param[in] mom Momentum
\return Lorentz factor
*/
SPECTRUM_DEVICE_FUNC inline double RelFactor1(double mom, unsigned int isp = SPECIES_PROTON_CORE)
{
   return sqrt(1.0 + Sqr(mom / (SpeciesMasses[isp] * c_code)));
};

/*!
\brief Relativistic factor based on parallel momentum and magnetic moment
\author Vladimir Florinski
\date 01/07/2021
\param[in] mom     Momentum (parallel)
\param[in] mag_mom Magnetic moment
\param[in] B       Magnetic field magnitude
\param[in] isp     Index of the species
\return Lorentz factor
*/
SPECTRUM_DEVICE_FUNC inline double RelFactor2(double mom, double mag_mom, double B, unsigned int isp = SPECIES_PROTON_CORE)
{
   return sqrt(1.0 + (2.0 * SpeciesMasses[isp] * mag_mom * B + Sqr(mom)) / Sqr(SpeciesMasses[isp] * c_code));
};

/*!
\brief Calculate particle momentum from its kinetic energy
\author Vladimir Florinski
\date 07/15/2020
\param[in] enr Kinetic energy
\param[in] isp Index of the species
\return Relativistic momentum
*/
SPECTRUM_DEVICE_FUNC inline double Mom(double enr, unsigned int isp = SPECIES_PROTON_CORE)
{
   return sqrt(enr * (enr + 2.0 * SpeciesMasses[isp] * c2_code)) / c_code;
};

/*!
\brief Calculate particle total energy from momentum magnitude
\author Vladimir Florinski
\date 07/15/2020
\param[in] mom Momentum
\param[in] isp Index of the species
\return Total relativistic energy
*/
SPECTRUM_DEVICE_FUNC inline double EnrTot(double mom, unsigned int isp = SPECIES_PROTON_CORE)
{
   return c_code * sqrt(Sqr(mom) + Sqr(SpeciesMasses[isp] * c_code));
};

/*!
\brief Calculate particle kinetic energy from momentum magnitude
\author Vladimir Florinski
\date 07/15/2020
\param[in] mom Momentum
\param[in] isp Index of the species
\return Relativistic kinetic energy
*/
SPECTRUM_DEVICE_FUNC inline double EnrKin(double mom, unsigned int isp = SPECIES_PROTON_CORE)
{
   return EnrTot(mom, isp) - SpeciesMasses[isp] * c2_code;
};

/*!
\brief Calculate particle velocity magnitude from momentum magnitude
\author Vladimir Florinski
\date 07/15/2020
\param[in] mom Momentum
\param[in] isp Index of the species
\return Velocity
*/
SPECTRUM_DEVICE_FUNC inline double Vel(double mom, unsigned int isp = SPECIES_PROTON_CORE)
{
   return mom * c2_code / EnrTot(mom, isp);
};

/*!
\brief Calculate particle velocity from momentum
\author Vladimir Florinski
\date 07/15/2020
\param[in] mom Momentum
\param[in] isp Index of the species
\return Velocity
*/
SPECTRUM_DEVICE_FUNC inline GeoVector Vel(const GeoVector& mom, unsigned int isp = SPECIES_PROTON_CORE)
{
   double mmag = mom.Norm();
   double vmag = Vel(mmag, isp);
   return (vmag / mmag) * mom;
};

/*!
\brief Calculate particle momentum from velocity
\author Vladimir Florinski
\date 12/24/2024
\param[in] vel Velocity
\param[in] isp Index of the species
\return Relativistic momentum
*/
SPECTRUM_DEVICE_FUNC inline GeoVector Mom(const GeoVector& vel, unsigned int isp = SPECIES_PROTON_CORE)
{
   return (RelFactor(vel.Norm()) * SpeciesMasses[isp]) * vel;
};

/*!
\brief Calculate particle rigidity
\author Vladimir Florinski
\date 07/14/2020
\param[in] mom Momentum
\param[in] isp Index of the species
\return Rigidity
*/
SPECTRUM_DEVICE_FUNC inline double Rigidity(double mom, unsigned int isp = SPECIES_PROTON_CORE)
{
   return mom * c_code / fabs(SpeciesCharges[isp]);
};

/*!
\brief Effective temperature (Maxwell) based on characteristic velocity
\author Vladimir Florinski
\date 07/14/2020
\param[in] v_th Thermal speed
\param[in] isp  Index of the species
\return Effective temperature
*/
SPECTRUM_DEVICE_FUNC inline double EffectiveTemperature(double v_th, unsigned int isp = SPECIES_PROTON_CORE)
{
   return SpeciesMasses[isp] * Sqr(v_th) / (2.0 * kb_code);
};

/*!
\brief Thermal speed (Maxwell) based on temperature
\author Vladimir Florinski
\date 07/14/2020
\param[in] T   Temperature
\param[in] isp Index of the species
\return Thermal speed
*/
SPECTRUM_DEVICE_FUNC inline double ThermalSpeed(double T, unsigned int isp = SPECIES_PROTON_CORE)
{
   return sqrt(2.0 * kb_code * T / SpeciesMasses[isp]);
};

/*!
\brief Cyclotron frequency
\author Vladimir Florinski
\date 07/14/2020
\param[in] vel Velocity
\param[in] B   Magnetic field magnitude
\param[in] isp Index of the species
\return Cyclotron frequency
*/
SPECTRUM_DEVICE_FUNC inline double CyclotronFrequency(double vel, double B, unsigned int isp = SPECIES_PROTON_CORE)
{
   return charge_mass_particle * SpeciesCharges[isp] * B / (RelFactor(vel) * SpeciesMasses[isp] * c_code);
};

/*!
\brief Larmor radius
\author Vladimir Florinski
\date 07/14/2020
\param[in] mom Momentum
\param[in] B   Magnetic field magnitude
\param[in] isp Index of the species
\return Larmor radius
*/
SPECTRUM_DEVICE_FUNC inline double LarmorRadius(double mom, double B, unsigned int isp = SPECIES_PROTON_CORE)
{
   return mom * c_code / (charge_mass_particle * fabs(SpeciesCharges[isp]) * B);
};

/*!
\brief Magnetic moment
\author Vladimir Florinski
\date 01/07/2021
\param[in] mom Momentum (perp. component)
\param[in] B   Magnetic field magnitude
\param[in] isp Index of the species
\return Magnetic moment
*/
SPECTRUM_DEVICE_FUNC inline double MagneticMoment(double mom, double B, unsigned int isp = SPECIES_PROTON_CORE)
{
   return Sqr(mom) / (2.0 * SpeciesMasses[isp] * B);
};

/*!
\brief Perpendicular momentum from magnetic moment
\author Vladimir Florinski
\date 03/30/2022
\param[in] mag_mom Magnetic moment
\param[in] B       Magnetic field magnitude
\param[in] isp     Index of the species
\return Momentum (perp. component)
*/
SPECTRUM_DEVICE_FUNC inline double PerpMomentum(double mag_mom, double B, unsigned int isp = SPECIES_PROTON_CORE)
{
   return sqrt(2.0 * SpeciesMasses[isp] * mag_mom * B);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Plasma physics routines
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma frequency
\author Vladimir Florinski
\date 01/21/2025
\param[in] den Density
\param[in] isp Index of the species
\return Plasma frequency
*/
SPECTRUM_DEVICE_FUNC inline double PlasmaFrequency(double den, unsigned int isp = SPECIES_PROTON_CORE)
{
   return sqrt(M_4PI * den) * charge_mass_particle * fabs(SpeciesCharges[isp]) / SpeciesMasses[isp];
};

/*!
\brief Collision frequency
\author Vladimir Florinski
\date 01/21/2025
\param[in] den Density
\param[in] T   Temperature
\param[in] Lam Coulomb logarithm
\param[in] isp Index of the species
\return Collision frequency
\note Our system of units (MHD centric) requires a coefficient to relate the number density unit to the inverse cube of distance unit.
*/
SPECTRUM_DEVICE_FUNC inline double CollisionFrequency(double den, double T, double Lam, unsigned int isp = SPECIES_PROTON_CORE)
{
   return 2.0 / 3.0 * M_SQRT2 / M_SQRTPI * Lam * SpeciesMasses[isp] / Cube(ThermalSpeed(T, isp) * unit_length_fluid) / unit_density_fluid
          * unit_mass_particle * Sqr(PlasmaFrequency(den, isp) * charge_mass_particle * SpeciesCharges[isp] / SpeciesMasses[isp]);
};

//! Print all units and constants
void PrintUnits(void);

};

#endif