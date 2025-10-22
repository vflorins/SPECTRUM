/*!
\file physics.hh
\brief Declares physical constants and conversion routines for various physical quantities
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_PHYSICS_HH
#define SPECTRUM_PHYSICS_HH

//#include <array>
//
//#include "config.h"
//
//#ifdef USE_GSL
//#include <gsl/gsl_const_cgsm.h>
//#endif

#include "common/vectors.hh"
#include "common/compiletime_lists.hh"


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
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] den Mass density
\param[in] prs Pressure
\return Sound speed
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double SoundSpeed(double den, double prs)
{
   return sqrt(specie.pt_idx * prs / den);
};


/*!
\brief Compute square of the sound speed
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] den Mass density
\param[in] prs Pressure
\return Sound speed squared
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double Sound2(double den, double prs)
{
   return specie.pt_idx * prs / den;
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
\author Lucius Schoenbaum
\date 09/13/2025
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
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] den Mass density
\param[in] u2  Square of velocity
\param[in] prs Gas pressure
\param[in] B2  Square of the magnetic field (if defined)
\return Energy per unit volume
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double Energy(double den, double u2, double prs, double B2 = 0.0)
{
   return den * u2 / 2.0 + prs / (specie.pt_idx - 1.0) + B2 / M_8PI;
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
\param[in] enr energy density
\param[in] B2  Square of the magnetic field (if defined)
\return gas pressure
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double Pressure(double den, double u2, double enr, double B2 = 0.0)
{
   return (enr - den * u2 / 2.0 - B2 / M_8PI) * (specie.pt_idx - 1.0);
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
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mom Momentum
\return Lorentz factor
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double RelFactor1(double mom)
{
   return sqrt(1.0 + Sqr(mom / (specie.mass * c_code)));
};


/*!
\brief Relativistic factor based on parallel momentum and magnetic moment
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mom     Momentum (parallel)
\param[in] mag_mom Magnetic moment
\param[in] B       Magnetic field magnitude
\return Lorentz factor
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double RelFactor2(double mom, double mag_mom, double B)
{
   return sqrt(1.0 + (2.0 * specie.mass * mag_mom * B + Sqr(mom)) / Sqr(specie.mass * c_code));
};


/*!
\brief Calculate particle momentum from its kinetic energy
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] enr Kinetic energy
\return Relativistic momentum
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double Mom(double enr)
{
   return sqrt(enr * (enr + 2.0 * specie.mass * c2_code)) / c_code;
};


/*!
\brief Calculate particle total energy from momentum magnitude
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mom Momentum
\return Total relativistic energy
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double EnrTot(double mom)
{
   return c_code * sqrt(Sqr(mom) + Sqr(specie.mass * c_code));
};


/*!
\brief Calculate particle kinetic energy from momentum magnitude
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mom Momentum
\return Relativistic kinetic energy
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double EnrKin(double mom)
{
   return EnrTot<specie>(mom) - specie.mass * c2_code;
};


/*!
\brief Calculate particle velocity magnitude from momentum magnitude
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mom Momentum
\return Velocity
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double Vel(double mom)
{
   return mom * c2_code / EnrTot<specie>(mom);
};


/*!
\brief Calculate particle velocity from momentum
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mom Momentum
\return Velocity
*/
template <Specie specie, CoordinateSystem Vel_Sys, CoordinateSystem Mom_Sys>
SPECTRUM_DEVICE_FUNC inline GeoVector Vel(const GeoVector& mom)
{
//   if constexpr (Vel_Sys == CoordinateSystem::)
   double mmag = mom.Norm();
   double vmag = Vel<specie>(mmag);
   return (vmag / mmag) * mom;
};


/*!
\brief Calculate particle momentum from velocity
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] vel Velocity
\return Relativistic momentum
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline GeoVector Mom(const GeoVector& vel)
{
   return (RelFactor(vel.Norm()) * specie.mass) * vel;
};


/*!
\brief Calculate particle rigidity
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mom Momentum
\return Rigidity
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double Rigidity(double mom)
{
   return mom * c_code / fabs(specie.charge);
};


/*!
\brief Effective temperature (Maxwell) based on characteristic velocity
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] v_th Thermal speed
\return Effective temperature
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double EffectiveTemperature(double v_th)
{
   return specie.mass * Sqr(v_th) / (2.0 * kb_code);
};


/*!
\brief Thermal speed (Maxwell) based on temperature
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] T   Temperature
\return Thermal speed
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double ThermalSpeed(double T)
{
   return sqrt(2.0 * kb_code * T / specie.mass);
};


/*!
\brief Cyclotron frequency
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] vel Velocity
\param[in] B   Magnetic field magnitude
\return Cyclotron frequency
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double CyclotronFrequency(double vel, double B)
{
   return charge_mass_particle * specie.charge * B / (RelFactor(vel) * specie.mass * c_code);
};


/*!
\brief Larmor radius
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mom Momentum
\param[in] B   Magnetic field magnitude
\return Larmor radius
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double LarmorRadius(double mom, double B)
{
   return mom * c_code / (charge_mass_particle * fabs(specie.charge) * B);
};


/*!
\brief Magnetic moment
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mom Momentum (perp. component)
\param[in] B   Magnetic field magnitude
\return Magnetic moment
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double MagneticMoment(double mom, double B)
{
   return Sqr(mom) / (2.0 * specie.mass * B);
};


/*!
\brief Perpendicular momentum from magnetic moment
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] mag_mom Magnetic moment
\param[in] B       Magnetic field magnitude
\return Momentum (perp. component)
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double PerpMomentum(double mag_mom, double B)
{
   return sqrt(2.0 * specie.mass * mag_mom * B);
};


//----------------------------------------------------------------------------------------------------------------------------------------------------
// Plasma physics routines
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\brief Plasma frequency
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] den Density
\return Plasma frequency
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double PlasmaFrequency(double den)
{
   return sqrt(M_4PI * den) * charge_mass_particle * fabs(specie.charge) / specie.mass;
};


/*!
\brief Collision frequency
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/13/2025
\param[in] den Density
\param[in] T   Temperature
\param[in] Lam Coulomb logarithm
\return Collision frequency
\note Our system of units (MHD centric) requires a coefficient to relate the number density unit to the inverse cube of distance unit.
*/
template <Specie specie>
SPECTRUM_DEVICE_FUNC inline double CollisionFrequency(double den, double T, double Lam)
{
   return 2.0 / 3.0 * M_SQRT2 / M_SQRTPI * Lam * specie.mass / Cube(ThermalSpeed<specie>(T) * unit_length_fluid) / unit_density_fluid
          * unit_mass_particle * Sqr(PlasmaFrequency<specie>(den) * charge_mass_particle * specie.charge / specie.mass);
};



//! Print all units and constants
void PrintUnits(void);

};

#endif