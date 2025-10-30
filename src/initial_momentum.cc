/*!
\file initial_momentum.cc
\brief Implements several classes to specify momentum initial conditions
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "initial_momentum.hh"

namespace Spectrum {

#if (TRAJ_TYPE == TRAJ_LORENTZ) || (TRAJ_TYPE == TRAJ_FIELDLINE)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumFixed methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
InitialMomentumFixed::InitialMomentumFixed(void)
                    : InitialBase(init_name_momentum_fixed, 0, INITIAL_MOMENTUM | INITIAL_POINT)
{
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialMomentumFixed::InitialMomentumFixed(const InitialMomentumFixed& other)
                    : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_MOMENTUM);
   RAISE_BITS(_status, INITIAL_POINT);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialMomentumFixed::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);

#if INITIAL_MOM_FIXED_COORD == 0

// Pre-assign "_mom" so that it never needs to change
   container.Read(initmom);
   _mom = initmom;

#else

   double theta0, phi0;

   container.Read(p0);
   container.Read(theta0);
   container.Read(phi0);

   mu0 = cos(theta0);
   st0 = sin(theta0);
   sp0 = sin(phi0);
   cp0 = cos(phi0);

#endif
};

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
void InitialMomentumFixed::EvaluateInitial(void)
{
// Nothing to do for FIXED_COORD - the value of "_mom" was assigned in "Setupinitial()"
#if INITIAL_MOM_FIXED_COORD != 0

   GeoVector e1, e2;

// Component parallel to "axis"
   _mom = p0 * mu0 * axis;

// Components normal to "axis". The gyrophase is not uniquely defined, so the component normal to "axis" effectively has a "random" azimuthal angle. However, if a spatial position is not changing and the field is time-independent, then the azimuthal angle is also fixed.
   e1 = GetSecondUnitVec(axis);
   e2 = axis ^ e1;
   _mom += p0 * st0 * (cp0 * e1 + sp0 * e2);

#endif
};

#endif

#if (TRAJ_TYPE != TRAJ_PARKER) && (TRAJ_TYPE != TRAJ_PARKER_SOURCE) && (TRAJ_TYPE != TRAJ_FIELDLINE)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumBeam methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 04/04/2023
*/
InitialMomentumBeam::InitialMomentumBeam(void)
                   : InitialBase(init_name_momentum_beam, 0, INITIAL_MOMENTUM | INITIAL_POINT)
{
};

/*!
\author Vladimir Florinski
\date 04/04/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialMomentumBeam::InitialMomentumBeam(const InitialMomentumBeam& other)
                   : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_MOMENTUM);
   RAISE_BITS(_status, INITIAL_POINT);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Vladimir Florinski
\date 04/04/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialMomentumBeam::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);
   container.Read(p0);

#if (TRAJ_TYPE == TRAJ_GUIDING) || (TRAJ_TYPE == TRAJ_GUIDING_SCATT) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF_SCATT)
   _mom = GeoVector(0.0, 0.0, p0);
#elif TRAJ_TYPE == TRAJ_FOCUSED
   _mom = GeoVector(p0, 0.0, 0.0);
#endif
};

/*!
\author Vladimir Florinski
\date 04/04/2023
*/
void InitialMomentumBeam::EvaluateInitial(void)
{
// For non-Lorentz trajectories the value of "_mom" is assigned in "SetupInitial()"
#if (TRAJ_TYPE == TRAJ_LORENTZ)
   _mom = p0 * axis;
#endif
};

#endif

#if (TRAJ_TYPE != TRAJ_PARKER) && (TRAJ_TYPE != TRAJ_PARKER_SOURCE) && (TRAJ_TYPE != TRAJ_FIELDLINE)

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumRing methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
InitialMomentumRing::InitialMomentumRing(void)
                   : InitialBase(init_name_momentum_ring, 0, INITIAL_MOMENTUM | INITIAL_CURVE)
{
};

/*!
\author Vladimir Florinski
\date 10/01/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialMomentumRing::InitialMomentumRing(const InitialMomentumRing& other)
                   : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_MOMENTUM);
   RAISE_BITS(_status, INITIAL_CURVE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Vladimir Florinski
\date 10/01/2021
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialMomentumRing::SetupInitial(bool construct)
{
   double theta0;

// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);
   container.Read(p0);
   container.Read(theta0);

   mu0 = cos(theta0);
   st0 = sin(theta0);

#if (TRAJ_TYPE == TRAJ_GUIDING) || (TRAJ_TYPE == TRAJ_GUIDING_SCATT) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF_SCATT)
   _mom = GeoVector(p0 * st0, 0.0, p0 * mu0);
#elif TRAJ_TYPE == TRAJ_FOCUSED
   _mom = GeoVector(p0, mu0, 0.0);
#endif
};

/*!
\author Vladimir Florinski
\date 06/21/2021
*/
void InitialMomentumRing::EvaluateInitial(void)
{
#if (TRAJ_TYPE == TRAJ_GUIDING) || (TRAJ_TYPE == TRAJ_GUIDING_SCATT) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF_SCATT) || (TRAJ_TYPE == TRAJ_FOCUSED)
// Nothing to do - the value of "_mom" was assigned in "SetupInitial()"
#elif (TRAJ_TYPE == TRAJ_LORENTZ)

   double phi = M_2PI * rng->GetUniform();
   GeoVector e1, e2;

// Component parallel to "axis"
   _mom = p0 * mu0 * axis;

// Components normal to "axis"
   e1 = GetSecondUnitVec(axis);
   e2 = axis ^ e1;
   _mom += p0 * st0 * (cos(phi) * e1 + sin(phi) * e2);

#endif
};

#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumShell methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
InitialMomentumShell::InitialMomentumShell(void)
                    : InitialBase(init_name_momentum_shell, 0, INITIAL_MOMENTUM | INITIAL_SURFACE)
{
};

/*!
\author Vladimir Florinski
\date 10/01/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialMomentumShell::InitialMomentumShell(const InitialMomentumShell& other)
                    : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_MOMENTUM);
   RAISE_BITS(_status, INITIAL_SURFACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Vladimir Florinski
\date 10/01/2021
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialMomentumShell::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);
   container.Read(p0);

#if (TRAJ_TYPE == TRAJ_PARKER) || (TRAJ_TYPE == TRAJ_PARKER_SOURCE)
   _mom = GeoVector(p0, 0.0, 0.0);
#elif TRAJ_TYPE == TRAJ_FIELDLINE
   _mom = GeoVector(0.0, 0.0, p0);
#endif
};

/*!
\author Vladimir Florinski
\date 06/21/2021
*/
void InitialMomentumShell::EvaluateInitial(void)
{
#if (TRAJ_TYPE == TRAJ_PARKER) || (TRAJ_TYPE == TRAJ_PARKER_SOURCE) || (TRAJ_TYPE == TRAJ_FIELDLINE)
// Nothing to do - the value of "_mom" was assigned in "SetupInitial()"
#elif (TRAJ_TYPE == TRAJ_GUIDING) || (TRAJ_TYPE == TRAJ_GUIDING_SCATT) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF_SCATT) || (TRAJ_TYPE == TRAJ_FOCUSED)

   double mu, st;

   mu = -1.0 + 2.0 * rng->GetUniform();
   st = sqrt(1.0 - Sqr(mu));

#if TRAJ_TYPE == TRAJ_FOCUSED
   _mom = GeoVector(p0, mu, 0.0);
#else
   _mom = GeoVector(p0 * st, 0.0, p0 * mu);
#endif

#elif TRAJ_TYPE == TRAJ_LORENTZ

   double mu, st, phi;
   mu = -1.0 + 2.0 * rng->GetUniform();
   st = sqrt(1.0 - Sqr(mu));
   phi = M_2PI * rng->GetUniform();
   _mom = GeoVector(p0 * st * cos(phi), p0 * st * sin(phi), p0 * mu);

#endif
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumThickShell methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/21/2021
*/
InitialMomentumThickShell::InitialMomentumThickShell(void)
                         : InitialBase(init_name_momentum_thickshell, 0, INITIAL_MOMENTUM | INITIAL_VOLUME)
{
};

/*!
\author Vladimir Florinski
\date 10/01/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialMomentumThickShell::InitialMomentumThickShell(const InitialMomentumThickShell& other)
                         : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_MOMENTUM);
   RAISE_BITS(_status, INITIAL_VOLUME);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Vladimir Florinski
\date 10/01/2021
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialMomentumThickShell::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);
   container.Read(p1);
   container.Read(p2);
   container.Read(log_bias);

   if (log_bias) {
      p1 = log10(p1);
      p2 = log10(p2);
   };
};

/*!
\author Vladimir Florinski
\date 06/21/2021
*/
void InitialMomentumThickShell::EvaluateInitial(void)
{
   double p;

   if (log_bias) p = pow(10.0, p1 + (p2 - p1) * rng->GetUniform());
   else p = p1 + (p2 - p1) * rng->GetUniform();
   
#if (TRAJ_TYPE == TRAJ_PARKER) || (TRAJ_TYPE == TRAJ_PARKER_SOURCE)
   _mom = GeoVector(p, 0.0, 0.0);
#elif TRAJ_TYPE == TRAJ_FIELDLINE
   _mom = GeoVector(0.0, 0.0, p);
#elif (TRAJ_TYPE == TRAJ_GUIDING) || (TRAJ_TYPE == TRAJ_GUIDING_SCATT) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF_SCATT) || (TRAJ_TYPE == TRAJ_FOCUSED)

   double mu, st;

   mu = -1.0 + 2.0 * rng->GetUniform();
   st = sqrt(1.0 - Sqr(mu));

#if TRAJ_TYPE == TRAJ_FOCUSED
   _mom = GeoVector(p, mu, 0.0);
#else
   _mom = GeoVector(p * st, 0.0, p * mu);
#endif

#elif TRAJ_TYPE == TRAJ_LORENTZ

   double mu, st, phi;
   mu = -1.0 + 2.0 * rng->GetUniform();
   st = sqrt(1.0 - Sqr(mu));
   phi = M_2PI * rng->GetUniform();
   _mom = GeoVector(p * st * cos(phi), p * st * sin(phi), p * mu);

#endif
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumTable methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 12/27/2023
*/
InitialMomentumTable::InitialMomentumTable(void)
                    : InitialTable(init_name_momentum_table, 0, INITIAL_MOMENTUM | INITIAL_POINT)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/27/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialMomentumTable::InitialMomentumTable(const InitialMomentumTable& other)
                    : InitialTable(other)
{
   RAISE_BITS(_status, INITIAL_MOMENTUM);
   RAISE_BITS(_status, INITIAL_POINT);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 12/27/2023
*/
void InitialMomentumTable::EvaluateInitial(void)
{
   if (random) {

// Generate random integer between 0 and initquant.size() - 1
      table_counter = rng->GetUniform() * initquant.size();

// Pull momentum in randomly selected place on the table
      _mom = initquant[table_counter];
   }
   else {

// Pull next momentum on the table
      _mom = initquant[table_counter++];

// If all momenta have been sampled, reset the counter
      if (table_counter == initquant.size()) table_counter = 0;
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumMaxwell methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 04/04/2023
*/
InitialMomentumMaxwell::InitialMomentumMaxwell(void)
                      : InitialBase(init_name_momentum_maxwell, 0, INITIAL_MOMENTUM | INITIAL_VOLUME)
{
};

/*!
\author Vladimir Florinski
\date 04/04/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialMomentumMaxwell::InitialMomentumMaxwell(const InitialMomentumMaxwell& other)
                      : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_MOMENTUM);
   RAISE_BITS(_status, INITIAL_VOLUME);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Vladimir Florinski
\date 04/04/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialMomentumMaxwell::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);
   container.Read(p0);
   container.Read(dp_para);
   container.Read(dp_perp);

// The RNG calculation is based on a standard deviation, which is smaller than the thermal speeed by a factor of sqrt(2).
   dp_para /= M_SQRT2;
   dp_perp /= M_SQRT2;
};

/*!
\author Vladimir Florinski
\date 04/04/2023
*/
void InitialMomentumMaxwell::EvaluateInitial(void)
{
#if TRAJ_TYPE == TRAJ_FIELDLINE
   _mom = GeoVector(0.0, 0.0, p0);
   return;
#endif

   double p_para = p0 + dp_para * rng->GetNormal();
   double p_perp = dp_perp * rng->GetRayleigh();

#if TRAJ_TYPE == TRAJ_LORENTZ

   double phi = M_2PI * rng->GetUniform();
   GeoVector e1, e2;
   e1 = GetSecondUnitVec(axis);
   e2 = axis ^ e1;
   _mom = p_para * axis + p_perp * (cos(phi) * e1 + sin(phi) * e2);

#elif (TRAJ_TYPE == TRAJ_GUIDING) || (TRAJ_TYPE == TRAJ_GUIDING_SCATT) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF) || (TRAJ_TYPE == TRAJ_GUIDING_DIFF_SCATT)
   _mom = GeoVector(p_perp, 0.0, p_para);

#else

   double p = sqrt(Sqr(p_para) + Sqr(p_perp));

#if TRAJ_TYPE == TRAJ_FOCUSED
   double mu = p_para / p;
   _mom = GeoVector(p, mu, 0.0);
#elif (TRAJ_TYPE == TRAJ_PARKER) || (TRAJ_TYPE == TRAJ_PARKER_SOURCE)
   _mom = GeoVector(p, 0.0, 0.0);
#endif

#endif
};

};
