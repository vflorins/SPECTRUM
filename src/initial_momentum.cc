/*!
\file initial_momentum.cc
\brief Implements several classes to specify momentum initial conditions
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "initial_momentum.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumFixed methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
template <typename HConfig>
InitialMomentumFixed<HConfig>::InitialMomentumFixed(void)
                    : InitialBase(init_name, INITIAL_MOMENTUM | INITIAL_POINT)
{
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialMomentumFixed<HConfig>::InitialMomentumFixed(const InitialMomentumFixed& other)
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
template <typename HConfig>
void InitialMomentumFixed<HConfig>::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);

// Pre-assign "_mom" so that it never needs to change
   container.Read(initmom);
   _coords.Mom() = initmom;

};

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
template <typename HConfig>
void InitialMomentumFixed<HConfig>::EvaluateInitial(void)
{
// Nothing to do for FIXED_COORD - the value of "_coords.Mom()" was assigned in "Setupinitial()"
};



//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumBeam methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 04/04/2023
*/
template <typename HConfig>
InitialMomentumBeam<HConfig>::InitialMomentumBeam(void)
                   : InitialBase(init_name, INITIAL_MOMENTUM | INITIAL_POINT)
{
};

/*!
\author Vladimir Florinski
\date 04/04/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialMomentumBeam<HConfig>::InitialMomentumBeam(const InitialMomentumBeam& other)
                   : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_MOMENTUM);
   RAISE_BITS(_status, INITIAL_POINT);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/10/2025
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void InitialMomentumBeam<HConfig>::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);
   container.Read(p0);

   if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Focused) {
      _coords.Mom() = GeoVector(p0, 0.0, 0.0);
   }
   else if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Guiding) {
      _coords.Mom() = GeoVector(0.0, 0.0, p0);
   }
   else if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Lorentz){
// Nothing to do.
   }
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 08/10/2025
*/
template <typename HConfig>
void InitialMomentumBeam<HConfig>::EvaluateInitial(void)
{
   if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Lorentz){
      _coords.Mom() = p0 * axis;
   }
   else {
// For non-Lorentz trajectories the value of "_coords.Mom()" is assigned in "SetupInitial()"
   }
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumRing methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
template <typename HConfig>
InitialMomentumRing<HConfig>::InitialMomentumRing(void)
                   : InitialBase(init_name, INITIAL_MOMENTUM | INITIAL_CURVE)
{
};

/*!
\author Vladimir Florinski
\date 10/01/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialMomentumRing<HConfig>::InitialMomentumRing(const InitialMomentumRing& other)
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
template <typename HConfig>
void InitialMomentumRing<HConfig>::SetupInitial(bool construct)
{
   double theta0;

// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);
   container.Read(p0);
   container.Read(theta0);

   mu0 = cos(theta0);
   st0 = sin(theta0);

   if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Focused) {
      _coords.Mom() = GeoVector(p0, mu0, 0.0);
   }
   else if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Guiding) {
      _coords.Mom() = GeoVector(p0 * st0, 0.0, p0 * mu0);
   }
   else if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Lorentz){
// Nothing to do.
   }

};

/*!
\author Vladimir Florinski
\date 06/21/2021
*/
template <typename HConfig>
void InitialMomentumRing<HConfig>::EvaluateInitial(void)
{
   if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Lorentz){
      double phi = M_2PI * rng->GetUniform();
      GeoVector e1, e2;

// Component parallel to "axis"
      _coords.Mom() = p0 * mu0 * axis;

// Components normal to "axis"
      e1 = GetSecondUnitVec(axis);
      e2 = axis ^ e1;
      _coords.Mom() += p0 * st0 * (cos(phi) * e1 + sin(phi) * e2);
   }
   else {
// Nothing to do - the value of "_coords.Mom()" was assigned in "SetupInitial()"
      ;
   }
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumShell methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
template <typename HConfig>
InitialMomentumShell<HConfig>::InitialMomentumShell(void)
                    : InitialBase(init_name, INITIAL_MOMENTUM | INITIAL_SURFACE)
{
};

/*!
\author Vladimir Florinski
\date 10/01/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialMomentumShell<HConfig>::InitialMomentumShell(const InitialMomentumShell& other)
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
template <typename HConfig>
void InitialMomentumShell<HConfig>::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);
   container.Read(p0);

   if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Parker) {
      _coords.Mom() = GeoVector(p0, 0.0, 0.0);
   }
   else if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Fieldline) {
      _coords.Mom() = GeoVector(0.0, 0.0, p0);
   }
   else {
// Nothing to do.
   }
};

/*!
\author Vladimir Florinski
\date 06/21/2021
*/
template <typename HConfig>
void InitialMomentumShell<HConfig>::EvaluateInitial(void)
{

   if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Guiding
   || HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Focused) {

      double mu, st;

      mu = -1.0 + 2.0 * rng->GetUniform();
      st = sqrt(1.0 - Sqr(mu));

      if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Guiding) {
         _coords.Mom() = GeoVector(p0 * st, 0.0, p0 * mu);
      }
      else {
         _coords.Mom() = GeoVector(p0, mu, 0.0);
      }
   }
   else if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Lorentz) {
      double mu, st, phi;
      mu = -1.0 + 2.0 * rng->GetUniform();
      st = sqrt(1.0 - Sqr(mu));
      phi = M_2PI * rng->GetUniform();
      _coords.Mom() = GeoVector(p0 * st * cos(phi), p0 * st * sin(phi), p0 * mu);
   }
   else {
// Parker, Fieldline
// Nothing to do - the value of "_coords.Mom()" was assigned in "SetupInitial()"
   }
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumThickShell methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/21/2021
*/
template <typename HConfig>
InitialMomentumThickShell<HConfig>::InitialMomentumThickShell(void)
                         : InitialBase(init_name, INITIAL_MOMENTUM | INITIAL_VOLUME)
{
};

/*!
\author Vladimir Florinski
\date 10/01/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialMomentumThickShell<HConfig>::InitialMomentumThickShell(const InitialMomentumThickShell& other)
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
template <typename HConfig>
void InitialMomentumThickShell<HConfig>::SetupInitial(bool construct)
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
template <typename HConfig>
void InitialMomentumThickShell<HConfig>::EvaluateInitial(void)
{
   double p;

   if (log_bias) p = pow(10.0, p1 + (p2 - p1) * rng->GetUniform());
   else p = p1 + (p2 - p1) * rng->GetUniform();

   if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Parker) {
      _coords.Mom() = GeoVector(p, 0.0, 0.0);
   }
   else if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Fieldline) {
      _coords.Mom() = GeoVector(0.0, 0.0, p);
   }
   else if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Guiding || HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Focused) {
      double mu, st;

      mu = -1.0 + 2.0 * rng->GetUniform();
      st = sqrt(1.0 - Sqr(mu));
      if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Guiding) {
         _coords.Mom() = GeoVector(p * st, 0.0, p * mu);
      }
      else {
         _coords.Mom() = GeoVector(p, mu, 0.0);
      }
   }
   else if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Lorentz){
      double mu, st, phi;
      mu = -1.0 + 2.0 * rng->GetUniform();
      st = sqrt(1.0 - Sqr(mu));
      phi = M_2PI * rng->GetUniform();
      _coords.Mom() = GeoVector(p * st * cos(phi), p * st * sin(phi), p * mu);
   }
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialMomentumTable methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 12/27/2023
*/
template <typename HConfig>
InitialMomentumTable<HConfig>::InitialMomentumTable(void)
                    : InitialTable(init_name, INITIAL_MOMENTUM | INITIAL_POINT)
{
};

/*!
\author Juan G Alonso Guzman
\date 12/27/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialMomentumTable<HConfig>::InitialMomentumTable(const InitialMomentumTable& other)
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
template <typename HConfig>
void InitialMomentumTable<HConfig>::EvaluateInitial(void)
{
   if (random) {

// Generate random integer between 0 and initquant.size() - 1
      table_counter = rng->GetUniform() * initquant.size();

// Pull momentum in randomly selected place on the table
      _coords.Mom() = initquant[table_counter];
   }
   else {

// Pull next momentum on the table
      _coords.Mom() = initquant[table_counter++];

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
template <typename HConfig>
InitialMomentumMaxwell<HConfig>::InitialMomentumMaxwell(void)
                      : InitialBase(init_name, INITIAL_MOMENTUM | INITIAL_VOLUME)
{
};

/*!
\author Vladimir Florinski
\date 04/04/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialMomentumMaxwell<HConfig>::InitialMomentumMaxwell(const InitialMomentumMaxwell& other)
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
template <typename HConfig>
void InitialMomentumMaxwell<HConfig>::SetupInitial(bool construct)
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
template <typename HConfig>
void InitialMomentumMaxwell<HConfig>::EvaluateInitial(void)
{

   if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Fieldline) {
      _coords.Mom() = GeoVector(0.0, 0.0, p0);
      return;
   }

   double p_para = p0 + dp_para * rng->GetNormal();
   double p_perp = dp_perp * rng->GetRayleigh();

   if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Lorentz){
      double phi = M_2PI * rng->GetUniform();
      GeoVector e1, e2;
      e1 = GetSecondUnitVec(axis);
      e2 = axis ^ e1;
      _coords.Mom() = p_para * axis + p_perp * (cos(phi) * e1 + sin(phi) * e2);
   }
   else if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Guiding) {
      _coords.Mom() = GeoVector(p_perp, 0.0, p_para);
   }
   else {
      double p = sqrt(Sqr(p_para) + Sqr(p_perp));
      if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Focused){
         double mu = p_para / p;
         _coords.Mom() = GeoVector(p, mu, 0.0);
      }
      else if constexpr (HConfig::TrajectoryConfig::trajectoryid == TrajectoryId::Parker) {
         _coords.Mom() = GeoVector(p, 0.0, 0.0);
      }
   }
};

};
