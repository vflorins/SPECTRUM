/*!
\file initial_space.cc
\brief Implements several classes to specify spatial initial conditions
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "initial_space.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceFixed methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
template <typename HConfig>
InitialSpaceFixed<HConfig>::InitialSpaceFixed(void)
                 : InitialBase(init_name, INITIAL_SPACE | INITIAL_POINT)
{
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialSpaceFixed<HConfig>::InitialSpaceFixed(const InitialSpaceFixed& other)
                 : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
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
void InitialSpaceFixed<HConfig>::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);
   container.Read(initpos);

// Pre-assign "_coords.Pos()" so that it never needs to change
   _coords.Pos() = initpos;
};

/*!
\author Vladimir Florinski
\date 12/16/2021
*/
template <typename HConfig>
void InitialSpaceFixed<HConfig>::EvaluateInitial(void)
{
// Nothing to do - the value of "_coords.Pos()" was assigned in "Setupinitial()"
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceLine methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
template <typename HConfig>
InitialSpaceLine<HConfig>::InitialSpaceLine(void)
                : InitialBase(init_name, INITIAL_SPACE | INITIAL_CURVE)
{
};

/*!
\author Vladimir Florinski
\date 09/16/2022
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
InitialSpaceLine<HConfig>::InitialSpaceLine(const std::string& name_in, uint16_t status_in)
                : InitialBase(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialSpaceLine<HConfig>::InitialSpaceLine(const InitialSpaceLine& other)
                : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_CURVE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Vladimir Florinski
\date 04/04/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void InitialSpaceLine<HConfig>::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);

   int n_intervals;
   container.Read(startpos);
   container.Read(endpos);
   container.Read(n_intervals);

   if (n_intervals <= 0) randompos = true;
   else {
      randompos = false;
      increment = (endpos - startpos) / (double)n_intervals;
      _coords.Pos() = startpos - 0.5 * increment;
   };
};

/*!
\author Vladimir Florinski
\date 04/04/2023
*/
template <typename HConfig>
void InitialSpaceLine<HConfig>::EvaluateInitial(void)
{
   if (randompos) _coords.Pos() = startpos + (endpos - startpos) * rng->GetUniform();
   else _coords.Pos() += increment;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceCircle methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 10/11/2023
*/
template <typename HConfig>
InitialSpaceCircle<HConfig>::InitialSpaceCircle(void)
                  : InitialBase(init_name, INITIAL_SPACE | INITIAL_CURVE)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/11/2023
\param[in] name_in   Readable name of the class
\param[in] status_in Initial status
*/
template <typename HConfig>
InitialSpaceCircle<HConfig>::InitialSpaceCircle(const std::string& name_in, uint16_t status_in)
                  : InitialBase(name_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 10/11/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialSpaceCircle<HConfig>::InitialSpaceCircle(const InitialSpaceCircle& other)
                  : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_CURVE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 10/11/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void InitialSpaceCircle<HConfig>::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);

   GeoVector normal;
   double radius;
   container.Read(origin);
   container.Read(normal);
   container.Read(radius);
   radius_x = radius * GetSecondUnitVec(normal.Normalize());
   radius_y = normal ^ radius_x;
};

/*!
\author Juan G Alonso Guzman
\date 10/11/2023
*/
template <typename HConfig>
void InitialSpaceCircle<HConfig>::EvaluateInitial(void)
{
   double phi;

   phi = M_2PI * rng->GetUniform();
   _coords.Pos() = origin + radius_x * cos(phi) + radius_y * sin(phi);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceBox methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 09/16/2022
*/
template <typename HConfig>
InitialSpaceBox<HConfig>::InitialSpaceBox(void)
               : InitialSpaceLine(init_name, INITIAL_SPACE | INITIAL_VOLUME)
{
};

/*!
\author Vladimir Florinski
\date 09/16/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialSpaceBox<HConfig>::InitialSpaceBox(const InitialSpaceBox& other)
               : InitialSpaceLine(other)
{
   LOWER_BITS(_status, INITIAL_CURVE);
   RAISE_BITS(_status, INITIAL_VOLUME);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Vladimir Florinski
\date 09/16/2022
*/
template <typename HConfig>
void InitialSpaceBox<HConfig>::EvaluateInitial(void)
{
   for (auto xyz = 0; xyz < 3; xyz++) _coords.Pos()[xyz] = startpos[xyz] + (endpos[xyz] - startpos[xyz]) * rng->GetUniform();
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceSphere methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
template <typename HConfig>
InitialSpaceSphere<HConfig>::InitialSpaceSphere(void)
                  : InitialBase(init_name, INITIAL_SPACE | INITIAL_SURFACE)
{
};

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
template <typename HConfig>
InitialSpaceSphere<HConfig>::InitialSpaceSphere(const std::string& name_in, uint16_t status_in)
                  : InitialBase(name_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialSpaceSphere<HConfig>::InitialSpaceSphere(const InitialSpaceSphere& other)
                  : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_SURFACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void InitialSpaceSphere<HConfig>::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);
   container.Read(origin);
   container.Read(radius);
};

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
template <typename HConfig>
void InitialSpaceSphere<HConfig>::EvaluateInitial(void)
{
   double st, ct, phi;

// The polar angle _cosine_ is uniformly distributed to achieve equal point density everywhere.
   ct = -1.0 + 2.0 * rng->GetUniform();
   st = sqrt(1.0 - Sqr(ct));
   phi = M_2PI * rng->GetUniform();

   _coords.Pos() = origin + GeoVector(radius * st * cos(phi), radius * st * sin(phi), radius * ct);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceSphereSector methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
template <typename HConfig>
InitialSpaceSphereSector<HConfig>::InitialSpaceSphereSector(void)
                        : InitialSpaceSphere(init_name, INITIAL_SPACE | INITIAL_SURFACE)
{
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialSpaceSphereSector<HConfig>::InitialSpaceSphereSector(const InitialSpaceSphereSector& other)
                        : InitialSpaceSphere(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_SURFACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/28/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void InitialSpaceSphereSector<HConfig>::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialSpaceSphere::SetupInitial(false);

   container.Read(theta1);
   if (theta1 < 0.0) theta1 = 0.0;
   costheta1 = cos(theta1);

   container.Read(theta2);
   if (theta2 > M_PI) theta2 = M_PI;
   if (theta2 < theta1) theta2 = theta1;
   costheta2 = cos(theta2);

   container.Read(phi1);
   // if (phi1 < 0.0) phi1 = 0.0;

   container.Read(phi2);
   if (phi2 > M_2PI) phi2 = M_2PI;
   if (phi2 < phi1) phi2 = phi1;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/28/2022
*/
template <typename HConfig>
void InitialSpaceSphereSector<HConfig>::EvaluateInitial(void)
{
   double st, ct, phi;

// The polar angle _cosine_ is uniformly distributed to achieve equal point density everywhere.
   ct = costheta1 + (costheta2 - costheta1) * rng->GetUniform();
   st = sqrt(1.0 - Sqr(ct));
   phi = phi1 + (phi2 - phi1) * rng->GetUniform();

   _coords.Pos() = origin + GeoVector(radius * st * cos(phi), radius * st * sin(phi), radius * ct);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceRankineHalfBody methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
template <typename HConfig>
InitialSpaceRankineHalfBody<HConfig>::InitialSpaceRankineHalfBody(void)
                           : InitialBase(init_name, INITIAL_SPACE | INITIAL_SURFACE)
{
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialSpaceRankineHalfBody<HConfig>::InitialSpaceRankineHalfBody(const InitialSpaceRankineHalfBody& other)
                           : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_SURFACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/28/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void InitialSpaceRankineHalfBody<HConfig>::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);

// Read parameters
   container.Read(origin);
   container.Read(z_nose);
   container.Read(radius);

// Find limiting cos(theta)
   cos_lim = (2.0 * Sqr(z_nose) - Sqr(radius)) / Sqr(radius);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/28/2022
*/
template <typename HConfig>
void InitialSpaceRankineHalfBody<HConfig>::EvaluateInitial(void)
{
   double costheta, r, z, s, phi;

// The polar angle _cosine_ is uniformly distributed from -1 to cos_lim which leads to a non-uniform distribution along the surface because the surface is non-spherical
   costheta = cos_lim + (1.0 - cos_lim) * rng->GetUniform();
   r = z_nose * sqrt(2.0 / (1.0 + costheta));
   z = r * costheta;
   s = sqrt(Sqr(r) - Sqr(z));
   phi = M_2PI * rng->GetUniform();

   _coords.Pos() = origin + GeoVector(s * cos(phi), s * sin(phi), z);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceTable methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 04/09/2023
*/
template <typename HConfig>
InitialSpaceTable<HConfig>::InitialSpaceTable(void)
                 : InitialTable(init_name, INITIAL_SPACE | INITIAL_POINT)
{
};

/*!
\author Juan G Alonso Guzman
\date 04/09/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialSpaceTable<HConfig>::InitialSpaceTable(const InitialSpaceTable& other)
                 : InitialTable(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_POINT);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 12/27/2023
*/
template <typename HConfig>
void InitialSpaceTable<HConfig>::EvaluateInitial(void)
{
   if (random) {

// Generate random integer between 0 and initquant.size() - 1
      table_counter = rng->GetUniform() * initquant.size();

// Pull position in randomly selected place on the table
      _coords.Pos() = initquant[table_counter];
   }
   else {

// Pull next position on the table
      _coords.Pos() = initquant[table_counter++];

// If all positions have been sampled, reset the counter
      if (table_counter == initquant.size()) table_counter = 0;
   };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceCylinder methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
*/
template <typename HConfig>
InitialSpaceCylinder<HConfig>::InitialSpaceCylinder(void)
                    : InitialBase(init_name, INITIAL_SPACE | INITIAL_SURFACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
*/
template <typename HConfig>
InitialSpaceCylinder<HConfig>::InitialSpaceCylinder(const std::string& name_in, uint16_t status_in)
                    : InitialBase(name_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialSpaceCylinder<HConfig>::InitialSpaceCylinder(const InitialSpaceCylinder& other)
                    : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_SURFACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void InitialSpaceCylinder<HConfig>::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialBase::SetupInitial(false);
   container.Read(origin);
   container.Read(height);
   container.Read(radius_x);

   radius_y = UnitVec(height) ^ radius_x;
};

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
*/
template <typename HConfig>
void InitialSpaceCylinder<HConfig>::EvaluateInitial(void)
{
   double z, phi;

   z = rng->GetUniform();
   phi = M_2PI * rng->GetUniform();

   _coords.Pos() = origin + radius_x * cos(phi) + radius_y * sin(phi) + height * z;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceCylinderSector methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
*/
template <typename HConfig>
InitialSpaceCylinderSector<HConfig>::InitialSpaceCylinderSector(void)
                          : InitialSpaceCylinder(init_name, INITIAL_SPACE | INITIAL_SURFACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
template <typename HConfig>
InitialSpaceCylinderSector<HConfig>::InitialSpaceCylinderSector(const InitialSpaceCylinderSector& other)
                          : InitialSpaceCylinder(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_SURFACE);
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
template <typename HConfig>
void InitialSpaceCylinderSector<HConfig>::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) InitialSpaceCylinder::SetupInitial(false);

   container.Read(phi1);
   // if (phi1 < 0.0) phi1 = 0.0;

   container.Read(phi2);
   if (phi2 > M_2PI) phi2 = M_2PI;
   if (phi2 < phi1) phi2 = phi1;
};

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
*/
template <typename HConfig>
void InitialSpaceCylinderSector<HConfig>::EvaluateInitial(void)
{
   double z, phi;

   z = rng->GetUniform();
   phi = phi1 + (phi2 - phi1) * rng->GetUniform();

   _coords.Pos() = origin + radius_x * cos(phi) + radius_y * sin(phi) + height * z;
};

};
