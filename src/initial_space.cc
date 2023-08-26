/*!
\file initial_space.cc
\brief Implements several classes to specify spatial initial conditions
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "initial_space.hh"
#include <iostream> // std::cout
#include <fstream> // std::ifsteam

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceFixed methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
InitialSpaceFixed::InitialSpaceFixed(void)
                 : InitialBase(init_name_space_fixed, 0, INITIAL_SPACE | INITIAL_POINT)
{
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialSpaceFixed::InitialSpaceFixed(const InitialSpaceFixed& other)
                 : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_POINT);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialSpaceFixed::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) InitialBase::SetupInitial(false);
   container.Read(initpos.Data());

// Pre-assign "_pos" so that it never needs to change
   _pos = initpos;
};

/*!
\author Vladimir Florinski
\date 12/16/2021
*/
void InitialSpaceFixed::EvaluateInitial(void)
{
// Nothing to do - the value of "_pos" was assigned in "Setupinitial()"
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceLine methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
InitialSpaceLine::InitialSpaceLine(void)
                : InitialBase(init_name_space_line, 0, INITIAL_SPACE | INITIAL_CURVE)
{
};

/*!
\author Vladimir Florinski
\date 09/16/2022
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
InitialSpaceLine::InitialSpaceLine(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                : InitialBase(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialSpaceLine::InitialSpaceLine(const InitialSpaceLine& other)
                : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_CURVE);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Vladimir Florinski
\date 04/04/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialSpaceLine::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) InitialBase::SetupInitial(false);

   int n_intervals;
   container.Read(startpos.Data());
   container.Read(endpos.Data());
   container.Read(&n_intervals);

   if(n_intervals < 0) randompos = true;
   else {
      randompos = false;
      increment = (endpos - startpos) / (double)n_intervals;
      _pos = startpos - 0.5 * increment;
   };
};

/*!
\author Vladimir Florinski
\date 04/04/2023
*/
void InitialSpaceLine::EvaluateInitial(void)
{
   if(randompos) _pos = startpos + (endpos - startpos) * rng->GetUniform();
   else _pos += increment;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceCube methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 09/16/2022
*/
InitialSpaceCube::InitialSpaceCube(void)
                : InitialSpaceLine(init_name_space_cube, 0, INITIAL_SPACE | INITIAL_VOLUME)
{
};

/*!
\author Vladimir Florinski
\date 09/16/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialSpaceCube::InitialSpaceCube(const InitialSpaceCube& other)
                : InitialSpaceLine(other)
{
   LOWER_BITS(_status, INITIAL_CURVE);
   RAISE_BITS(_status, INITIAL_VOLUME);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Vladimir Florinski
\date 09/16/2022
*/
void InitialSpaceCube::EvaluateInitial(void)
{
   for(auto xyz = 0; xyz < 3; xyz++) _pos[xyz] = startpos[xyz] + (endpos[xyz] - startpos[xyz]) * rng->GetUniform();
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceSphere methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
InitialSpaceSphere::InitialSpaceSphere(void)
                  : InitialBase(init_name_space_sphere, 0, INITIAL_SPACE | INITIAL_SURFACE)
{
};

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
InitialSpaceSphere::InitialSpaceSphere(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                  : InitialBase(name_in, specie_in, status_in)
{
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialSpaceSphere::InitialSpaceSphere(const InitialSpaceSphere& other)
                  : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_SURFACE);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialSpaceSphere::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) InitialBase::SetupInitial(false);
   container.Read(origin.Data());
   container.Read(&radius);
};

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
void InitialSpaceSphere::EvaluateInitial(void)
{
   double st, ct, phi;

// The polar angle _cosine_ is uniformly distributed to achieve equal point density everywhere.
   ct = -1.0 + 2.0 * rng->GetUniform();
   st = sqrt(1.0 - Sqr(ct));
   phi = twopi * rng->GetUniform();

   _pos = origin + GeoVector(radius * st * cos(phi), radius * st * sin(phi), radius * ct);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceSphereSector methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
InitialSpaceSphereSector::InitialSpaceSphereSector(void)
                        : InitialSpaceSphere(init_name_space_sphere_sector, 0, INITIAL_SPACE | INITIAL_SURFACE)
{
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialSpaceSphereSector::InitialSpaceSphereSector(const InitialSpaceSphereSector& other)
                        : InitialSpaceSphere(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_SURFACE);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/28/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialSpaceSphereSector::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) InitialSpaceSphere::SetupInitial(false);

   container.Read(&theta1);
   if(theta1 < 0.0) theta1 = 0.0;
   costheta1 = cos(theta1);

   container.Read(&theta2);
   if(theta2 > M_PI) theta2 = M_PI;
   if(theta2 < theta1) theta2 = theta1;
   costheta2 = cos(theta2);

   container.Read(&phi1);
   // if(phi1 < 0.0) phi1 = 0.0;

   container.Read(&phi2);
   if(phi2 > twopi) phi2 = twopi;
   if(phi2 < phi1) phi2 = phi1;
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/28/2022
*/
void InitialSpaceSphereSector::EvaluateInitial(void)
{
   double st, ct, phi;

// The polar angle _cosine_ is uniformly distributed to achieve equal point density everywhere.
   ct = costheta1 + (costheta2 - costheta1) * rng->GetUniform();
   st = sqrt(1.0 - Sqr(ct));
   phi = phi1 + (phi2 - phi1) * rng->GetUniform();

   _pos = origin + GeoVector(radius * st * cos(phi), radius * st * sin(phi), radius * ct);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceRankineHalfBody methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 06/14/2021
*/
InitialSpaceRankineHalfBody::InitialSpaceRankineHalfBody(void)
                           : InitialBase(init_name_space_rankine_half_body, 0, INITIAL_SPACE | INITIAL_SURFACE)
{
};

/*!
\author Vladimir Florinski
\date 09/30/2021
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialSpaceRankineHalfBody::InitialSpaceRankineHalfBody(const InitialSpaceRankineHalfBody& other)
                           : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_SURFACE);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/28/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialSpaceRankineHalfBody::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) InitialBase::SetupInitial(false);

// Read parameters
   container.Read(origin.Data());
   container.Read(axis.Data());
   container.Read(&z_nose);
   container.Read(&radius);

// Find limiting cos(theta)
   cos_lim = (2.0 * Sqr(z_nose) - Sqr(radius)) / Sqr(radius);
};

/*!
\author Juan G Alonso Guzman
\author Vladimir Florinski
\date 06/28/2022
*/
void InitialSpaceRankineHalfBody::EvaluateInitial(void)
{
   double costheta, r, z, s, phi;

// The polar angle _cosine_ is uniformly distributed from -1 to cos_lim which leads to a non-uniform distribution along the surface because the surface is non-spherical
   costheta = cos_lim + (1.0 - cos_lim) * rng->GetUniform();
   r = z_nose * sqrt(2.0 / (1.0 + costheta));
   z = r * costheta;
   s = sqrt(Sqr(r) - Sqr(z));
   phi = twopi * rng->GetUniform();

   _pos = origin + GeoVector(s * cos(phi), s * sin(phi), z);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceTable methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 04/09/2023
*/
InitialSpaceTable::InitialSpaceTable(void)
                 : InitialBase(init_name_space_table, 0, INITIAL_SPACE | INITIAL_POINT)
{
};

/*!
\author Juan G Alonso Guzman
\date 04/09/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialSpaceTable::InitialSpaceTable(const InitialSpaceTable& other)
                 : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_POINT);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 04/09/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialSpaceTable::SetupInitial(bool construct)
{
   std::string initpos_file_name, coord_type;
   std::ifstream initpos_file;
   int i, size, coord;
   GeoVector entry;

// The parent version must be called explicitly if not constructing
   if(!construct) InitialBase::SetupInitial(false);
   container.Read(&initpos_file_name);

// Input initial positions from file
   initpos_file.open(initpos_file_name.c_str());

   initpos_file >> coord_type;
   if(coord_type == "RTP") {
      std::cerr << "Reading initial positions file in spherical coordinates." << std::endl;
      coord = 1;
   }
   else if(coord_type == "XYZ") {
      std::cerr << "Reading initial positions file in cartesian coordinates." << std::endl;
      coord = 0;
   }
   else {
      std::cerr << "Reading initial positions file in unrecognized coordinate type. "
                << "Defaulting to cartesian coordinates." << std::endl;
      coord = 0;
   };
   initpos_file >> size;
   initpos.resize(size);
   for(i = 0; i < size; i++) {
      initpos_file >> entry[0];
      initpos_file >> entry[1];
      initpos_file >> entry[2];
      if(coord) entry.RTP_XYZ();
      initpos[i] = entry;
   };
   initpos_file.close();

// Initialize table counter
   table_counter = 0;
};

/*!
\author Juan G Alonso Guzman
\date 04/09/2023
*/
void InitialSpaceTable::EvaluateInitial(void)
{
// Pull next position on the table
   _pos = initpos[table_counter++];
// If all positions have been sampled, reset the counter
   if(table_counter == initpos.size()) table_counter = 0;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceCylinder methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
*/
InitialSpaceCylinder::InitialSpaceCylinder(void)
                    : InitialBase(init_name_space_cylinder, 0, INITIAL_SPACE | INITIAL_SURFACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
*/
InitialSpaceCylinder::InitialSpaceCylinder(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
                    : InitialBase(name_in, specie_in, status_in)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialSpaceCylinder::InitialSpaceCylinder(const InitialSpaceCylinder& other)
                    : InitialBase(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_SURFACE);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialSpaceCylinder::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) InitialBase::SetupInitial(false);
   container.Read(origin.Data());
   container.Read(height.Data());
   container.Read(radius_x.Data());

   radius_y = UnitVec(height) ^ radius_x;
};

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
*/
void InitialSpaceCylinder::EvaluateInitial(void)
{
   double z, phi;

   z = rng->GetUniform();
   phi = twopi * rng->GetUniform();

   _pos = origin + radius_x * cos(phi) + radius_y * sin(phi) + height * z;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// InitialSpaceCylinderSector methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
*/
InitialSpaceCylinderSector::InitialSpaceCylinderSector(void)
                          : InitialSpaceCylinder(init_name_space_cylinder_sector, 0, INITIAL_SPACE | INITIAL_SURFACE)
{
};

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupInitial()" with the argument of "true".
*/
InitialSpaceCylinderSector::InitialSpaceCylinderSector(const InitialSpaceCylinderSector& other)
                          : InitialSpaceCylinder(other)
{
   RAISE_BITS(_status, INITIAL_SPACE);
   RAISE_BITS(_status, INITIAL_SURFACE);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupInitial(true);
};

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void InitialSpaceCylinderSector::SetupInitial(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) InitialSpaceCylinder::SetupInitial(false);

   container.Read(&phi1);
   // if(phi1 < 0.0) phi1 = 0.0;

   container.Read(&phi2);
   if(phi2 > twopi) phi2 = twopi;
   if(phi2 < phi1) phi2 = phi1;
};

/*!
\author Juan G Alonso Guzman
\date 05/16/2023
*/
void InitialSpaceCylinderSector::EvaluateInitial(void)
{
   double z, phi;

   z = rng->GetUniform();
   phi = phi1 + (phi2 - phi1) * rng->GetUniform();

   _pos = origin + radius_x * cos(phi) + radius_y * sin(phi) + height * z;
};

};
