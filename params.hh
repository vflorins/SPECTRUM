/*!
\file params.hh
\brief Declares a simple class for entering parameters
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_PARAMS_HH
#define SPECTRUM_PARAMS_HH

// This includes (algorithm, cmath, cstdint, cstring, fstream, vector), definitions, multi_index
#include "vectors.hh"
#include "data_container.hh"
#include "random.hh"
#include <exception>
#include <memory>

namespace Spectrum {

//! Zero state (for initialization)
const uint16_t STATE_NONE = 0x0000;

//! The internal state is invalid
const uint16_t STATE_INVALID = 0x0001;

//! Setup was completed
const uint16_t STATE_SETUP_COMPLETE = 0x0002;

//! The model has no time dependence
const uint16_t MODEL_STATIC = 0x0004;

//! The model is mesh based
const uint16_t MODEL_MESH_BASED = 0x0008;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Exceptions
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Exception if operation is performed on an uninitialized object
\author Vladimir Florinski
*/
class ExUninitialized : public std::exception {

public:

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Vladimir Florinski
\date 12/09/2021
\return Text describing the error
*/
inline const char* ExUninitialized::what(void) const noexcept
{
   return "Operation on an uninitialized object";
};

/*!
\brief Exception if coordinates are incompatible
\author Vladimir Florinski
*/
class ExCoordinates : public std::exception {

public:

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Vladimir Florinski
\date 12/09/2021
\return Text describing the error
*/
inline const char* ExCoordinates::what(void) const noexcept
{
   return "Incompatible coordinates";
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Params class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Base class for background, trajectory, boundary, distribution, and initial classes
\author Vladimir Florinski
*/
class Params {

protected:

//! Readable name of the class (persistent)
   std::string class_name = "";

//! Particle's specie (persistent)
   unsigned int specie = 0;

//! Random number generator object (persistent)
   std::shared_ptr<RNG> rng = nullptr;

//! Parameter storage (persistent)
   DataContainer container;

//! Status
   uint16_t _status = STATE_NONE;

//! Time (transient)
   double _t;

//! Spatial position (transient)
   GeoVector _pos;

//! Velocity vector (transient)
   GeoVector _vel;

//! Momentum vector (transient)
   GeoVector _mom;

//! Default constructor (protected, class not designed to be instantiated)
   Params(void) = default;

//! Constructor with arguments (to speed up construction of derived classes)
   Params(const std::string& name_in, unsigned int specie_in, uint16_t status_in);

//! Copy constructor (protected, class not designed to be instantiated)
   Params(const Params& other);

public:

//! Destructor
   ~Params() = default;

//! Return the name of the class
   std::string GetName(void) const;

//! Connect to an existing RNG object
   void ConnectRNG(const std::shared_ptr<RNG> rng_in);

//! Copy the user-supplied data container into "container"
   void SetContainer(const DataContainer& cont_in);

//! Return a copy of the data container
   DataContainer GetContainer(void) const;

//! Set the particle specie
   void SetSpecie(unsigned int specie_in);

//! Return the particle specie
   unsigned int GetSpecie(void) const;

//! Set the internal phase space position
   void SetState(double t_in, const GeoVector& pos_in, const GeoVector& mom_in = gv_zeros);

//! Return the internal phase space position
   void GetState(double& t_out, GeoVector& pos_out, GeoVector& mom_out) const;

//! Return the status
   uint16_t GetStatus(void) const;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Params inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 11/25/2020
*/
inline std::string Params::GetName(void) const
{
   return class_name;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 05/27/2022
\param[in] rng_in A pointer to an RNG object
*/
inline void Params::ConnectRNG(const std::shared_ptr<RNG> rng_in)
{
   rng = rng_in;
};

/*!
\author Vladimir Florinski
\date 12/03/2020
\param[in] cont_in Container with parameters
*/
inline void Params::SetContainer(const DataContainer& cont_in)
{
   container = cont_in;
};

/*!
\author Juan G Alonso Guzman
\date 04/28/2021
\return Params data container
*/
inline DataContainer Params::GetContainer(void) const
{
   return container;
};

/*!
\author Vladimir Florinski
\date 05/24/2022
\param[in] specie_in Index of the particle species defined in physics.hh
*/
inline void Params::SetSpecie(unsigned int specie_in)
{
   specie = specie_in;
};

/*!
\author Vladimir Florinski
\date 05/24/2022
\return Index of the particle species defined in physics.hh
*/
inline unsigned int Params::GetSpecie(void) const
{
   return specie;
};

/*!
\author Vladimir Florinski
\date 11/30/2020
\return Internal status
*/
inline uint16_t Params::GetStatus(void) const
{
   return _status;
};

};

#endif
