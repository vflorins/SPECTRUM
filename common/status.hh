/*!
\file status.hh
\brief Define status code flags
\author Vladimir Florinski
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/


#ifndef SPECTRUM_STATUS_HH
#define SPECTRUM_STATUS_HH

#include <cstdint>
#include <exception>

namespace Spectrum {

using status_t = uint16_t;

//! Zero state (for initialization)
const status_t STATE_NONE = 0x0000;

//! The internal state is invalid
const status_t STATE_INVALID = 0x0001;

//! Setup was completed
const status_t STATE_SETUP_COMPLETE = 0x0002;

//! The model has no time dependence
const status_t MODEL_STATIC = 0x0004;

//! The model is mesh based
const status_t MODEL_MESH_BASED = 0x0008;

//! The trajectory will end after the step is completed
constexpr status_t TRAJ_FINISH = 0x0010;

//! Time boundary was crossed
constexpr status_t TRAJ_TIME_CROSSED = 0x0020;

//! Spatial boundary was crossed
constexpr status_t TRAJ_SPATIAL_CROSSED = 0x0040;

//! Momentum boundary was crossed
constexpr status_t TRAJ_MOMENTUM_CROSSED = 0x0080;

//! Trajectory is invalid and must be discarded
constexpr status_t TRAJ_DISCARD = 0x0100;

//! Background requires setup but setup has not occured
constexpr status_t BACKGROUND_SETUP_AWAIT = 0x0400;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Exceptions
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Exception if field evaluation failed
\author Vladimir Florinski
*/
struct ExFieldError : public std::exception {

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Vladimir Florinski
\date 12/09/2021
\return Text describing the error
*/
inline const char* ExFieldError::what(void) const noexcept
{
   return "Field evaluation error";
};

/*!
\brief Exception if server functions failed
\author Juan G Alonso Guzman
*/
class ExServerError : public std::exception {

public:

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Juan G Alonso Guzman
\date 23/08/2023
\return Text describing the error
*/
inline const char* ExServerError::what(void) const noexcept
{
   return "Server error";
};


//----------------------------------------------------------------------------------------------------------------------------------------------------
// Exceptions
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Exception if maximum number of steps is reached in trajectory
\author Juan G Alonso Guzman
*/
class ExMaxStepsReached : public std::exception {

public:

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2022
\return Text describing the error
*/
inline const char* ExMaxStepsReached::what(void) const noexcept
{
   return "Maximum number of steps in a trajectory reached.";
};

/*!
\brief Exception if maximum number of time step adaptations in a single step is reached
\author Juan G Alonso Guzman
*/
class ExMaxTimeAdaptsReached : public std::exception {

public:

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2022
\return Text describing the error
*/
inline const char* ExMaxTimeAdaptsReached::what(void) const noexcept
{
   return "Maximum number of time adaptations in a single step reached.";
};

/*!
\brief Exception if time step becomes too small
\author Juan G Alonso Guzman
*/
class ExTimeStepTooSmall : public std::exception {

public:

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Juan G Alonso Guzman
\date 07/14/2022
\return Text describing the error
*/
inline const char* ExTimeStepTooSmall::what(void) const noexcept
{
   return "Time step became too small.";
};

/*!
\brief Exception if time step becomes nan
\author Juan G Alonso Guzman
*/
class ExTimeStepNan : public std::exception {

public:

//! Return explanatory string
   const char* what(void) const noexcept override;
};

/*!
\author Juan G Alonso Guzman
\date 07/20/2022
\return Text describing the error
*/
inline const char* ExTimeStepNan::what(void) const noexcept
{
   return "Time step became nan.";
};



};

#endif
