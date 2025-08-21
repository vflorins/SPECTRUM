/*!
\file random.hh
\brief Declares random number generator class
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_RANDOM_HH
#define SPECTRUM_RANDOM_HH

#include <config.h>

#ifdef USE_GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RNG class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief A simple wrapper for the GSL random number generator
\author Vladimir Florinski
*/
class RNG {

private:

//! An RNG from the GSL library
   gsl_rng* rng_internal = nullptr;

public:

//! Default constructor
   RNG(void);

//! Constructor with arguments
   RNG(int seed);

//! Destructor
   ~RNG(void);

//! Return uniformly distributed number on 0..1
   double GetUniform(void) const;

//! Return normally distributed number with variance 1
   double GetNormal(void) const;

//! Return a normally distributed radius in polar coordinates with variance 1
   double GetRayleigh(void) const;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// RNG inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/20/2022
*/
inline RNG::RNG(void)
{
   rng_internal = gsl_rng_alloc(gsl_rng_default);
};

/*!
\author Vladimir Florinski
\date 05/20/2022
\param[in] rng_in A pointer to a previously initialized and seeded GSL random number generator
*/
inline RNG::RNG(int seed)
{
   rng_internal = gsl_rng_alloc(gsl_rng_default);
   gsl_rng_set(rng_internal, seed);
};

/*!
\author Vladimir Florinski
\date 05/20/2022
*/
inline RNG::~RNG()
{
   if (rng_internal) gsl_rng_free(rng_internal);
};

/*!
\author Vladimir Florinski
\date 03/04/2022
\return A random number uniformly distributed between 0 and 1
*/
inline double RNG::GetUniform(void) const
{
   return gsl_rng_uniform(rng_internal);
};

/*!
\author Vladimir Florinski
\date 03/04/2022
\return A random number normally distributed with mean 0 and variance 1
*/
inline double RNG::GetNormal(void) const
{
   return gsl_ran_ugaussian(rng_internal);
};

/*!
\author Vladimir Florinski
\date 04/03/2023
\return A random number distributed as a radial coordinate on a disk with variance 1
*/
inline double RNG::GetRayleigh(void) const
{
   return  gsl_ran_rayleigh(rng_internal, 1.0);
};

};

#endif
