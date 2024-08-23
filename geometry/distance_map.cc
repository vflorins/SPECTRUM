/*!
\file distance_map.cc
\brief Mapping between reference and physical distances
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific nThis file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "geometry/distance_map.hh"

#ifdef GEO_DEBUG
#include <iostream>
#include <iomanip>
#endif

#include <gsl/gsl_sf_lambert.h>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistanceBase methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/14/2024
*/
DistanceBase::DistanceBase(void)
            : Params("", 0, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 05/14/2024
\param[in] name_in   Readable name of the class
\param[in] specie_in Particle's specie
\param[in] status_in Initial status
*/
DistanceBase::DistanceBase(const std::string& name_in, unsigned int specie_in, uint16_t status_in)
            : Params(name_in, specie_in, status_in)
{
};

/*!
\brief Copy constructor (protected, class not designed to be instantiated)
\author Vladimir Florinski
\date 05/14/2024
\param[in] other Object to initialize from
*/
DistanceBase::DistanceBase(const DistanceBase& other)
            : Params(other)
{
// Params' constructor sets the state to "STATE_NONE"
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistance(true);
};

/*!
\author Vladimir Florinski
\date 05/14/2024
\param[in] cont_in Container with parameters
*/
void DistanceBase::SetupObject(const DataContainer& cont_in)
{
   Params::SetContainer(cont_in);
   SetupDistance(false);
};

/*!
\author Vladimir Florinski
\date 05/14/2024
\param [in] construct Whether called from a copy constructor or separately
*/
void DistanceBase::SetupDistance(bool construct)
{
// Only needed in the parent version
   container.Reset();
   container.Read(&rmin);
   container.Read(&rmax);
   
   RAISE_BITS(_status, STATE_SETUP_COMPLETE);
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\date 05/14/2024
\note This is only a stub
*/
void DistanceBase::EvaluateDistance(void)
{
   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\date 05/15/2024
\param[in] ref Reference variable
\return Physical variable
*/
double DistanceBase::GetPhysical(double ref)
{
   _pos[1] = ref;
   selector = 0;
   EvaluateDistance();
   return _pos[0];
};

/*!
\author Vladimir Florinski
\date 05/15/2024
\param[in] phys Physical variable
\return Reference variable
*/
double DistanceBase::GetReference(double phys)
{
   _pos[0] = phys;
   selector = 1;
   EvaluateDistance();
   return _pos[1];
};

/*!
\author Vladimir Florinski
\date 05/15/2024
\param[in] phys Reference variable
\return Derivative of the physical variable
*/
double DistanceBase::GetDerivative(double ref)
{
   _pos[1] = ref;
   selector = 2;
   EvaluateDistance();
   return _pos[2];
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistanceExponential methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/15/2024
*/
DistanceExponential::DistanceExponential(void)
                   : DistanceBase(distance_name_exponential, 0, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 05/15/2024
\param[in] other Object to initialize from
*/
DistanceExponential::DistanceExponential(const DistanceExponential& other)
                   : DistanceBase(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistance(true);
};

/*!
\author Vladimir Florinski
\date 05/15/2024
\param [in] construct Whether called from a copy constructor or separately
*/
void DistanceExponential::SetupDistance(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistanceBase::SetupDistance(false);

   ratio = rmax / rmin;
   log_ratio = log(ratio);
};

/*!
\author Vladimir Florinski
\date 05/15/2024
*/
void DistanceExponential::EvaluateDistance(void)
{
   if (selector == 0) _pos[0] = rmin * pow(ratio, _pos[1]);
   else if (selector == 1) _pos[1] = log(_pos[0] / rmin) / log_ratio;
   else _pos[2] = rmin * log_ratio * pow(ratio, _pos[1]);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistancePowerLaw methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/14/2024
*/
DistancePowerLaw::DistancePowerLaw(void)
                : DistanceBase(distance_name_power_law, 0, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 05/14/2024
\param[in] other Object to initialize from
*/
DistancePowerLaw::DistancePowerLaw(const DistancePowerLaw& other)
                : DistanceBase(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistance(true);
};

/*!
\author Vladimir Florinski
\date 05/14/2024
\param [in] construct Whether called from a copy constructor or separately
*/
void DistancePowerLaw::SetupDistance(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistanceBase::SetupDistance(false);

   container.Read(&powl);
   pow_inv = 1.0 / powl;
   ref0 = 1.0 / (pow(rmax / rmin, pow_inv) - 1.0);
   c = powl * rmin / ref0;
};

/*!
\author Vladimir Florinski
\date 05/14/2024
*/
void DistancePowerLaw::EvaluateDistance(void)
{
   if (selector == 0) _pos[0] = rmin * pow(1.0 + _pos[1] / ref0, powl);
   else if (selector == 1) _pos[1] = ref0 * (pow(_pos[0] / rmin, pow_inv) - 1.0);
   else _pos[2] = c * pow(1.0 + _pos[1] / ref0, powl - 1.0);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// DistanceLinearExp methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 08/22/2024
*/
DistanceLinearExp::DistanceLinearExp(void)
                 : DistanceBase(distance_name_linear_exp, 0, STATE_NONE)
{
};

/*!
\author Vladimir Florinski
\date 08/22/2024
\param[in] other Object to initialize from
*/
DistanceLinearExp::DistanceLinearExp(const DistanceLinearExp& other)
                 : DistanceBase(other)
{
   if (BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistance(true);
};

/*!
\author Vladimir Florinski
\date 08/22/2024
\param [in] construct Whether called from a copy constructor or separately
*/
void DistanceLinearExp::SetupDistance(bool construct)
{
// The parent version must be called explicitly if not constructing
   if (!construct) DistanceBase::SetupDistance(false);

   double frac, rmed;

   container.Read(&rmed);
   container.Read(&frac);
   container.Read(&logy);
   
   base = exp(logy);
   C = (rmed / rmin - 1.0 - frac * (rmax / rmin - 1.0)) / (pow(base, frac) - 1.0 - frac * (base - 1.0));
   ref0 = 1.0 / (rmax / rmin - 1.0 - C * (base - 1.0));
   c = C * rmin * logy;
};

/*!
\author Vladimir Florinski
\date 08/22/2024
*/
void DistanceLinearExp::EvaluateDistance(void)
{
   double chi;
   if (selector == 0) _pos[0] = rmin * (1.0 - C + _pos[1] / ref0 + C * pow(base, _pos[1]));
   else if (selector == 1) {
      chi = ref0 * (_pos[0] / rmin - 1.0 + C);
      _pos[1] = chi - gsl_sf_lambert_W0(ref0 * C * logy * pow(base, chi)) / logy;
   }
   else _pos[2] = rmin / ref0 + c * pow(base, _pos[1]);
};

#ifdef GEO_DEBUG

/*!
\author Vladimir Florinski
\date 08/23/2024
\param [in] dist_map An object derived from "DistanceBase"
\param [in] as_phys  Print the physical distance as first column
*/
void PrintDistanceMaps(DistanceBase* dist_map, bool as_phys)
{
   const double n = 640;
   const double rmin = 30.0;
   const double rmed = 250.0;
   const double rmax = 1000.0;
   const double powl = 2.0;
   const double frac = 0.4;
   const double logy = 10.0;

   DataContainer container;
   container.Clear();
   container.Insert(rmin);
   container.Insert(rmax);
   if (dist_map->GetName() == distance_name_power_law) {
      container.Insert(powl);
   }
   else if (dist_map->GetName() == distance_name_linear_exp) {
      container.Insert(rmed);
      container.Insert(frac);
      container.Insert(logy);
   };
   dist_map->SetupObject(container);

   double xi, r, dxi = 1.0 / n;
   std::cout << std::setprecision(6);

   for (auto i = 0; i <= n; i++) {
      xi = i * dxi;

// Check the forward map. The first column to print is either "xi" or "r" depending on what the user wants.
      r = dist_map->GetPhysical(xi);
      std::cout << std::setw(16) << (as_phys ? r : xi);
      std::cout << std::setw(16) << (as_phys ? xi : r);

// Print the mesh spacing
      std::cout << std::setw(16) << dist_map->GetDerivative(xi) * dxi;

// Check if the inverse map works
      xi = dist_map->GetReference(r);
      std::cout << std::setw(16) << xi;
      std::cout << std::endl;
   };
};

#endif

};
