/*!
\file distance_map.cc
\brief Mapping between reference and physical distances
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific nThis file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifdef GEO_DEBUG
#include <iostream>
#include <iomanip>
#endif
#include "geometry/distance_map.hh"

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
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistance(true);
};

/*!
\author Vladimir Florinski
\date 05/15/2024
\param [in] construct Whether called from a copy constructor or separately
*/
void DistanceExponential::SetupDistance(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) DistanceBase::SetupDistance(false);

   ratio = rmax / rmin;
   log_ratio = log(ratio);
};

/*!
\author Vladimir Florinski
\date 05/15/2024
*/
void DistanceExponential::EvaluateDistance(void)
{
   if(selector == 0) _pos[0] = rmin * pow(ratio, _pos[1]);
   else if(selector == 1) _pos[1] = log(_pos[0] / rmin) / log_ratio;
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
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupDistance(true);
};

/*!
\author Vladimir Florinski
\date 05/14/2024
\param [in] construct Whether called from a copy constructor or separately
*/
void DistancePowerLaw::SetupDistance(bool construct)
{
// The parent version must be called explicitly if not constructing
   if(!construct) DistanceBase::SetupDistance(false);

   container.Read(&pow_law);
   pow_inv = 1.0 / pow_law;
   ref0 = 1.0 / (pow(rmax / rmin, pow_inv) - 1.0);
   c = pow_law * rmin / ref0;
};

/*!
\author Vladimir Florinski
\date 05/14/2024
*/
void DistancePowerLaw::EvaluateDistance(void)
{
   if(selector == 0) _pos[0] = rmin * pow(1.0 + _pos[1] / ref0, pow_law);
   else if(selector == 1) _pos[1] = ref0 * (pow(_pos[0] / rmin, pow_inv) - 1.0);
   else _pos[2] = c * pow(1.0 + _pos[1] / ref0, pow_law - 1.0);
};

#ifdef GEO_DEBUG

/*!
\author Vladimir Florinski
\date 05/14/2024
\param[in] pwl_in Powr law for the DistancePowerLaw class
*/
void PrintDistanceMaps(double pwl_in)
{
   double n = 60;
   double rmin = 1.0;
   double rmax = 100.0;
   double pwl = pwl_in;
   double xi, r, drdxi;
   DataContainer container;

   DistanceExponential exp_map;
   DistancePowerLaw pwl_map;

   container.Clear();
   container.Insert(rmin);
   container.Insert(rmax);
   exp_map.SetupObject(container);

   container.Insert(pwl);
   pwl_map.SetupObject(container);

   std::setprecision(6);
   for(auto i = 0; i <= n; i++) {

// Test the exponential map
      xi = i / (double)n;
      std::cerr << std::setw(16) << xi;
      r = exp_map.GetPhysical(xi);
      std::cerr << std::setw(16) << r;
      drdxi = exp_map.GetDerivative(xi);
      std::cerr << std::setw(16) << drdxi;
      xi = exp_map.GetReference(r);
      std::cerr << std::setw(16) << xi;

// Test the power law map
      xi = i / (double)n;
      r = pwl_map.GetPhysical(xi);
      std::cerr << std::setw(16) << r;
      drdxi = pwl_map.GetDerivative(xi);
      std::cerr << std::setw(16) << drdxi;
      xi = pwl_map.GetReference(r);
      std::cerr << std::setw(16) << xi;

      std::cerr << std::endl;
   };
};

#endif

};
