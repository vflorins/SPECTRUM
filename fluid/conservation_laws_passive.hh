/*!
\file conservation_laws_passive.hh
\brief Declares rules for a passively advected scalar (PAS)
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_CONSERVATION_LAWS_PASSIVE_HH
#define SPECTRUM_CONSERVATION_LAWS_PASSIVE_HH

#include "common/vectors.hh"

namespace Spectrum {

//! Number of fluid variables
#define CL_PASSIVE_NVARS 1

/*!
\brief Codifies physics of passively advected scalar
\author Vladimir Florinski
*/
class Passive
{
public:

   struct ConservedState;
   struct FluxFunction;

//! A trait to be used in template specializations
   static constexpr bool is_cons_law_active = false;

//! Return the number of variables
   SPECTRUM_DEVICE_FUNC static constexpr int Nvars(void) {return CL_PASSIVE_NVARS;};

/*!
\brief Codifies physics of passively advected scalar for primitive variables
\author Vladimir Florinski
*/
   struct PrimitiveState : public SimpleArray<double, CL_PASSIVE_NVARS>
   {
      using SimpleArray<double, CL_PASSIVE_NVARS>::data;

   //! Alias for PAS (RO)
      const double& pas(void) const {return data[0];};

   //! Alias for PAS (RW)
      double& pas(void) {return data[0];};

   //! Default constructor
      SPECTRUM_DEVICE_FUNC PrimitiveState(void) = default;

   //! Constructor from components
      SPECTRUM_DEVICE_FUNC PrimitiveState(double pas_in);

   //! Constructor from SimpleArray
      SPECTRUM_DEVICE_FUNC PrimitiveState(const SimpleArray<double, CL_PASSIVE_NVARS>& other);

   //! Calculate conserved state
      template <typename ActiveConsLaw, std::enable_if_t<ActiveConsLaw::is_cons_law_active, bool> = true>
      SPECTRUM_DEVICE_FUNC ConservedState ToConserved(const ActiveConsLaw::PrimitiveState& active_prim) const;

   //! Calculate flux function
      template <typename ActiveConsLaw, std::enable_if_t<ActiveConsLaw::is_cons_law_active, bool> = true>
      SPECTRUM_DEVICE_FUNC FluxFunction ToFlux(const ActiveConsLaw::PrimitiveState& active_prim) const;
   };

/*!
\brief Codifies physics of passively advected scalar for conserved variables
\author Vladimir Florinski
*/
   struct ConservedState : public SimpleArray<double, CL_PASSIVE_NVARS>
   {
      using SimpleArray<double, CL_PASSIVE_NVARS>::data;

   //! Alias for primary density*PAS (RO)
      const double& dpas(void) const {return data[0];};

   //! Alias for primary density*PAS (RW)
      double& dpas(void) {return data[0];};

   //! Default constructor
      SPECTRUM_DEVICE_FUNC ConservedState(void) = default;

   //! Constructor from components
      SPECTRUM_DEVICE_FUNC ConservedState(double dpas_in);

   //! Constructor from SimpleArray
      SPECTRUM_DEVICE_FUNC ConservedState(const SimpleArray<double, CL_PASSIVE_NVARS>& other);

   //! Calculate primitive state
      template <typename ActiveConsLaw, std::enable_if_t<ActiveConsLaw::is_cons_law_active, bool> = true>
      SPECTRUM_DEVICE_FUNC PrimitiveState ToPrimitive(const ActiveConsLaw::PrimitiveState& active_prim) const;

   //! Calculate flux function
      template <typename ActiveConsLaw, std::enable_if_t<ActiveConsLaw::is_cons_law_active, bool> = true>
      SPECTRUM_DEVICE_FUNC FluxFunction ToFlux(const ActiveConsLaw::PrimitiveState& active_prim) const;
   };

/*!
\brief Codifies physics of passively advected scalar for the fluxes
\author Vladimir Florinski
*/
   struct FluxFunction : public SimpleArray<double, CL_PASSIVE_NVARS>
   {
      using SimpleArray<double, CL_PASSIVE_NVARS>::data;

   //! Alias for primary density*PAS flux (RO)
      const double& dpasf(void) const {return data[0];};

   //! Alias for primary density*PAS flux (RW)
      double& dpasf(void) {return data[0];};

   //! Default constructor
      SPECTRUM_DEVICE_FUNC FluxFunction(void) = default;

   //! Constructor from components
      SPECTRUM_DEVICE_FUNC FluxFunction(double dpasf_in);

   //! Constructor from SimpleArray
      SPECTRUM_DEVICE_FUNC FluxFunction(const SimpleArray<double, CL_PASSIVE_NVARS>& other);
   };

};

};

#endif
