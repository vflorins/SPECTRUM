/*!
\file conservation_laws_passive.hh
\brief Declares rules for a passively advected scalar (PAS)
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_CONSERVATION_LAWS_PASSIVE_HH
#define SPECTRUM_CONSERVATION_LAWS_PASSIVE_HH

namespace Spectrum {

//! Number of fluid variables
#define CL_PASSIVE_NVARS 1

struct PrimitiveStatePassive;
struct ConservedStatePassive;
struct FluxFunctionPassive;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// PrimitiveStatePassive class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Codifies physics for a passively advected scalar variable
\author Vladimir Florinski
*/
struct PrimitiveStatePassive : public SimpleArray<double, CL_PASSIVE_NVARS>
{
   using SimpleArray<double, CL_PASSIVE_NVARS>::data;

//! A trait to be used in template specializations
   static constexpr bool is_cons_law_active = false;

//! Alias for PAS (RO)
   const double& pas(void) const {return data[0];};

//! Alias for PAS (RW)
   double& pas(void) {return data[0];};

//! Return the number of variables
   SPECTRUM_DEVICE_FUNC static constexpr int Nvars(void) {return CL_PASSIVE_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC PrimitiveStatePassive(void) = default;

//! Constructor from components
   SPECTRUM_DEVICE_FUNC PrimitiveStatePassive(double pas_in);

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC PrimitiveStatePassive(const SimpleArray<double, CL_PASSIVE_NVARS>& other);

//! Calculate conserved state
   template <typename cl_prim, std::enable_if_t<cl_prim::is_cons_law_active, bool> = true>
   SPECTRUM_DEVICE_FUNC ConservedStatePassive ToConserved(cl_prim active_prim) const;

//! Calculate flux function
   template <typename cl_prim, std::enable_if_t<cl_prim::is_cons_law_active, bool> = true>
   SPECTRUM_DEVICE_FUNC FluxFunctionPassive ToFlux(cl_prim active_prim) const;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ConservedStatePassive class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

struct ConservedStatePassive : public SimpleArray<double, CL_PASSIVE_NVARS>
{
   using SimpleArray<double, CL_PASSIVE_NVARS>::data;

//! A trait to be used in template specializations
   static constexpr bool is_cons_law_active = false;

//! Alias for primary density*PAS (RO)
   const double& dpas(void) const {return data[0];};

//! Alias for primary density*PAS (RW)
   double& dpas(void) {return data[0];};

//! Return the number of variables
   SPECTRUM_DEVICE_FUNC static constexpr int Nvars(void) {return CL_PASSIVE_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC ConservedStatePassive(void) = default;

//! Constructor from components
   SPECTRUM_DEVICE_FUNC ConservedStatePassive(double dpas_in);

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC ConservedStatePassive(const SimpleArray<double, CL_PASSIVE_NVARS>& other);

//! Calculate primitive state
   template <typename cl_prim, std::enable_if_t<cl_prim::is_cons_law_active, bool> = true>
   SPECTRUM_DEVICE_FUNC PrimitiveStatePassive ToPrimitive(cl_prim active_prim) const;

//! Calculate flux function
   template <typename cl_prim, std::enable_if_t<cl_prim::is_cons_law_active, bool> = true>
   SPECTRUM_DEVICE_FUNC FluxFunctionPassive ToFlux(cl_prim active_prim) const;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// FluxFunctionPassive class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

struct FluxFunctionPassive : public SimpleArray<double, CL_PASSIVE_NVARS>
{
   using SimpleArray<double, CL_PASSIVE_NVARS>::data;

//! A trait to be used in template specializations
   static constexpr bool is_cons_law_active = false;

//! Alias for primary density*PAS flux (RO)
   const double& dpasf(void) const {return data[0];};

//! Alias for primary density*PAS flux (RW)
   double& dpasf(void) {return data[0];};

//! Return the number of variables
   SPECTRUM_DEVICE_FUNC static constexpr int Nvars(void) {return CL_PASSIVE_NVARS;};

//! Default constructor
   SPECTRUM_DEVICE_FUNC FluxFunctionPassive(void) = default;

//! Constructor from components
   SPECTRUM_DEVICE_FUNC FluxFunctionPassive(double dpasf_in);

//! Constructor from the base class
   SPECTRUM_DEVICE_FUNC FluxFunctionPassive(const SimpleArray<double, CL_PASSIVE_NVARS>& other);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// PrimitiveStatePassive inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/31/2024
\param[in] pas_in PAS
*/
SPECTRUM_DEVICE_FUNC inline PrimitiveStatePassive::PrimitiveStatePassive(double pas_in)
{
   pas() = pas_in;
};

/*!
\author Vladimir Florinski
\date 05/31/2024
\param[in] other Object to initialize from
*/
SPECTRUM_DEVICE_FUNC inline PrimitiveStatePassive::PrimitiveStatePassive(const SimpleArray<double, CL_PASSIVE_NVARS>& other)
                                                 : SimpleArray<double, CL_PASSIVE_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 05/31/2024
\return Vector of conserved variables
*/
template <typename cl_prim, std::enable_if_t<cl_prim::is_cons_law_active, bool>>
SPECTRUM_DEVICE_FUNC inline ConservedStatePassive PrimitiveStatePassive::ToConserved(cl_prim active_prim) const
{
   return ConservedStatePassive(active_prim.den() * pas());
};

/*!
\author Vladimir Florinski
\date 05/31/2024
\return Flux vector
*/
template <typename cl_prim, std::enable_if_t<cl_prim::is_cons_law_active, bool>>
SPECTRUM_DEVICE_FUNC inline FluxFunctionPassive PrimitiveStatePassive::ToFlux(cl_prim active_prim) const
{
   return FluxFunctionPassive(active_prim.den() * pas() * active_prim.vel()[0]);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// ConservedStatePassive inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/31/2024
\param[in] dpas_in primary density*PAS
*/
SPECTRUM_DEVICE_FUNC inline ConservedStatePassive::ConservedStatePassive(double dpas_in)
{
   dpas() = dpas_in;
};

/*!
\author Vladimir Florinski
\date 05/31/2024
\param[in] other Object to initialize from
*/
SPECTRUM_DEVICE_FUNC inline ConservedStatePassive::ConservedStatePassive(const SimpleArray<double, CL_PASSIVE_NVARS>& other)
                                                 : SimpleArray<double, CL_PASSIVE_NVARS>(other)
{
};

/*!
\author Vladimir Florinski
\date 05/31/2024
\return Vector of primitive variables
*/
template <typename cl_prim, std::enable_if_t<cl_prim::is_cons_law_active, bool>>
SPECTRUM_DEVICE_FUNC inline PrimitiveStatePassive ConservedStatePassive::ToPrimitive(cl_prim active_prim) const
{
   return PrimitiveStatePassive(dpas() / active_prim.den());
};

/*!
\author Vladimir Florinski
\date 05/31/2024
\return Flux vector
*/
template <typename cl_prim, std::enable_if_t<cl_prim::is_cons_law_active, bool>>
SPECTRUM_DEVICE_FUNC inline FluxFunctionPassive ConservedStatePassive::ToFlux(cl_prim active_prim) const
{
   return FluxFunctionPassive(dpas() * active_prim.vel()[0]);
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// FluxFunctionPassive inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/31/2024
\param[in] dvarf_in Flux of primary density*PAS
*/
SPECTRUM_DEVICE_FUNC inline FluxFunctionPassive::FluxFunctionPassive(double dpasf_in)
{
   dpasf() = dpasf_in;
};

/*!
\author Vladimir Florinski
\date 05/31/2024
\param[in] other Object to initialize from
*/
SPECTRUM_DEVICE_FUNC inline FluxFunctionPassive::FluxFunctionPassive(const SimpleArray<double, CL_PASSIVE_NVARS>& other)
                                               : SimpleArray<double, CL_PASSIVE_NVARS>(other)
{
};

};

#endif
