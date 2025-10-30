/*!
\file arithmetic.hh
\brief Defines common binary operators acting on SimpleArray and derived classes
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_ARITHMETIC_HH
#define SPECTRUM_ARITHMETIC_HH

#include <type_traits>
#include <complex>

#include "config.h"
#include "common/gpu_config.hh"

namespace Spectrum {

// Ref: https://gist.github.com/sighingnow/505d3d5c82237741b4a18147b2f84811
template <template <typename...> typename T, typename U>
struct is_specialization_of : std::false_type {};

// Ref: https://gist.github.com/sighingnow/505d3d5c82237741b4a18147b2f84811
template <template <typename...> typename T, typename... Us>
struct is_specialization_of<T, T<Us...>> : std::true_type {};

/*!
\brief A template tester for arithmetic type extending it to include std::complex
\author Vladimir Florinski
*/
template <typename T>
struct is_arithmetic_or_complex : std::disjunction<std::is_arithmetic<T>, is_specialization_of<std::complex, T>> {};

/*!
\brief Add two objects together
\author Vladimir Florinski
\date 04/25/2024
\param[in] left  Left operand \f$\mathbf{v}_1\f$
\param[in] right Right operand \f$\mathbf{v}_2\f$
\return \f$\mathbf{v}_1+\mathbf{v}_2\f$
*/
template <typename T, std::enable_if_t<T::is_simple_array, bool> = true>
SPECTRUM_DEVICE_FUNC inline T operator +(const T& left, const T& right)
{
   T retval(left);
   retval += right;
   return retval;
};

/*!
\brief Increment each component by the same amount (as left operand)
\author Vladimir Florinski
\date 04/25/2024
\param[in] left  Left operand \f$a\f$
\param[in] right Right operand \f$\mathbf{v}\f$
\return \f$a(1,..,1)+\mathbf{v}\f$
*/
template <typename T, typename arithm, std::enable_if_t<T::is_simple_array && is_arithmetic_or_complex<arithm>::value, bool> = true>
SPECTRUM_DEVICE_FUNC inline T operator +(const T& left, arithm right)
{
   T retval(left);
   retval += right;
   return retval;
};

/*!
\brief Increment each component by the same amount (as right operand)
\author Vladimir Florinski
\date 04/25/2024
\param[in] left  Left operand \f$\mathbf{v}\f$
\param[in] right Right operand \f$a\f$
\return \f$\mathbf{v}+a(1,..,1)\f$
*/
template <typename T, typename arithm, std::enable_if_t<T::is_simple_array && is_arithmetic_or_complex<arithm>::value, bool> = true>
SPECTRUM_DEVICE_FUNC inline T operator +(arithm left, const T& right)
{
   T retval(left);
   retval += right;
   return retval;
};

/*!
\brief Subtract one object from another
\author Vladimir Florinski
\date 04/25/2024
\param[in] left  Left operand \f$\mathbf{v}_1\f$
\param[in] right Right operand \f$\mathbf{v}_2\f$
\return \f$\mathbf{v}_1-\mathbf{v}_2\f$
*/
template <typename T, std::enable_if_t<T::is_simple_array, bool> = true>
SPECTRUM_DEVICE_FUNC inline T operator -(const T& left, const T& right)
{
   T retval(left);
   retval -= right;
   return retval;
};

/*!
\brief Invert object and increment each component by the same amount (as left operand)
\author Vladimir Florinski
\date 04/25/2024
\param[in] left  Left operand \f$a\f$
\param[in] right Right operand \f$\mathbf{v}\f$
\return \f$a(1,..,1)-\mathbf{v}\f$
*/
template <typename T, typename arithm, std::enable_if_t<T::is_simple_array && is_arithmetic_or_complex<arithm>::value, bool> = true>
SPECTRUM_DEVICE_FUNC inline T operator -(arithm left, const T& right)
{
   T retval(left);
   retval -= right;
   return retval;
};

/*!
\brief Decrement each component by the same amount (as right operand)
\author Vladimir Florinski
\date 04/25/2024
\param[in] left  Left operand \f$\mathbf{v}\f$
\param[in] right Right operand \f$a\f$
\return \f$\mathbf{v}-a(1,..,1)\f$
*/
template <typename T, typename arithm, std::enable_if_t<T::is_simple_array && is_arithmetic_or_complex<arithm>::value, bool> = true>
SPECTRUM_DEVICE_FUNC inline T operator -(const T& left, arithm right)
{
   T retval(left);
   retval -= right;
   return retval;
};

/*!
\brief Multiply each component by the same amount (as left operand)
\author Vladimir Florinski
\date 04/25/2024
\param[in] left  Left operand \f$a\f$
\param[in] right Right operand \f$\mathbf{v}\f$
\return \f$a\mathbf{v}\f$
*/
template <typename T, typename arithm, std::enable_if_t<T::is_simple_array && is_arithmetic_or_complex<arithm>::value, bool> = true>
SPECTRUM_DEVICE_FUNC inline T operator *(arithm left, const T& right)
{
   T retval(left);
   retval *= right;
   return retval;
};

/*!
\brief Multiply each component by the same amount (as right operand)
\author Vladimir Florinski
\date 04/25/2024
\param[in] left  Left operand \f$\mathbf{v}\f$
\param[in] right Right operand \f$a\f$
\return \f$a\mathbf{v}\f$
*/
template <typename T, typename arithm, std::enable_if_t<T::is_simple_array && is_arithmetic_or_complex<arithm>::value, bool> = true>
SPECTRUM_DEVICE_FUNC inline T operator *(const T& left, arithm right)
{
   T retval(left);
   retval *= right;
   return retval;
};

/*!
\brief Divide an object by another component-wise
\author Vladimir Florinski
\date 04/25/2024
\param[in] left  Left operand \f$\mathbf{v}_1\f$
\param[in] right Right operand \f$\mathbf{v}_2\f$
\return Result of the component-wise division
*/
template <typename T, std::enable_if_t<T::is_simple_array, bool> = true>
SPECTRUM_DEVICE_FUNC inline T operator /(const T& left, const T& right)
{
   T retval(left);
   retval /= right;
   return retval;
};

/*!
\brief Divide each component by the same amount
\author Vladimir Florinski
\date 04/25/2024
\param[in] left  Left operand \f$\mathbf{v}\f$
\param[in] right Right operand \f$a\f$
\return \f$a^{-1}\mathbf{v}\f$
*/
template <typename T, typename arithm, std::enable_if_t<T::is_simple_array && is_arithmetic_or_complex<arithm>::value, bool> = true>
SPECTRUM_DEVICE_FUNC inline T operator /(const T& left, arithm right)
{
   T retval(left);
   retval /= right;
   return retval;
};

/*!
\brief Modulo divide an object by another component-wise
\author Vladimir Florinski
\date 04/25/2024
\param[in] left  Left operand \f$\mathbf{v}_1\f$
\param[in] right Right operand \f$\mathbf{v}_2\f$
\return Result of the component-wise modulo division
*/
template <typename T, std::enable_if_t<T::is_simple_array, bool> = true>
SPECTRUM_DEVICE_FUNC inline T operator %(const T& left, const T& right)
{
   T retval(left);
   retval %= right;
   return retval;
};

/*!
\brief Modulo divide each component by the same amount
\author Vladimir Florinski
\date 04/25/2024
\param[in] left  Left operand \f$\mathbf{v}\f$
\param[in] right Right operand \f$a\f$
\return First object modulo divided by a number
*/
template <typename T, typename arithm, std::enable_if_t<T::is_simple_array && is_arithmetic_or_complex<arithm>::value, bool> = true>
SPECTRUM_DEVICE_FUNC inline T operator %(const T& left, arithm right)
{
   T retval(left);
   retval %= right;
   return retval;
};

};

#endif
