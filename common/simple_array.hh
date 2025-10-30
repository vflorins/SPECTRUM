/*!
\file simple_array.hh
\brief Declares and defines a base template for a contiguous storage class
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SIMPLE_ARRAY_HH
#define SPECTRUM_SIMPLE_ARRAY_HH

#include <cstring>
#include <ostream>

#include "common/definitions.hh"
#include "common/arithmetic.hh"

namespace Spectrum {

/*!
\brief An embedded type containing different aliases for the data
\author Vladimir Florinski
*/
template <typename data_type, int n_vars, std::enable_if_t<(n_vars > 0), bool> = true>
struct SimpleArrayBase
{
   union {
      data_type data[n_vars];
      struct {
         data_type x, y, z;
      };
      data_type ijk[n_vars];
      struct {
         data_type i, j, k;
      };
   };

//! Default constructor, required because gcc makes it deleted
   SPECTRUM_DEVICE_FUNC constexpr SimpleArrayBase(void) {};
};

// FIXME: nvcc bug prevents the following from compiling. For now, don't use "n_vars" of 1 or 2 in CUDA code.
#ifndef __CUDACC__

/*!
\brief Specialization of "SimpleArrayBase" for a one-component array
\author Vladimir Florinski
*/
template <typename data_type>
struct SimpleArrayBase<data_type, 1>
{
   union {
      data_type data[1];
      struct {
         data_type x;
      };
      data_type ijk[1];
      struct {
         data_type i;
      };
   };

//! Default constructor, required because gcc makes it deleted
   SPECTRUM_DEVICE_FUNC constexpr SimpleArrayBase(void) {};
};

/*!
\brief Specialization of "SimpleArrayBase" for a two-component array
\author Vladimir Florinski
*/
template <typename data_type>
struct SimpleArrayBase<data_type, 2>
{
   union {
      data_type data[2];
      struct {
         data_type x, y;
      };
      data_type ijk[2];
      struct {
         data_type i, j;
      };
   };

//! Default constructor, required because gcc makes it deleted
   SPECTRUM_DEVICE_FUNC constexpr SimpleArrayBase(void) {};
};

#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimpleArray class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief An array of arbitrary type and fixed size with arithmetic operations
\author Vladimir Florinski
\author Juan G Alonso Guzman
\note Copy constructor and "operator=" are compiler provided (same for derived classes), so they are not implemented here. Binary operators are external to the class. Derived classes must provide a conversion constructor from the base class to use those.
*/
template <typename data_type, int n_vars>
struct SimpleArray : SimpleArrayBase<data_type, n_vars>
{
//! A trait to be used in template specializations
   static constexpr bool is_simple_array = true;
   
   using SimpleArrayBase<data_type, n_vars>::data;

//! Default constructor
   SPECTRUM_DEVICE_FUNC constexpr SimpleArray(void) {};

//! Constructor from a single value
   SPECTRUM_DEVICE_FUNC explicit constexpr SimpleArray(data_type val) {operator =(val);};

//! Constructor from an array
   SPECTRUM_DEVICE_FUNC explicit constexpr SimpleArray(const data_type* other);

//! Return the number of components
   SPECTRUM_DEVICE_FUNC static constexpr int size(void) {return n_vars;};

//! Access to the data for reading
   SPECTRUM_DEVICE_FUNC const data_type* Data(void) const {return data;};

//! Access to the data for writing
   SPECTRUM_DEVICE_FUNC data_type* Data(void) {return data;};

//! Access to components for reading
   SPECTRUM_DEVICE_FUNC const data_type& operator [](int i) const {return data[i];};

//! Access to components for writing
   SPECTRUM_DEVICE_FUNC data_type& operator [](int i) {return data[i];};

//! Store the content of the simple array into an array
   SPECTRUM_DEVICE_FUNC void Store(data_type* other) const;

//! Set all components to the same value
   SPECTRUM_DEVICE_FUNC constexpr SimpleArray& operator =(data_type a);

//! Assignment operator from an array
   SPECTRUM_DEVICE_FUNC SimpleArray& operator =(const data_type* other);

//! Add a number to each components
   SPECTRUM_DEVICE_FUNC SimpleArray& operator +=(data_type a);

//! Add another simple array to this
   SPECTRUM_DEVICE_FUNC SimpleArray& operator +=(const SimpleArray& other);

//! Subtract a number from each components
   SPECTRUM_DEVICE_FUNC SimpleArray& operator -=(data_type a);

//! Subtract another simple array from this
   SPECTRUM_DEVICE_FUNC SimpleArray& operator -=(const SimpleArray& other);

//! Multiply each component by a number
   SPECTRUM_DEVICE_FUNC SimpleArray& operator *=(data_type a);

//! Compute a result of component-wise multiplication of two simple arrays
   SPECTRUM_DEVICE_FUNC SimpleArray& operator *=(const SimpleArray& other);

//! Divide each component by a number
   SPECTRUM_DEVICE_FUNC SimpleArray& operator /=(data_type a);

//! Compute a result of component-wise division of two simple arrays
   SPECTRUM_DEVICE_FUNC SimpleArray& operator /=(const SimpleArray& other);

//! Modulo divide each component by a number
   SPECTRUM_DEVICE_FUNC SimpleArray& operator %=(data_type a);

//! Compute a result of component-wise modulo division of two simple arrays
   SPECTRUM_DEVICE_FUNC SimpleArray& operator %=(const SimpleArray& other);

//! Computes the square of the norm of this simple array
   SPECTRUM_DEVICE_FUNC data_type Norm2(void) const;

//! Sum of all components
   SPECTRUM_DEVICE_FUNC data_type Sum(void) const;

//! Product of all components
   SPECTRUM_DEVICE_FUNC data_type Prod(void) const;

//! Smallest component
   SPECTRUM_DEVICE_FUNC data_type Smallest(void) const;

//! Largest component
   SPECTRUM_DEVICE_FUNC data_type Largest(void) const;

//! Negate all components
   SPECTRUM_DEVICE_FUNC void Negate(void);

//! Compute a scalar product with another simple array
   SPECTRUM_DEVICE_FUNC data_type ScalarProd(const SimpleArray& other) const;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// SimpleArray inline methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 03/08/2024
\param[in] other Array to initialize from
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline constexpr SimpleArray<data_type, n_vars>::SimpleArray(const data_type* other)
{
   memcpy(data, other, n_vars * sizeof(data_type));
};

/*!
\author Vladimir Florinski
\date 03/07/2024
\param[out] other Array to store the simple array in
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline void SimpleArray<data_type, n_vars>::Store(data_type* other) const
{
   memcpy(other, data, n_vars * sizeof(data_type));
};

/*!
\author Vladimir Florinski
\date 03/08/2024
\param[in] val Value to be assigned to each component
\return Reference to this object
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline constexpr SimpleArray<data_type, n_vars>& SimpleArray<data_type, n_vars>::operator =(data_type val)
{
   for (auto i = 0; i < n_vars; i++) data[i] = val;
   return *this;
};

/*!
\author Vladimir Florinski
\date 03/07/2024
\param[in] other Array to copy into this simple array
\return Reference to this object
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline SimpleArray<data_type, n_vars>& SimpleArray<data_type, n_vars>::operator =(const data_type* other)
{
   memcpy(data, other, n_vars * sizeof(data_type));
   return *this;
};

/*!
\author Vladimir Florinski
\date 03/12/2024
\param[in] a Right operand \fa\f$
\return \f$\mathbf{v}+a(1,..,1)\f$
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline SimpleArray<data_type, n_vars>& SimpleArray<data_type, n_vars>::operator +=(data_type a)
{
   for (auto i = 0; i < n_vars; i++) data[i] += a;
   return *this;
};

/*!
\author Vladimir Florinski
\date 03/08/2024
\param[in] other Right operand \f$\mathbf{v}_1\f$
\return \f$\mathbf{v}+\mathbf{v}_1\f$
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline SimpleArray<data_type, n_vars>& SimpleArray<data_type, n_vars>::operator +=(const SimpleArray<data_type, n_vars>& other)
{
   for (auto i = 0; i < n_vars; i++) data[i] += other.data[i];
   return *this;
};

/*!
\author Vladimir Florinski
\date 03/12/2024
\param[in] a Right operand \fa\f$
\return \f$\mathbf{v}-a(1,..,1)\f$
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline SimpleArray<data_type, n_vars>& SimpleArray<data_type, n_vars>::operator -=(data_type a)
{
   for (auto i = 0; i < n_vars; i++) data[i] -= a;
   return *this;
};

/*!
\author Vladimir Florinski
\date 03/13/2024
\param[in] other Right operand \f$\mathbf{v}_1\f$
\return \f$\mathbf{v}-\mathbf{v}_1\f$
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline SimpleArray<data_type, n_vars>& SimpleArray<data_type, n_vars>::operator -=(const SimpleArray<data_type, n_vars>& other)
{
   for (auto i = 0; i < n_vars; i++) data[i] -= other.data[i];
   return *this;
};

/*!
\author Vladimir Florinski
\date 03/13/2024
\param[in] a Right operand \f$a\f$
\return \f$a\mathbf{v}\f$
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline SimpleArray<data_type, n_vars>& SimpleArray<data_type, n_vars>::operator *=(data_type a)
{
   for (auto i = 0; i < n_vars; i++) data[i] *= a;
   return *this;
};

/*!
\author Vladimir Florinski
\date 03/13/2024
\param[in] other Right operand \f$\mathbf{v}_1\f$
\return This simple array scaled by the other simple array
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline SimpleArray<data_type, n_vars>& SimpleArray<data_type, n_vars>::operator *=(const SimpleArray<data_type, n_vars>& other)
{
   for (auto i = 0; i < n_vars; i++) data[i] *= other.data[i];
   return *this;
};

/*!
\author Vladimir Florinski
\date 03/13/2024
\param[in] a Right operand \f$a\f$
\return \f$a^{-1}\mathbf{v}\f$
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline SimpleArray<data_type, n_vars>& SimpleArray<data_type, n_vars>::operator /=(data_type a)
{
   for (auto i = 0; i < n_vars; i++) data[i] /= a;
   return *this;
};

/*!
\author Vladimir Florinski
\date 03/13/2024
\param[in] other Right operand \f$\mathbf{v}_1\f$
\return This simple array scaled by the other simple array
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline SimpleArray<data_type, n_vars>& SimpleArray<data_type, n_vars>::operator /=(const SimpleArray& other)
{
   for (auto i = 0; i < n_vars; i++) data[i] /= other.data[i];
   return *this;
};

/*!
\author Vladimir Florinski
\date 03/13/2024
\param[in] a Right operand \f$a\f$
\return This simple array modulo divided by a number
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline SimpleArray<data_type, n_vars>& SimpleArray<data_type, n_vars>::operator %=(data_type a)
{
   for (auto i = 0; i < n_vars; i++) data[i] %= a;
   return *this;
};

/*!
\author Vladimir Florinski
\date 03/13/2024
\param[in] other Right operand \f$\mathbf{v}_1\f$
\return This simple array modulo divided by the other simple array
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline SimpleArray<data_type, n_vars>& SimpleArray<data_type, n_vars>::operator %=(const SimpleArray& other)
{
   for (auto i = 0; i < n_vars; i++) data[i] %= other.data[i];
   return *this;
};

/*!
\author Vladimir Florinski
\date 03/08/2024
\return \f$|v|^2\f$
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline data_type SimpleArray<data_type, n_vars>::Norm2(void) const
{
   data_type sum = data[0] * data[0];
   for (auto i = 1; i < n_vars; i++) sum += data[i] * data[i];
   return sum;
};

/*!
\author Vladimir Florinski
\date 03/08/2024
\return Sum of components
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline data_type SimpleArray<data_type, n_vars>::Sum(void) const
{
   data_type sum = data[0];
   for (auto i = 1; i < n_vars; i++) sum += data[i];
   return sum;
};

/*!
\author Vladimir Florinski
\date 03/08/2024
\return Product of components
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline data_type SimpleArray<data_type, n_vars>::Prod(void) const
{
   data_type prod = data[0];
   for (auto i = 1; i < n_vars; i++) prod *= data[i];
   return prod;
};

/*!
\author Vladimir Florinski
\date 03/08/2024
\return Smallest component
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline data_type SimpleArray<data_type, n_vars>::Smallest(void) const
{
   data_type smallest = data[0];
   for (auto i = 1; i < n_vars; i++) smallest = std::min(smallest, data[i]);
   return smallest;
};

/*!
\author Vladimir Florinski
\date 03/08/2024
\return Largest component
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline data_type SimpleArray<data_type, n_vars>::Largest(void) const
{
   data_type largest = data[0];
   for (auto i = 1; i < n_vars; i++) largest = std::max(largest, data[i]);
   return largest;
};

/*!
\author Vladimir Florinski
\date 04/25/2024
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline void SimpleArray<data_type, n_vars>::Negate(void)
{
   for (auto i = 0; i < n_vars; i++) data[i] = -data[i];
};

/*!
\author Vladimir Florinski
\date 03/13/2024
\param[in] other right operand \f$\mathbf{v}_1\f$
\return \f$\mathbf{v}\cdot\mathbf{v_1}\f$
*/
template <typename data_type, int n_vars>
SPECTRUM_DEVICE_FUNC inline data_type SimpleArray<data_type, n_vars>::ScalarProd(const SimpleArray& other) const
{
   data_type sum = data[0] * other.data[0];
   for (auto i = 1; i < n_vars; i++) sum += data[i] * other.data[i];
   return sum;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Stream insertion operator (host only)
\author Vladimir Florinski
\date 03/10/2024
\param[in,out] os     Output stream
\param[in]     sarr_o Simple array to output
\return Modified output stream 
*/
template<typename data_type, int n_vars>
inline std::ostream& operator <<(std::ostream& os, const SimpleArray<data_type, n_vars>& sarr_o)
{
   os << "(" << sarr_o[0];
   for (auto i = 1; i < n_vars; i++) os << ", " << sarr_o[i];
   os << ")";
   return os;
};

/*!
\brief Stream extraction operator (host only)
\author Juan G Alonso Guzman
\date 08/11/2024
\param[in,out] is     Input stream
\param[in]     sarr_i Simple array to input
\return Modified input stream 
*/
template<typename data_type, int n_vars>
inline std::istream& operator >>(std::istream& is, SimpleArray<data_type, n_vars>& sarr_i)
{
   for (auto i = 0; i < n_vars; i++) is >> sarr_i[i];
   return is;
};

};

#endif
