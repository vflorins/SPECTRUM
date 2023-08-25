/*!
\file multi_index.cc
\brief Defines a three-component index class
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "common/multi_index.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Friend methods of MultiIndex
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\date 06/16/2020
\param[in] midx_l Left operand
\param[in] midx_r Right operand
\return Component-wise sum of the two multi-indices
*/
SPECTRUM_DEVICE_FUNC MultiIndex operator +(const MultiIndex& midx_l, const MultiIndex& midx_r)
{
   MultiIndex midx_tmp(midx_l);
   midx_tmp.i += midx_r.i;
   midx_tmp.j += midx_r.j;
   midx_tmp.k += midx_r.k;
   return midx_tmp;
};

/*!
\author Vladimir Florinski
\date 06/09/2020
\param[in] midx_l Multi-index to be incremented
\param[in] sclr_r    Value of the increment
\return Incremented multi-index
*/
SPECTRUM_DEVICE_FUNC MultiIndex operator +(const MultiIndex& midx_l, int sclr_r)
{
   MultiIndex midx_tmp(midx_l);
   midx_tmp.i += sclr_r;
   midx_tmp.j += sclr_r;
   midx_tmp.k += sclr_r;
   return midx_tmp;
};

/*!
\author Vladimir Florinski
\date 06/09/2020
\param[in] sclr_l    Value of the increment
\param[in] midx_r Multi-index to be incremented
\return Incremented multi-index
*/
SPECTRUM_DEVICE_FUNC MultiIndex operator +(int sclr_l, const MultiIndex& midx_r)
{
   MultiIndex midx_tmp(midx_r);
   midx_tmp.i += sclr_l;
   midx_tmp.j += sclr_l;
   midx_tmp.k += sclr_l;
   return midx_tmp;
};

/*!
\author Vladimir Florinski
\date 06/16/2020
\param[in] midx_l Left operand
\param[in] midx_r Right operand
\return Component-wise difference of the two multi-indices
*/
SPECTRUM_DEVICE_FUNC MultiIndex operator -(const MultiIndex& midx_l, const MultiIndex& midx_r)
{
   MultiIndex midx_tmp(midx_l);
   midx_tmp.i -= midx_r.i;
   midx_tmp.j -= midx_r.j;
   midx_tmp.k -= midx_r.k;
   return midx_tmp;
};

/*!
\author Vladimir Florinski
\date 06/09/2020
\param[in] midx_l Left operand
\param[in] sclr_r Right operand
\return Decremented multi-index
*/
SPECTRUM_DEVICE_FUNC MultiIndex operator -(const MultiIndex& midx_l, int sclr_r)
{
   MultiIndex midx_tmp(midx_l);
   midx_tmp.i -= sclr_r;
   midx_tmp.j -= sclr_r;
   midx_tmp.k -= sclr_r;
   return midx_tmp;
};

/*!
\author Vladimir Florinski
\date 06/16/2020
\param[in] midx_l Left operand
\param[in] sclr_r Right operand
\return Multi-index with each component multiplied by "sclr_r"
*/
SPECTRUM_DEVICE_FUNC MultiIndex operator *(const MultiIndex& midx_l, int sclr_r)
{
   MultiIndex midx_tmp(midx_l);
   midx_tmp.i *= sclr_r;
   midx_tmp.j *= sclr_r;
   midx_tmp.k *= sclr_r;
   return midx_tmp;
};

/*!
\author Vladimir Florinski
\date 06/16/2020
\param[in] sclr_l Left operand
\param[in] midx_r Right operand
\return Multi-index with each component multiplied by "sclr_l"
*/
SPECTRUM_DEVICE_FUNC MultiIndex operator *(int sclr_l, const MultiIndex& midx_r)
{
   MultiIndex midx_tmp(midx_r);
   midx_tmp.i *= sclr_l;
   midx_tmp.j *= sclr_l;
   midx_tmp.k *= sclr_l;
   return midx_tmp;
};

/*!
\author Vladimir Florinski
\date 06/16/2020
\param[in] midx_l Left operand
\param[in] sclr_r Right operand
\return Multi-index with each component divided by "sclr_r"
*/
SPECTRUM_DEVICE_FUNC MultiIndex operator /(const MultiIndex& midx_l, int sclr_r)
{
   MultiIndex midx_tmp(midx_l);
   midx_tmp.i /= sclr_r;
   midx_tmp.j /= sclr_r;
   midx_tmp.k /= sclr_r;
   return midx_tmp;
};

/*!
\author Vladimir Florinski
\date 06/16/2020
\param[in] midx_l The divisor
\param[in] midx_l The divis
\return Modulo-divided multi-index

\note Suitable for negative values of "midx_l".
*/
SPECTRUM_DEVICE_FUNC MultiIndex operator %(const MultiIndex& midx_l, const MultiIndex& midx_r)
{
   MultiIndex midx_tmp = midx_l + midx_r;
   midx_tmp.i %= midx_r.i;
   midx_tmp.j %= midx_r.j;
   midx_tmp.k %= midx_r.k;
   return midx_tmp;
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] midx_l Left operand
\param[in] midx_r Right operand
\return True if all three components of "midx_l" are larger than those of "midx_r"
*/
SPECTRUM_DEVICE_FUNC bool operator >(const MultiIndex& midx_l, const MultiIndex& midx_r)
{
   return ((midx_l.i > midx_r.i) && (midx_l.j > midx_r.j) && (midx_l.k > midx_r.k));
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] midx_l Left operand
\param[in] midx_r Right operand
\return True if all three components of "midx_l" are smaller than those of "midx_r"
*/
SPECTRUM_DEVICE_FUNC bool operator <(const MultiIndex& midx_l, const MultiIndex& midx_r)
{
   return ((midx_l.i < midx_r.i) && (midx_l.j < midx_r.j) && (midx_l.k < midx_r.k));
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] midx_l Left operand
\param[in] midx_r Right operand
\return True if all three components of "midx_l" are greater or equal to those of "midx_r"
*/
SPECTRUM_DEVICE_FUNC bool operator >=(const MultiIndex& midx_l, const MultiIndex& midx_r)
{
   return ((midx_l.i >= midx_r.i) && (midx_l.j >= midx_r.j) && (midx_l.k >= midx_r.k));
};

/*!
\author Vladimir Florinski
\date 06/14/2021
\param[in] midx_l Left operand
\param[in] midx_r Right operand
\return True if all three components of "midx_l" are smaller or equal to those of "midx_r"
*/
SPECTRUM_DEVICE_FUNC bool operator <=(const MultiIndex& midx_l, const MultiIndex& midx_r)
{
   return ((midx_l.i <= midx_r.i) && (midx_l.j <= midx_r.j) && (midx_l.k <= midx_r.k));
};

/*!
\author Vladimir Florinski
\date 06/15/2021
\param[in] midx_l Left operand
\param[in] midx_r Right operand
\return True if all three components of "midx_l" are equal to those of "midx_r"
*/
SPECTRUM_DEVICE_FUNC bool operator ==(const MultiIndex& midx_l, const MultiIndex& midx_r)
{
   return ((midx_l.i == midx_r.i) && (midx_l.j == midx_r.j) && (midx_l.k == midx_r.k));
};

/*!
\author Vladimir Florinski
\date 06/15/2021
\param[in] midx_l Left operand
\param[in] midx_r Right operand
\return True if any of the three components of "midx_l" is not equal to those of "midx_r"
*/
SPECTRUM_DEVICE_FUNC bool operator !=(const MultiIndex& midx_l, const MultiIndex& midx_r)
{
   return ((midx_l.i != midx_r.i) || (midx_l.j != midx_r.j) || (midx_l.k != midx_r.k));
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Other methods operating on multi-indices
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 05/11/2022
\param[in] midx_r right operand
\return Modified "ostream" object
*/
std::ostream& operator <<(std::ostream& os, const MultiIndex& midx_r)
{
   os << "(" << midx_r.i << ", " << midx_r.j << ", " << midx_r.k << ")";
   return os;
};

};
