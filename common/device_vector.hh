/*!
\file device_vector.hh
\brief A simple std::vector replacement for the device
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_DEVICE_VECTOR_HH
#define SPECTRUM_DEVICE_VECTOR_HH

#include "common/gpu_config.hh"

namespace Spectrum {

#ifdef __CUDACC__

template <typename T> class device_vector {

private:

//! Beginning address
   T* m_begin;

//! End address
   T* m_end;

//! Current capacity
   size_t capacity;

//! Current size
   size_t length;

    __device__ void expand() {
        capacity *= 2;
        size_t tempLength = (m_end - m_begin);
        T* tempBegin = new T[capacity];

        memcpy(tempBegin, m_begin, tempLength * sizeof(T));
        delete[] m_begin;
        m_begin = tempBegin;
        m_end = m_begin + tempLength;
        length = static_cast<size_t>(m_end - m_begin);
    }

public:

//! Default constructor
   __device__ device_vector(void);

//! Element access
   __device__ T& operator[](unsigned int index);

//! Size of the vector
   __device__ size_type size(void);

    __device__ T* begin() {
        return m_begin;
    }
    __device__ T* end() {
        return m_end;
    }
    __device__ ~LocalVector()
    {
        delete[] m_begin;
        m_begin = nullptr;
    }

    __device__ void add(T t) {

        if ((m_end - m_begin) >= capacity) {
            expand();
        }

        new (m_end) T(t);
        m_end++;
        length++;
    }
    __device__ T pop() {
        T endElement = (*m_end);
        delete m_end;
        m_end--;
        return endElement;
    }

    __device__ size_t getSize() {
        return length;
    }
};

#endif

};

#endif
