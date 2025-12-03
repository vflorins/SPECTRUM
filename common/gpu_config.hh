/*!
\file gpu_config.hh
\brief Defines some macros for CUDA (and HIP in the future)
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_GPU_CONFIG_HH
#define SPECTRUM_GPU_CONFIG_HH

#ifdef __CUDACC__

#include <cuda_runtime.h>
#include <iostream>

//! Every function marked so is compiled for both host and device
#define SPECTRUM_DEVICE_FUNC __host__ __device__
#define SPECTRUM_CONSTEXPR __device__ constexpr

/*!
\brief Print the GPU information (from the host)
\author Vladimir Florinski
\date 10/28/2025
\param[in] device The device to print information about
*/
__host__ inline void PrintDeviceInfo(int device = 0)
{
   cudaDeviceProp prop;
   cudaGetDeviceProperties(&prop, device);
   std::cerr << "Device name: " << prop.name << "\n";
   std::cerr << "Compute capability: " << prop.major << "." << prop.minor << "\n";
   std::cerr << "Number of MPs is " << prop.multiProcessorCount << "\n";
   std::cerr << "Device memmory is " << prop.totalGlobalMem / 1024 / 1024 << " MB\n";
   std::cerr << "L2 cache memory is " << prop.l2CacheSize / 1024 / 1024 << " MB\n";
   std::cerr << "Maximum number of threads: " << prop.maxThreadsPerBlock << " per block, "
                                              << prop.maxThreadsPerMultiProcessor << " per MP\n";
   std::cerr << "Register memory: " << prop.regsPerBlock / 256 << " kB per block, "
                                    << prop.regsPerMultiprocessor / 256 << " kB per MP\n";
   std::cerr << "Shared memory: " << prop.sharedMemPerBlock / 1024 << " kB per block, "
                                  << prop.sharedMemPerMultiprocessor / 1024 << " kB per MP\n";
};

#else

//! Define as empty when compiling with a non-GPU compiler
#define SPECTRUM_DEVICE_FUNC
#define SPECTRUM_CONSTEXPR constexpr

/*!
\brief Empty function
\author Vladimir Florinski
\date 04/17/2023
\param[in] device Unused
*/
inline void PrintDeviceInfo(int device = 0)
{
};

#endif

#endif
