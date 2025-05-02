/*!
\file spherical_slab.hh
\brief Declares slab class, a range of spherical shells
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_SPHERICAL_SLAB_HH
#define SPECTRUM_SPHERICAL_SLAB_HH

#include <common/print_warn.hh>

namespace Spectrum {

//! Smallest height of a block
#define min_block_height 2

//! Largest height of a block
#define max_block_height 1000

//! Smallest number of ghost shells
#define min_ghost_height 1

//! Largest number of ghost shells
#define max_ghost_height 5

/*!
\brief A class describing slabs in the mesh
\author Vladimir Florinski

This is a base class for grid blocks that provides shell and interface indexing.
*/
class SphericalSlab
{
protected:

//! Number of interior r-shells
   int n_shells = -1;

//! Total number of r-shells, including ghost
   int n_shells_withghost;

//! Number of interior r-ifaces
   int n_ifaces;

//! Total number of r-ifaces, including ghost
   int n_ifaces_withghost;

//! Number of ghost shells on each side of the slab
   int ghost_height;

//! Determine whether a shell index falls in a given range
   SPECTRUM_DEVICE_FUNC bool IsInteriorShell(int kstart, int height, int k) const;

//! Determine whether a shell belongs to the interior of the slab
   SPECTRUM_DEVICE_FUNC bool IsInteriorShellOfSlab(int k) const;

public:

//! Default constructor
   SPECTRUM_DEVICE_FUNC SphericalSlab(void);

//! Copy constructor
   SPECTRUM_DEVICE_FUNC SphericalSlab(const SphericalSlab& other);

//! Move constructor
   SPECTRUM_DEVICE_FUNC SphericalSlab(SphericalSlab&& other) noexcept;

//! Constructor with parameters
   SPECTRUM_DEVICE_FUNC SphericalSlab(int height, int hghost);

//! Destructor
   SPECTRUM_DEVICE_FUNC ~SphericalSlab(void);

//! Set the slab dimensions
   SPECTRUM_DEVICE_FUNC void SetDimensions(int height, int hghost, bool construct);

//! Return the number of r-shells excluding ghost
   SPECTRUM_DEVICE_FUNC int InteriorShells(void) const {return n_shells;};

//! Return the number of r-shells including ghost
   SPECTRUM_DEVICE_FUNC int TotalShells(void) const {return n_shells_withghost;};
};

/*!
\author Vladimir Florinski
\date 03/07/2025
*/
SPECTRUM_DEVICE_FUNC inline SphericalSlab::SphericalSlab(void)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Default constructing a SphericalSlab\n";
#endif
#endif
};

/*!
\author Vladimir Florinski
\date 01/08/2025
\param[in] other Object to initialize from
*/
SPECTRUM_DEVICE_FUNC inline SphericalSlab::SphericalSlab(const SphericalSlab& other)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Copy constructing a SphericalSlab\n";
#endif
#endif
   if(other.n_shells != -1) SetDimensions(other.n_shells, other.ghost_height, true);
};

/*!
\author Vladimir Florinski
\date 01/08/2025
\param[in] other Object to move into this
*/
SPECTRUM_DEVICE_FUNC inline SphericalSlab::SphericalSlab(SphericalSlab&& other) noexcept
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Move constructing a SphericalSlab\n";
#endif
#endif

   if(other.n_shells == -1) return;

   SetDimensions(other.n_shells, other.ghost_height, true);
   other.n_shells = -1;
};

/*!
\author Vladimir Florinski
\date 07/22/2019
\param[in] width  Height of the slab, without ghost cells
\param[in] wgohst Height of the ghost cell layer outside the slab
*/
SPECTRUM_DEVICE_FUNC inline SphericalSlab::SphericalSlab(int height, int hghost)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Argument constructing a SphericalSlab\n";
#endif
#endif

   SetDimensions(height, hghost, true);
};

/*!
\author Vladimir Florinski
\date 03/07/2025
*/
SPECTRUM_DEVICE_FUNC inline SphericalSlab::~SphericalSlab(void)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Destructing a SphericalSlab\n";
#endif
#endif
};

/*!
\author Vladimir Florinski
\date 02/19/2020
\param[in] width     Height of the slab, without ghost cells
\param[in] hghost    Height of the ghost cell layer outside the slab
\param[in] construct Set to true when called from a constructor
*/
SPECTRUM_DEVICE_FUNC inline void SphericalSlab::SetDimensions(int height, int hghost, bool construct)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Setting dimension of " << height << " for a SphericalSlab\n";
#endif
#endif

   if ((height < min_block_height) || (height > max_block_height) || (hghost < min_ghost_height) || (hghost > max_ghost_height)) {
      PrintError(__FILE__, __LINE__, "Cannot allocate a slab with these dimensions", true);
      return;
   };

// Set up the grid dimensions
   n_shells = height;
   n_ifaces = n_shells + 1;
   ghost_height = hghost;
   n_shells_withghost = n_shells + 2 * ghost_height;
   n_ifaces_withghost = n_shells_withghost + 1;
};

/*!
\author Vladimir Florinski
\date 07/22/2019
\param[in] kstart Starting shell of the sub-block
\param[in] height Height of the sub-block
\param[in] k      Shell index
\return True if the shell is between "kstart" and "kstart+height"
*/
SPECTRUM_DEVICE_FUNC inline bool SphericalSlab::IsInteriorShell(int kstart, int height, int k) const
{
   return (k >= kstart) && (k < kstart + height);
};

/*!
\author Vladimir Florinski
\date 07/22/2019
\param[in] k Shell index
\return True if the shell is interior to the slab
*/
SPECTRUM_DEVICE_FUNC inline bool SphericalSlab::IsInteriorShellOfSlab(int k) const
{
   return (k >= ghost_height) && (k <= n_shells_withghost - 1 - ghost_height);
};

};

#endif
