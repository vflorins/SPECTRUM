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

//! Smallest number of ghost shells
#define min_ghost_height 1

//! Largest number of ghost shells
#define max_ghost_height 4

//! Smallest ratio between block height and ghost height
#define min_height_to_ghost 2

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
   bool IsInteriorShell(int kstart, int height, int k) const;

//! Determine whether a shell belongs to the interior of the slab
   bool IsInteriorShellOfSlab(int k) const;

public:

//! Default constructor
   SphericalSlab(void);

//! Copy constructor
   SphericalSlab(const SphericalSlab& other);

//! Move constructor
   SphericalSlab(SphericalSlab&& other) noexcept;

//! Constructor with parameters
   SphericalSlab(int height, int hghost);

//! Destructor
   ~SphericalSlab(void);

//! Set the slab dimensions
   void SetDimensions(int height, int hghost, bool construct);

//! Return the number of r-shells excluding ghost
   int InteriorShells(void) const {return n_shells;};

//! Return the number of r-shells including ghost
   int TotalShells(void) const {return n_shells_withghost;};
};

/*!
\author Vladimir Florinski
\date 03/07/2025
*/
inline SphericalSlab::SphericalSlab(void)
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
inline SphericalSlab::SphericalSlab(const SphericalSlab& other)
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
inline SphericalSlab::SphericalSlab(SphericalSlab&& other) noexcept
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Move constructing a SphericalSlab\n";
#endif
#endif

   if (other.n_shells == -1) {
      PrintMessage(__FILE__, __LINE__, "Move constructor called, but the dimension of the moved object was not set", true);
      return;
   };

   SetDimensions(other.n_shells, other.ghost_height, true);
};

/*!
\author Vladimir Florinski
\date 07/22/2019
\param[in] width  Height of the slab, without ghost cells
\param[in] wgohst Height of the ghost cell layer outside the slab
*/
inline SphericalSlab::SphericalSlab(int height, int hghost)
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
inline SphericalSlab::~SphericalSlab(void)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Destructing a SphericalSlab\n";
#endif
#endif
};

/*!
\author Vladimir Florinski
\date 03/21/2025
\param[in] width     Height of the slab, without ghost cells
\param[in] hghost    Height of the ghost cell layer outside the slab
\param[in] construct Set to true when called from a constructor
*/
inline void SphericalSlab::SetDimensions(int height, int hghost, bool construct)
{
#ifdef GEO_DEBUG
#if GEO_DEBUG_LEVEL >= 3
   std::cerr << "Setting dimension of " << height << " for a SphericalSlab\n";
#endif
#endif

   if ((hghost < min_ghost_height) || (hghost > max_ghost_height) || (height < min_height_to_ghost * hghost)) {
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
inline bool SphericalSlab::IsInteriorShell(int kstart, int height, int k) const
{
   return (k >= kstart) && (k < kstart + height);
};

/*!
\author Vladimir Florinski
\date 07/22/2019
\param[in] k Shell index
\return True if the shell is interior to the slab
*/
inline bool SphericalSlab::IsInteriorShellOfSlab(int k) const
{
   return (k >= ghost_height) && (k <= n_shells_withghost - 1 - ghost_height);
};

};

#endif