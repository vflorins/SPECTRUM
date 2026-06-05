/*!
\file batl_spectrum_interface.hh
\brief C interface for BATL library calls
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef BATL_SPECTRUM_INTERFACE_HH
#define BATL_SPECTRUM_INTERFACE_HH

namespace Spectrum {

#ifdef __cplusplus
extern "C" {
#endif

void wrapamr_clean();
void wrapamr_read_header(const char*, int, int);
void wrapamr_read_file(const char*, int, int, int);
void wrapamr_read_file_partial(const char*, int, int, int, int, int*);
void wrapamr_get_ndim(int*);
void wrapamr_get_nvar(int*);
void wrapamr_get_domain(double*, double*);
void wrapamr_get_block_size(int*, int*, int*, int*);
void wrapamr_get_data_serial(const double* pos, double* variables, int*);

void spectrum_init_mpi(int comm);
void spectrum_get_node(const double* pos, int* node);
void spectrum_get_neighbor_node(int node, int i, int j, int k, int* neighbor_node, int* neighbor_level);
void spectrum_get_interpolation_stencil(const double* pos, int* n_nodes, int* stencil_nodes, int* stencil_zones, double* stencil_weights);

//! Read the parameter file and generate and store a basic non-refined tree
void batl_init_tree_file(const char* filename, int filename_len);

// FIXME debug only
void spectrum_read_header(const char*, int, int);
//void spectrum_write_file(const char*, int);

#ifdef __cplusplus
}
#endif

};

#endif
