/*!
\file background_base_visual.cc
\brief Implements visualization methods of the base background class
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_base.hh"
#include <iostream>
#include <iomanip>

#ifdef USE_SILO

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundBase methods for visualization
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 08/31/2021
\param[in] xyz_min One corner of the box
\param[in] xyz_max Opposite corner of the box
\param[in] dim_z   Number of zones in each direction
\param[in] normal  Normal to the cut plane
\param[in] right   Direction defining the x-axis in the cut plane
*/
void BackgroundBase::SetBox(const GeoVector& xyz_min_in, const GeoVector& xyz_max_in, const MultiIndex& dims_z_in,
                            const GeoVector& normal, const GeoVector& right)
{
   xyz_min = xyz_min_in;
   xyz_max = xyz_max_in;
   dims_z = dims_z_in;
   plane = normal;
   plane.Normalize();
   x_dir = GetSecondUnitVec(plane, right);
};

/*!
\author Vladimir Florinski
\date 08/31/2021
\param[in] box_fname  Name of the SILO file
\param[in] phys_units Use physical units for output
*/
void BackgroundBase::BoxPlot3DMesh(const std::string box_fname, bool phys_units)
{
   if(BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

   int i, xyz;
   GeoVector incr;
   MultiIndex dims_n;
   double* coords[3];

// Memory for nodes
   dims_n = dims_z + 1;
   for(xyz = 0; xyz < 3; xyz++) coords[xyz] = new double[dims_n[xyz]];

// Compute coordinates of mesh points
   incr = (xyz_max - xyz_min) / dims_z;
   for(xyz = 0; xyz < 3; xyz++) {
      for(i = 0; i < dims_n[xyz]; i++) coords[xyz][i] = (xyz_min[xyz] + i * incr[xyz]) * (phys_units ? unit_length_fluid : 1.0);
   };

// Generate SILO output 
   silofile = DBCreate(box_fname.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
   if(!silofile) return;
   DBPutQuadmesh(silofile, mesh3d_name.c_str(), NULL, coords, dims_n.ijk, 3, DB_DOUBLE, DB_COLLINEAR, NULL);

// Clean up
   for(xyz = 0; xyz < 3; xyz++) delete[] coords[xyz];
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
\param[in] var_name   Name of the variable to be plotted
\param[in] phys_units Use physical units for output
*/
void BackgroundBase::BoxPlot3DScalar(const std::string var_name, bool phys_units)
{
   if(BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

   int ix, iy, iz, idx;
   double var_unit;
   bool is_vector;
   long pos;

   GeoVector* field_ptr;
   double* scl_ptr;
   double* scl_field;

// Parse the user input
   if(var_name == "Umag") {
      _spdata._mask = BACKGROUND_U;
      field_ptr = &_spdata.Uvec;
      var_unit = unit_velocity_fluid;
      is_vector = true;
   }
   else if(var_name == "Bmag") {
      _spdata._mask = BACKGROUND_B;
      field_ptr = &_spdata.Bvec;
      var_unit = unit_magnetic_fluid;
      is_vector = true;
   }
   else if(var_name == "Region1") {
      scl_ptr = &_spdata.region[0];
      is_vector = false;
   }
   else if(var_name == "Region2") {
      scl_ptr = &_spdata.region[1];
      is_vector = false;
   }
   else if(var_name == "Region3") {
      scl_ptr = &_spdata.region[2];
      is_vector = false;
   }
   else if(var_name == "n_dens") {
      scl_ptr = &_spdata.n_dens;
      var_unit = unit_number_density_fluid;
      is_vector = false;
   }
   else if(var_name == "p_ther") {
      scl_ptr = &_spdata.p_ther;
      var_unit = unit_pressure_fluid;
      is_vector = false;
   }
   else return;

// Memory for zone var
   scl_field = new double[dims_z.Prod()];
   GeoVector incr = (xyz_max - xyz_min) / dims_z;

   std::cerr << "Calculating the variable " << var_name << ":     ";

// Mesh based output
   idx = 0;
   for(iz = 0; iz < dims_z[2]; iz++) {
      _pos[2] = (xyz_min[2] + (iz + 0.5) * incr[2]);
      for(iy = 0; iy < dims_z[1]; iy++) {
         _pos[1] = (xyz_min[1] + (iy + 0.5) * incr[1]);
         for(ix = 0; ix < dims_z[0]; ix++) {
            _pos[0] = (xyz_min[0] + (ix + 0.5) * incr[0]);
            EvaluateBackground();
            if(is_vector) scl_field[idx] = field_ptr->Norm() * (phys_units ? var_unit : 1.0);
            else scl_field[idx] = *scl_ptr;
            idx++;
         };
      };
      
      std::cerr << "\e[4D";
      std::cerr << std::setw(3) << int(double(iz) / double(dims_z[2]) * 100.0) << "%";
   };
   std::cerr << "\e[4D100%\n";

// Generate SILO output 
   DBPutQuadvar1(silofile, var_name.c_str(), mesh3d_name.c_str(), scl_field, dims_z.ijk, 3, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);

// Clean up
   delete[] scl_field;
};

/*!
\author Vladimir Florinski
\date 10/23/2020
\param[in] var_name   Name of the variable to be plotted
\param[in] phys_units Use physical units for output
*/
void BackgroundBase::BoxPlot3DVector(const std::string var_name, bool phys_units)
{
   if(BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

   int xyz, ix, iy, iz, idx;
   double var_unit;
   GeoVector incr, field;
   GeoVector* field_ptr;
   double* vec_field[3];

// Parse the user input
   if(var_name == "Uvec") {
      _spdata._mask = BACKGROUND_U;
      field_ptr = &_spdata.Uvec;
      var_unit = unit_velocity_fluid;
   }
   else if(var_name == "Bvec") {
      _spdata._mask = BACKGROUND_B;
      field_ptr = &_spdata.Bvec;
      var_unit = unit_magnetic_fluid;
   }
   else return;

// Memory for zone var
   for(xyz = 0; xyz < 3; xyz++) vec_field[xyz] = new double[dims_z.Prod()];
   incr = (xyz_max - xyz_min) / dims_z;

   std::cerr << "Calculating the variable " << var_name << ":     ";

// Mesh based output
   idx = 0;
   for(iz = 0; iz < dims_z[2]; iz++) {
      _pos[2] = (xyz_min[2] + (iz + 0.5) * incr[2]);
      for(iy = 0; iy < dims_z[1]; iy++) {
         _pos[1] = (xyz_min[1] + (iy + 0.5) * incr[1]);
         for(ix = 0; ix < dims_z[0]; ix++) {
            _pos[0] = (xyz_min[0] + (ix + 0.5) * incr[0]);
            EvaluateBackground();
            field = *field_ptr * (phys_units ? var_unit : 1.0);
            vec_field[0][idx] = field[0];
            vec_field[1][idx] = field[1];
            vec_field[2][idx] = field[2];
            idx++;
         };
      };

      std::cerr << "\e[4D";
      std::cerr << std::setw(3) << int(double(iz) / double(dims_z[2]) * 100.0) << "%";
   };
   std::cerr << "\e[4D100%\n";

// Generate SILO output 
   const char* sub_names[] = {(var_name + "_x").c_str(), (var_name + "_y").c_str(), (var_name + "_z").c_str()};
   DBPutQuadvar(silofile, var_name.c_str(), mesh3d_name.c_str(), 3, sub_names, vec_field, dims_z.ijk, 3, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);

// Clean up
   for(xyz = 0; xyz < 3; xyz++) delete[] vec_field[xyz];
};

/*!
\author Vladimir Florinski
\date 08/31/2021
\param[in] box_fname  Name of the SILO file
\param[in] phys_units Use physical units for output
*/
void BackgroundBase::BoxPlot2DMesh(const std::string box_fname, bool phys_units)
{
   if(BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

   int i, xyz;
   GeoVector incr;
   MultiIndex dims_n;
   double* coords[2];

// Memory for nodes
   dims_n = dims_z + 1;
   for(xyz = 0; xyz < 2; xyz++) coords[xyz] = new double[dims_n[xyz]];

// Compute coordinates of mesh points
   incr = (xyz_max - xyz_min) / dims_z;
   for(xyz = 0; xyz < 2; xyz++) {
      for(i = 0; i < dims_n[xyz]; i++) coords[xyz][i] = (xyz_min[xyz] + i * incr[xyz]) * (phys_units ? unit_length_fluid : 1.0);
   };

// Generate SILO output 
   silofile = DBCreate(box_fname.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
   if(!silofile) return;
   DBPutQuadmesh(silofile, mesh2d_name.c_str(), NULL, coords, dims_n.ijk, 2, DB_DOUBLE, DB_COLLINEAR, NULL);

// Clean up
   for(xyz = 0; xyz < 2; xyz++) delete[] coords[xyz];
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 01/04/2024
\param[in] var_name   Name of the variable to be plotted
\param[in] phys_units Use physical units for output
*/
void BackgroundBase::BoxPlot2DScalar(const std::string var_name, bool phys_units)
{
   if(BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

   int xyz, ix, iy, idx;
   double var_unit;
   bool is_vector;

   GeoVector* field_ptr;
   double* scl_ptr;
   GeoVector incr, y_dir, local_coord;
   double* scl_field;

// Parse the user input
   if(var_name == "Umag") {
      _spdata._mask = BACKGROUND_U;
      field_ptr = &_spdata.Uvec;
      var_unit = unit_velocity_fluid;
      is_vector = true;
   }
   else if(var_name == "Bmag") {
      _spdata._mask = BACKGROUND_B;
      field_ptr = &_spdata.Bvec;
      var_unit = unit_magnetic_fluid;
      is_vector = true;
   }
   else if(var_name == "Region1") {
      scl_ptr = &_spdata.region[0];
      is_vector = false;
   }
   else if(var_name == "Region2") {
      scl_ptr = &_spdata.region[1];
      is_vector = false;
   }
   else if(var_name == "Region3") {
      scl_ptr = &_spdata.region[2];
      is_vector = false;
   }
   else if(var_name == "n_dens") {
      scl_ptr = &_spdata.n_dens;
      var_unit = unit_number_density_fluid;
      is_vector = false;
   }
   else if(var_name == "p_ther") {
      scl_ptr = &_spdata.p_ther;
      var_unit = unit_pressure_fluid;
      is_vector = false;
   }
   else return;

// Memory for zone var
   scl_field = new double[dims_z[0] * dims_z[1]];
   incr = (xyz_max - xyz_min) / dims_z;
   y_dir = plane ^ x_dir;

   std::cerr << "Calculating the variable " << var_name << ":     ";

// Mesh based output
   idx = 0;
   for(iy = 0; iy < dims_z[1]; iy++) {
      local_coord[1] = (xyz_min[1] + (iy + 0.5) * incr[1]);
      for(ix = 0; ix < dims_z[0]; ix++) {
         local_coord[0] = (xyz_min[0] + (ix + 0.5) * incr[0]);
         for(xyz = 0; xyz < 3; xyz++) _pos[xyz] = local_coord[0] * x_dir[xyz] + local_coord[1] * y_dir[xyz];
         EvaluateBackground();
         if(is_vector) scl_field[idx] = field_ptr->Norm() * (phys_units ? var_unit : 1.0);
         else scl_field[idx] = *scl_ptr;
         idx++;
      };

      std::cerr << "\e[4D";
      std::cerr << std::setw(3) << int(double(iy) / double(dims_z[1]) * 100.0) << "%";
   };
   std::cerr << "\e[4D100%\n";

// Generate SILO output 
   DBPutQuadvar1(silofile, var_name.c_str(), mesh2d_name.c_str(), scl_field, dims_z.ijk, 2, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);

// Clean up
   delete[] scl_field;
};

/*!
\author Vladimir Florinski
\date 11/02/2020
\param[in] var_name   Name of the variable to be plotted
\param[in] phys_units Use physical units for output
*/
void BackgroundBase::BoxPlot2DVector(const std::string var_name, bool phys_units)
{
   if(BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;

   int xyz, ix, iy, idx;
   double var_unit;
   GeoVector incr, field, y_dir, local_coord;
   GeoVector* field_ptr;
   double* vec_field[2];

// Parse the user input
   if(var_name == "Uvec") {
      _spdata._mask = BACKGROUND_U;
      field_ptr = &_spdata.Uvec;
      var_unit = unit_velocity_fluid;
   }
   else if(var_name == "Bvec") {
      _spdata._mask = BACKGROUND_B;
      field_ptr = &_spdata.Bvec;
      var_unit = unit_magnetic_fluid;
   }
   else return;

// Memory for zone var
   for(xyz = 0; xyz < 2; xyz++) vec_field[xyz] = new double[dims_z[0] * dims_z[1]];
   incr = (xyz_max - xyz_min) / dims_z;
   y_dir = plane ^ x_dir;

   std::cerr << "Calculating the variable " << var_name << ":     ";

// Mesh based output
   idx = 0;
   for(iy = 0; iy < dims_z[1]; iy++) {
      local_coord[1] = (xyz_min[1] + (iy + 0.5) * incr[1]);
      for(ix = 0; ix < dims_z[0]; ix++) {
         local_coord[0] = (xyz_min[0] + (ix + 0.5) * incr[0]);
         for(xyz = 0; xyz < 3; xyz++) _pos[xyz] = local_coord[0] * x_dir[xyz] + local_coord[1] * y_dir[xyz];
         EvaluateBackground();
         field = *field_ptr * (phys_units ? var_unit : 1.0);
         vec_field[0][idx] += field * x_dir;
         vec_field[1][idx] += field * y_dir;
         idx++;
      };

      std::cerr << "\e[4D";
      std::cerr << std::setw(3) << int(double(iy) / double(dims_z[1]) * 100.0) << "%";
   };
   std::cerr << "\e[4D100%\n";

// Generate SILO output 
   const char* sub_names[] = {(var_name + "_x").c_str(), (var_name + "_y").c_str()};
   DBPutQuadvar(silofile, var_name.c_str(), mesh2d_name.c_str(), 2, sub_names, vec_field, dims_z.ijk, 2, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);

// Clean up
   for(xyz = 0; xyz < 2; xyz++) delete[] vec_field[xyz];
};

/*!
\author Vladimir Florinski
\date 10/23/2020
*/
void BackgroundBase::BoxPlotFinalize(void)
{
   if(BITS_LOWERED(_status, STATE_SETUP_COMPLETE)) return;
   if(silofile) DBClose(silofile);
};

};

#endif
