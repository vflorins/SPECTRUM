#include <iostream>
#include <iomanip>
#include <src/background_vlism_bochum.hh>

using namespace Spectrum;

int main(void)
{
   BackgroundVLISMBochum background;

   DataContainer container;
   container.Clear();

// Start time
   double t0 = 0.0;
   container.Insert(t0);

// Position origin
   container.Insert(gv_zeros);

// Velocity reference
   double umag = 25.0E5 / unit_velocity_fluid;
   double theta_u = DegToRad(180.0);
   double phi_u = DegToRad(0.0);
   GeoVector u0(umag, theta_u, phi_u);
   u0.RTP_XYZ();
   container.Insert(u0);

// Magnetic field reference 
   double Bmag = 3.0E-6 / unit_magnetic_fluid;
   double theta_b = DegToRad(40.0);
   double phi_b = DegToRad(0.0);
   GeoVector B0(Bmag, theta_b, phi_b);
   B0.RTP_XYZ();
   container.Insert(B0);

// Required resolution
   double dmax = 0.1 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(dmax);

// Distance to nose
   double z_nose = 110.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(z_nose);

   background.SetupObject(container);

//----------------------------------------------------------------------------------------------------------------------------------------------------

#ifdef USE_SILO

// 3D box plot
   GeoVector xyz_min(-300.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid,
                     -300.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid,
                     -300.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid);
   GeoVector xyz_max( 300.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid,
                      300.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid,
                      300.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid);
   MultiIndex dims_z(200, 200, 200);
   GeoVector normal(0.0, 1.0, 0.0);
   GeoVector right(0.0, 0.0, 1.0);

   background.SetBox(xyz_min, xyz_max, dims_z, normal, right);
   background.BoxPlot3DMesh("box3d.silo", false);
//   background.BoxPlot3DScalar("Bmag", false);
   background.BoxPlot3DVector("Bvec", false);
   background.BoxPlotFinalize();

// 2D box plot
   xyz_min = GeoVector(-400.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid,
                       -400.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid, 0.0);
   xyz_max = GeoVector( 400.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid,
                        400.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid, 0.0);
   dims_z = MultiIndex(2000, 2000, 1);

   background.SetBox(xyz_min, xyz_max, dims_z, normal, right);
   background.BoxPlot2DMesh("box2d.silo", false);
   background.BoxPlot2DScalar("Bmag", false);
//   background.BoxPlot2DVector("Bvec", true);
   background.BoxPlotFinalize();

#endif
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Manual line output
   GeoVector xyz1(0.0, 0.0, 0.0);
   GeoVector xyz2(0.0, 0.0, -500.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid);
   int npts = 500;
   GeoVector pos, dxyz = (xyz2 - xyz1) / npts;
   SpatialData spdata;
   spdata._mask = BACKGROUND_B;
   
   std::cout << std::setprecision(8);
   for (auto i = 0; i <= npts; i++) {
      pos = xyz1 + i * dxyz;

      try {
         background.GetFields(0.0, pos, gv_zeros, spdata);
      }
      catch (ExFieldError& exception) {
         continue;
      };

      std::cout << std::setw(16) << (pos - xyz1).Norm();
      std::cout << std::setw(16) << spdata.Bmag;
      std::cout << std::endl;
   };

   return 0;
};
