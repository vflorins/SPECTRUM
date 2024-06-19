/*!
\file background_waves.cc
\brief Implements a background consisting of a superposition of waves
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_waves.hh"
#include <iostream>

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundDipole methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 08/29/2022
*/
BackgroundWaves::BackgroundWaves(void)
               : BackgroundBase(bg_name_waves, 0, MODEL_STATIC)
{
};

/*!
\author Vladimir Florinski
\date 08/29/2022
\param[in] other Object to initialize from

A copy constructor should first first call the Params' version to copy the data container and then check whether the other object has been set up. If yes, it should simply call the virtual method "SetupBackground()" with the argument of "true".
*/
BackgroundWaves::BackgroundWaves(const BackgroundWaves& other)
               : BackgroundBase(other)
{
   RAISE_BITS(_status, MODEL_STATIC);
   if(BITS_RAISED(other._status, STATE_SETUP_COMPLETE)) SetupBackground(true);
};

/*!
\author Vladimir Florinski
\date 08/30/2022
\param [in] construct Whether called from a copy constructor or separately

This method's main role is to unpack the data container and set up the class data members and status bits marked as "persistent". The function should assume that the data container is available because the calling function will always ensure this.
*/
void BackgroundWaves::SetupBackground(bool construct)
{
   int dim, wave;
   turb_type t_type;
   double kn, phi, sint, cost, sinp, cosp, alpha, dlnk, psd, psd_tot;
   GeoMatrix basis_b, basis_k;
   TurbProp properties;

// The parent version must be called explicitly if not constructing
   if(!construct) BackgroundBase::SetupBackground(false);

// Build the field-aliogned coordinate system
   basis_b.AxisymmetricBasis(B0);

// Two auxilliary coordinate systems are used here. The first is the B-frame, where e_3 is parallel to B, e_1 is an arbitrary vector normal to e_3, and e_2=e_3^e_1. The second is the K-frame, where e_3 is parallel to k, e_2 is parallel to B^k, and e_1=e_2^e_3. The fluctuating field lies in the K-plane to satisfy the divergence-free condition. The end result is the matrix "basis" that performs a transformation from the K-fram to the global frame.
   shortest_wave = dmax0 * sp_large;
   for(t_type = turb_alfven; t_type <= turb_isotropic; GEO_INCR(t_type, turb_type)) {
      container.Read(&properties);
      dlnk = log(properties.kmax / properties.kmin) / (properties.n_waves - 1.0);
      n_waves[t_type] = properties.n_waves;
      shortest_wave = fmin(M_2PI / properties.kmax, shortest_wave);

// Make all vectors zero length
      psd_tot = 0.0;
      Ampl[t_type].clear();
      k[t_type].clear();
      cosa[t_type].clear();
      sina[t_type].clear();
      phase[t_type].clear();
      basis[t_type].clear();

// Calculate the properties of each wave mode
      for(wave = 0; wave < n_waves[t_type]; wave++) {
         kn = properties.kmin * pow(properties.kmax / properties.kmin, wave / (properties.n_waves - 1.0));
         k[t_type].push_back(kn);
         phase[t_type].push_back(M_2PI * rng->GetUniform());

// "cost" is the cosine of the angle between k and B, "phi" is the angle of k_perp in the B-plane, "alpha" is the angle of deltaB in the K-plane.
         switch(t_type) {

// Parallel propagation (same as anti-parallel for static field), phi is not relevant, alpha is random in the K-plane
         case turb_alfven:
            dim = 0;
            cost = 1.0;
            phi = 0.0;
            alpha = M_2PI * rng->GetUniform();
            break;

// Perpendicular propagation, phi is random in the B-plane, alpha is pi/2 or -pi/2 in the K-plane
         case turb_transverse:
            dim = 1;
            cost = 0.0;
            phi = M_2PI * rng->GetUniform();

// The randomness of phase implies this can be fixed
//            alpha = 0.5 * M_PI * (rng->GetUniform() > 0.5 ? 1.0 : -1.0);
            alpha = 0.5 * M_PI;
            break;

// Perpendicular propagation, phi is random in the B-plane, alpha is 0 in the K-plane
         case turb_longitudinal:
            dim = 1;
            cost = 0.0;
            phi = M_2PI * rng->GetUniform();
            alpha = 0.0;
            break;

// Arbitrary propagation, alpha is random in the K-plane
         case turb_isotropic:
            dim = 2;
            cost = 2.0 * (rng->GetUniform() - 0.5);
            phi = M_2PI * rng->GetUniform();
            alpha = M_2PI * rng->GetUniform();
            break;

         default:
            break;         
         };

// No need to normalize PSD because later we calculate the weights as fractions of the total. The factor of "k" comes from logarithmic spacing of wavenumbers (regardless of geometry).
//         psd = pow(kn, 1.0 - properties.slope);
         psd = pow(kn, dim + 1) / (1.0 + pow(kn * properties.l0, properties.slope + dim));
         psd_tot += psd;

// Precompute the sine and cosine of the polarization angle
         cosa[t_type].push_back(cos(alpha));
         sina[t_type].push_back(sin(alpha));
         sint = sqrt(1.0 - Sqr(cost));
         cosp = cos(phi);
         sinp = sin(phi);

// Compute the unit vectors of the K-frame in the B-frame
         basis_k[0][0] = cost * cosp;
         basis_k[0][1] = cost * sinp;
         basis_k[0][2] = -sint;

         basis_k[1][0] = -sinp;
         basis_k[1][1] = cosp;
         basis_k[1][2] = 0.0;

         basis_k[2][0] = sint * cosp;
         basis_k[2][1] = sint * sinp;
         basis_k[2][2] = cost;

// Convert back to the global frame
         basis_k.ChangeFromBasis(basis_b);
         basis[t_type].push_back(basis_k);

// Unnormalized amplitude
         Ampl[t_type].push_back(sqrt(psd));
      };

// Normalize the amplitudes. The factor of 2 comes from using only cosine for the waves (variance is 1/2)
      for(wave = 0; wave < n_waves[t_type]; wave++) Ampl[t_type][wave] *= sqrt(2.0 * properties.variance / psd_tot);
   };
};

/*!
\author Vladimir Florinski
\date 10/14/2022
*/
void BackgroundWaves::EvaluateBackground(void)
{
   int wave;
   double arg, z_rot;
   turb_type t_type;
   GeoVector posprime, B_rot;

   posprime = _pos - r0;

   if(BITS_RAISED(_spdata._mask, BACKGROUND_U)) _spdata.Uvec = 0.0;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_B)) {
      _spdata.Bvec = B0;

      for(t_type = turb_alfven; t_type <= turb_isotropic; GEO_INCR(t_type, turb_type)) {
         for(wave = 0; wave < n_waves[t_type]; wave++) {

// Compute the z coordinate in the rotated (wave) frame. This is the direction of the wave vector.
            z_rot = posprime * basis[t_type][wave][2];

// Compute the magnetic field components in the rotated frame. Because the divergence of B is zero, we only have B_x and B_y in this frame.
            arg = k[t_type][wave] * z_rot + phase[t_type][wave];
            B_rot[0] =  Ampl[t_type][wave] * cosa[t_type][wave] * cos(arg);
            B_rot[1] = -Ampl[t_type][wave] * sina[t_type][wave] * sin(arg);

// Project the field back into the global frame. This is faster than calling "ChangeFromBasis()" because it saves 6 ops out of 15
            _spdata.Bvec[0] += B_rot[0] * basis[t_type][wave][0][0] + B_rot[1] * basis[t_type][wave][1][0];
            _spdata.Bvec[1] += B_rot[0] * basis[t_type][wave][0][1] + B_rot[1] * basis[t_type][wave][1][1];
            _spdata.Bvec[2] += B_rot[0] * basis[t_type][wave][0][2] + B_rot[1] * basis[t_type][wave][1][2];
         };
      };
   };
   if(BITS_RAISED(_spdata._mask, BACKGROUND_E)) _spdata.Evec = 0.0;
   _spdata.region = 1.0;

   LOWER_BITS(_status, STATE_INVALID);
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/14/2022
*/
void BackgroundWaves::EvaluateBackgroundDerivatives(void)
{
   int wave, xyz;
   double arg, z_rot;
   turb_type t_type;
   GeoVector posprime, dB_rot_dz, dB_rot;

   posprime = _pos - r0;

   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradU)) _spdata.gradUvec = gm_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradB)) {
      _spdata.gradBvec = gm_zeros;

      for(t_type = turb_alfven; t_type <= turb_isotropic; GEO_INCR(t_type, turb_type)) {
         for(wave = 0; wave < n_waves[t_type]; wave++) {

// Compute the z coordinate in the rotated (wave) frame. This is the direction of the wave vector.
            z_rot = posprime * basis[t_type][wave][2];

// Compute the z derivatives of the magnetic field components in the rotated frame.
            arg = k[t_type][wave] * z_rot + phase[t_type][wave];
            dB_rot_dz[0] = -Ampl[t_type][wave] * cosa[t_type][wave] * k[t_type][wave] * sin(arg);
            dB_rot_dz[1] = -Ampl[t_type][wave] * sina[t_type][wave] * k[t_type][wave] * cos(arg);

            dB_rot[0] = dB_rot_dz[0] * basis[t_type][wave][0][0] + dB_rot_dz[1] * basis[t_type][wave][1][0];
            dB_rot[1] = dB_rot_dz[0] * basis[t_type][wave][0][1] + dB_rot_dz[1] * basis[t_type][wave][1][1];
            dB_rot[2] = dB_rot_dz[0] * basis[t_type][wave][0][2] + dB_rot_dz[1] * basis[t_type][wave][1][2];

// Transform the derivataives to the global frame
            for(xyz = 0; xyz < 3; xyz++) {
               _spdata.gradBvec[0][xyz] += dB_rot[xyz] * basis[t_type][wave][3][0];
               _spdata.gradBvec[1][xyz] += dB_rot[xyz] * basis[t_type][wave][3][1];
               _spdata.gradBvec[2][xyz] += dB_rot[xyz] * basis[t_type][wave][3][2];
            };
         };
      };
   };
   if(BITS_RAISED(_spdata._mask, BACKGROUND_gradE)) _spdata.gradEvec = gm_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_dUdt)) _spdata.dUvecdt = gv_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_dBdt)) _spdata.dBvecdt = gv_zeros;
   if(BITS_RAISED(_spdata._mask, BACKGROUND_dEdt)) _spdata.dEvecdt = gv_zeros;
};

/*!
\author Vladimir Florinski
\date 10/14/2022
*/
void BackgroundWaves::EvaluateDmax(void)
{
   _spdata.dmax = fmin(shortest_wave, dmax0);
   LOWER_BITS(_status, STATE_INVALID);
};

};
