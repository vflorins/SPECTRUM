# File config.physical.py Created by Lucius Schoenbaum October 20, 2025

"""
The following dict of dicts define the parameters
and types present in the respective Config class.
Each parameter is also assigned a reasonable default
for the particular trajectory/background/diffusion type.
This file is the unique location in the code where
these default values are defined.

There is no inheritance mechanism, so
each class repeats parameters defined in base classes.
In some cases, this is indicated by commented-out dividers.

The dict of dicts has the following structure:
```
physical_defaults = {
    'Background': ..., # dict of backgrounds
    'Trajectory': ..., # dict of trajectories
    'Diffusion': ..., # dict of diffusions
}
```

"""


"""

TODO: data background parameter definitions::::::

Number of variables per zone 
    int n_variables ...... derived from DataFields
Size of the block in each dimension 
    MultiIndex block_size
Number of neighbors per dimension (depends only on the refinement ratio)
    int max_neighbors_per_dim 
Largest possible number of neighbors (depends only on the refinement ratio)
    int max_neighbors
Largest possible number of neighbor levels (depends only on the refinement ratio)
    int max_neighbor_levels
Number of elements in each stencil
    int stencil_n_elements
the size of the local block cache
    max_cache_size
the number of ghost cells per grid block
    num_ghost_cells
the number of servers per node
    n_servers_per_node
    
n_variables_cartesian 9
block_size_cartesian (4, 4, 4)
max_neighbors_per_dim_cartesian 3
max_neighbors_cartesian 27
max_neighbor_levels_cartesian 27
stencil_n_elements 8
max_cache_size 100
num_ghost_cells 2
n_servers_per_node 1
---------------------------------------------------
n_variables_batl 10
block_size_batl (8, 8, 8)
max_neighbors_per_dim_batl 4
max_neighbors_batl 64
max_neighbor_levels_batl 27
stencil_n_elements 8
max_cache_size 100
num_ghost_cells 2
n_servers_per_node 1



//! Flag to always request stencil directly from BATL for 1st order interpolation (0 = no, 1 = yes)
// todo config
//#define REQUEST_STENCIL_FROM_BATL 0




////! Index of the mass density variable
//// #define SERVER_VAR_INDEX_RHO 0
//
////! Index of the number density variable
//// #define SERVER_VAR_INDEX_DEN 0
//
////! Index of the momentum variable
////#define SERVER_VAR_INDEX_MOM 1
//
////! Index of the bulk flow variable
//#define SERVER_VAR_INDEX_FLO 0
//
////! Index of the magnetic field variable
//#define SERVER_VAR_INDEX_MAG 3
//
////! Index of the electric field variable
//#define SERVER_VAR_INDEX_ELE 6
//
////! Index and number of the regions variable
//// #define SERVER_VAR_INDEX_REG 7
//// #define SERVER_NUM_INDEX_REG 2
//
////! Index of thermal pressure
//// #define SERVER_VAR_INDEX_PTH 9

.......... ??? these values are not self-consistent ......
 a matching Field type::::::











/*
 *
 * Documentation (Work in Progress)
 *
 * The Background methods Evaluate, EvaluateDerivatives, and EvaluateDmax
 * are valid for any Coordinate type that includes position and time (Pos_t, Time_t) and
 * magnitude of momentum (AbsMom_t or at least Mom_t).
 * The fields type (Fields) must always contain magnetic field (Mag_t).
 *
 */

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 09/08/2025
\param[in] coords coordinates, any coordinate system providing (time, position, p*) with p* the magnitude of momentum (access via AbsMom), and position in cartesian system.
\param[out] fields All fields data requested by caller. Optional type RequestedFields specifies which fields to evaluate, if only a subset is needed.
\note This is a common routine that the derived classes should not change.
\note This public method is valid for any Coordinate type that includes
position and time (Pos_t, Time_t) and magnitude of momentum (AbsMom_t or at least Mom_t).
The fields type must always contain magnetic field (Mag_t).
Magnetic field magnitude and/or direction can also be tracked but magnetic field is sufficient.
*/








"""


physical_defaults = {
    'Background':{
        'Uniform': {
            'derivative_method': 'analytic',
            'num_numeric_grad_evals': 1,
            'incr_dmax_ratio': 0.0001,
            'dmax0': 0.1, # todo
        },
        'CylindricalObstacle': {
            'derivative_method': 'analytic',
            'num_numeric_grad_evals': 1,
            'incr_dmax_ratio': 0.0001,
            'dmax0': 0.1, # todo
        },
        'Dipole': {
            'derivative_method': 'analytic',
            'num_numeric_grad_evals': 1,
            'incr_dmax_ratio': 0.0001,
            'dmax0': 0.1, # todo
            'dmax_fraction': 0.1, # todo
            'r0': "{0.0, 0.0, 0.0}", # todo
            'u0': "{0.0, 0.0, 0.0}", # todo
            'B0': "{1.0, 1.0, 1.0}", # todo
            'r_ref': 1.0, # todo
        },
        'MagnetizedCylinder': {
            'derivative_method': 'analytic',
            'num_numeric_grad_evals': 1,
            'incr_dmax_ratio': 0.0001,
            'dmax0': 0.1, # todo
        },
        'MagnetizedSphere': {
            'derivative_method': 'analytic',
            'num_numeric_grad_evals': 1,
            'incr_dmax_ratio': 0.0001,
            'dmax0': 0.1, # todo
        },
        'SphericalObstacle': {
            'derivative_method': 'analytic',
            'num_numeric_grad_evals': 1,
            'incr_dmax_ratio': 0.0001,
            'dmax0': 0.1, # todo
            'r_ref': 1.0, # todo
        },
        'Discontinuity': {
            'derivative_method': 'analytic',
            'num_numeric_grad_evals': 1,
            'incr_dmax_ratio': 0.0001,
            'dmax0': 0.1, # todo
        },
        'Shock': {
            'derivative_method': 'analytic',
            'num_numeric_grad_evals': 1,
            'incr_dmax_ratio': 0.0001,
            'dmax0': 0.1, # todo
            ####
            'tanh_width_factor': 4.0,
        },
        'SmoothDiscontinuity': {
            'derivative_method': 'analytic',
            'num_numeric_grad_evals': 1,
            'incr_dmax_ratio': 0.0001,
            'dmax0': 0.1, # todo
            ####
            'smooth_discontinuity_order': 4,
            'tanh_width_factor': 4.0,
        },
        'SmoothShock': {
            'derivative_method': 'analytic',
            'num_numeric_grad_evals': 1,
            'incr_dmax_ratio': 0.0001,
            'dmax0': 0.1, # todo
            ####
            'smooth_discontinuity_order': 4,
            'tanh_width_factor': 4.0,
        },
        'SolarWind': {
            'derivative_method': 'numeric',
            'num_numeric_grad_evals': 1,
            # todo discuss
            # 'dmax0':  "0.1 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid",
            # 'dmax_fraction': 0.1,
            ####
            'solarwind_current_sheet': 'disabled',
            'solarwind_sectored_region': 'nowhere',
            'solarwind_polar_correction': 'none',
            'solarwind_speed_latitude_profile': 'constant',
            'with_termination_shock': False,
            'termshock_speed_exponent': 'square',
        },
        'VLISMBochum': {
            'derivative_method': 'numeric',
            'num_numeric_grad_evals': 1,
            'incr_dmax_ratio': 0.0001,
            'dmax0': 0.1, # todo
            ####
            'mod_type': 'scaled',
            'mod_rpos': 'scale_rel_zero',
            'z_nose': 1.0, # todo
        },
        'Waves': {
            'derivative_method': 'analytic',
            'num_numeric_grad_evals': 1,
            'incr_dmax_ratio': 0.0001,
            'dmax0': 0.1, # todo
        },
        'Server': {
            'derivative_method': 'numeric',
            'num_numeric_grad_evals': 1,
            'incr_dmax_ratio': 0.0001,
            'dmax0': 0.1, # todo
            ####
            'server_num_ghost_cells': 2,
            'server_interpolation_order': 1,
        },
        'DataBATL': {
            'derivative_method': 'numeric',
            'num_numeric_grad_evals': 1,
            'incr_dmax_ratio': 0.0001,
            'dmax0': 0.1, # todo
            ####
            'server_num_ghost_cells': 2,
            'server_interpolation_order': 1,
        },
        'DataCartesian': {
            'derivative_method': 'numeric',
            'num_numeric_grad_evals': 1,
            'incr_dmax_ratio': 0.0001,
            'dmax0': 0.1, # todo
            ####
            'server_num_ghost_cells': 2,
            'server_interpolation_order': 1,
        },
    },
    'Trajectory': {
        'Fieldline': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::anisotropic>, Pos_t, Time_t, Mom_t, Vel_t>",
            'RecordCoordinates': "Fields<FConfig<>, Pos_t, Time_t>",
            'Fields': "Fields<FConfig<specieid_>>",
            'FieldlineField_t': "Mag_t",
            'time_flow': 'forward',
            'rk_integrator': "DormandPrince_54E",
            'record_mag_extrema': False,
            'record_trajectory': False,
            'record_trajectory_segment_presize': 10000,
            'advance_safety_level': 'low',
            'max_trajectory_steps': 100000,
            'max_time_adaptations': 1,
            'n_max_calls': -1,
            ####
        },
        'Focused': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Mom_t, Vel_t>",
            'RecordCoordinates': "Fields<FConfig<>, Pos_t, Time_t>",
            'Fields': "Fields<FConfig<specieid_>, Fluv_t, Mag_t, AbsMag_t, HatMag_t, DelFluv_t, DelAbsMag_t, DelMag_t, DotFluv_t>",
            'time_flow': 'forward',
            'rk_integrator': "DormandPrince_54E",
            'record_mag_extrema': False,
            'record_trajectory': False,
            'record_trajectory_segment_presize': 10000,
            'advance_safety_level': 'low',
            'max_trajectory_steps': 100000,
            'max_time_adaptations': 1,
            'n_max_calls': -1,
            ####
            'pperp_method': "scheme",
            'use_B_drifts': "none",
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
        },
        'Guiding': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::anisotropic>, Pos_t, Time_t, Mom_t, Vel_t>",
            'RecordCoordinates': "Fields<FConfig<>, Pos_t, Time_t>",
            'Fields': "Fields<FConfig<specieid_>, Fluv_t, Mag_t, Elc_t, AbsMag_t, HatMag_t, DelMag_t, DelAbsMag_t, DotAbsMag_t>",
            'time_flow': 'forward',
            'rk_integrator': "DormandPrince_54E",
            'record_mag_extrema': False,
            'record_trajectory': False,
            'record_trajectory_segment_presize': 10000,
            'advance_safety_level': 'low',
            'max_trajectory_steps': 100000,
            'max_time_adaptations': 1,
            'n_max_calls': -1,
            ####
            'pperp_method': "scheme",
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
        },
        'GuidingDiff': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::anisotropic>, Pos_t, Time_t, Mom_t, Vel_t>",
            'RecordCoordinates': "Fields<FConfig<>, Pos_t, Time_t>",
            'Fields': "Fields<FConfig<specieid_>, Fluv_t, Mag_t, Elc_t, AbsMag_t, HatMag_t, DelMag_t, DelAbsMag_t, DotAbsMag_t>",
            'time_flow': 'forward',
            'rk_integrator': "DormandPrince_54E",
            'record_mag_extrema': False,
            'record_trajectory': False,
            'record_trajectory_segment_presize': 10000,
            'advance_safety_level': 'low',
            'max_trajectory_steps': 100000,
            'max_time_adaptations': 1,
            'n_max_calls': -1,
            ####
            'pperp_method': "scheme",
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
        },
        'GuidingScatt': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::anisotropic>, Pos_t, Time_t, Mom_t, Vel_t>",
            'RecordCoordinates': "Fields<FConfig<>, Pos_t, Time_t>",
            'Fields': "Fields<FConfig<specieid_>, Fluv_t, Mag_t, Elc_t, AbsMag_t, HatMag_t, DelMag_t, DelAbsMag_t, DotAbsMag_t>",
            'time_flow': 'forward',
            'rk_integrator': "DormandPrince_54E",
            'record_mag_extrema': False,
            'record_trajectory': False,
            'record_trajectory_segment_presize': 10000,
            'advance_safety_level': 'low',
            'max_trajectory_steps': 100000,
            'max_time_adaptations': 1,
            'n_max_calls': -1,
            ####
            'pperp_method': "scheme",
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
            'split_scatt_fraction': 0.0,
            'const_dmumax': "constant_dtheta_max",
            'stochastic_method': "Euler",
            'cfl_pitchangle': 0.5,
        },
        'GuidingDiffScatt': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::anisotropic>, Pos_t, Time_t, Mom_t, Vel_t>",
            'RecordCoordinates': "Fields<FConfig<>, Pos_t, Time_t>",
            'Fields': "Fields<FConfig<specieid_>, Fluv_t, Mag_t, Elc_t, AbsMag_t, HatMag_t, DelMag_t, DelAbsMag_t, DotAbsMag_t>",
            'time_flow': 'forward',
            'rk_integrator': "DormandPrince_54E",
            'record_mag_extrema': False,
            'record_trajectory': False,
            'record_trajectory_segment_presize': 10000,
            'advance_safety_level': 'low',
            'max_trajectory_steps': 100000,
            'max_time_adaptations': 1,
            'n_max_calls': -1,
            ####
            'pperp_method': "scheme",
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
            'split_scatt_fraction': 0.0,
            'const_dmumax': "constant_dtheta_max",
            'stochastic_method_mu': "Euler",
            'stochastic_method_perp': "Euler",
            'cfl_pitchangle': 0.5,
        },
        'Lorentz': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::cartesian>, Pos_t, Time_t, Mom_t, Vel_t>",
            'RecordCoordinates': "Fields<FConfig<>, Pos_t, Time_t>",
            'Fields': "Fields<FConfig<specieid_>, Mag_t, Elc_t>",
            'time_flow': 'forward',
            'rk_integrator': "DormandPrince_54E",
            'record_mag_extrema': False,
            'record_trajectory': False,
            'record_trajectory_segment_presize': 10000,
            'advance_safety_level': 'low',
            'max_trajectory_steps': 100000,
            'max_time_adaptations': 1,
            'n_max_calls': -1,
            ####
            'pperp_method': "scheme",
            'cfl_advection': 0.1,
            'drift_safety': 0.5,
            'mirror_threshold': 300,
            'steps_per_orbit': 100,
        },
        'Parker': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Mom_t, Vel_t>",
            'RecordCoordinates': "Fields<FConfig<>, Pos_t, Time_t>",
            'Fields': "Fields<FConfig<specieid_>, Fluv_t, Mag_t, AbsMag_t, HatMag_t, DelMag_t, DelAbsMag_t>",
            'time_flow': 'forward',
            'rk_integrator': "DormandPrince_54E",
            'record_mag_extrema': False,
            'record_trajectory': False,
            'record_trajectory_segment_presize': 10000,
            'advance_safety_level': 'low',
            'max_trajectory_steps': 100000,
            'max_time_adaptations': 1,
            'n_max_calls': -1,
            ####
            'pperp_method': "scheme",
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
            'stochastic_method': "Euler",
            'use_B_drifts': "none",
            'divk_method': "direct",
            'cfl_diffusion': 0.5,
            'cfl_acceleration': 0.5,
            'dlnp_max': 0.01,
        },
    },
    'Diffusion': {
        'None': {
            'Coordinates': "Fields<FConfig<>>",
            'Fields': "Fields<FConfig<>>",
        },
        # todo coords and fields for diffusion types need to be optimized
        # todo defaults for diffusion parameters need to be set
        # todo Flum or AbsFlum in DiffusionFlowMomentumPowerLaw
        # todo those that have indicator fields
        'IsotropicConstant': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>",
            'Fields': "Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>",
            'D0': 1.0,
        },
        'ParaConstant': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>",
            'Fields': "Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>",
            'D0': 1.0,
        },
        'PerpConstant': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>",
            'Fields': "Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>",
            'D0': 1.0,
        },
        'FullConstant': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>",
            'Fields': "Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>",
            'Dperp': 1.0,
            'Dpara': 1.0,
        },
        'QLTConstant': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>",
            'Fields': "Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>",
            'A2A': 1.0,
            'l_max': 1.0,
            'k_min': 1.0,
            'ps_index': 1.0,
            'ps_minus': 1.0,
        },
        'WNLTConstant': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>",
            'Fields': "Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>",
            'use_qlt_scatt': False,
            'A2A': 1.0,
            'l_max': 1.0,
            'k_min': 1.0,
            'ps_index': 1.0,
            'ps_minus': 1.0,
            'A2T': 1.0,
            'A2L': 1.0,
            'ps_plus': 1.0,
        },
        'WNLTRampVLISM': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>",
            'Fields': "Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>",
            'use_qlt_scatt': False,
            'A2A': 1.0,
            'l_max': 1.0,
            'k_min': 1.0,
            'ps_index': 1.0,
            'ps_minus': 1.0,
            'A2T': 1.0,
            'A2L': 1.0,
            'ps_plus': 1.0,
            #
            'k_min_ref': 1.0,
            'A2A_ref': 1.0,
            'A2T_ref': 1.0,
            'A2L_ref': 1.0,
            'l_max_HP': 1.0,
            'dl_max': 1.0,
            'z_nose': 1.0,
            'z_sheath': 1.0,
            'nose_dz': 1.0,
        },
        'FlowMomentumPowerLaw': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>",
            'Fields': "Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t, Fluv_t, AbsFluv_t, DotFluv_t, DelFluv_t>",
            'kappa0': 1.0,
            'U0': 1.0,
            'pow_law_U': 1.0,
            'p0': 1.0,
            'pow_law_p': 1.0,
            'kappa_ratio': 1.0,
        },
        'KineticEnergyRadialDistancePowerLaw': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>",
            'Fields': "Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>",
            'kappa0': 1.0,
            'T0': 1.0,
            'r0_nf': 1.0,
            'pow_law_T': 1.0,
            'pow_law_r': 1.0,
            'kappa_ratio': 1.0,
            'stream_dep_idx': 1,
            'u_upstream': 1.0,
            'w_sh': 1.0,
            's_sh': 1.0,
            'dn_up_ratio': 1.0,
        },
        'RigidityMagneticFieldPowerLaw': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>",
            'Fields': "Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>",
            'lam0': 1.0,
            'R0_nf': 1.0,
            'B0_nf': 1.0,
            'pow_law_R': 1.0,
            'pow_law_B': 1.0,
            'kappa_ratio': 1.0,
        },
        'StraussEtAl2013': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>",
            'Fields': "Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>",
            'LISM_idx': 1,
            'lam_inner': 1.0,
            'lam_outer': 1.0,
            'R0_nf': 1.0,
            'B0_nf': 1.0,
            'kappa_ratio_inner': 1.0,
            'kappa_ratio_outer': 1.0,
        },
        'GuoEtAl2014': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>",
            'Fields': "Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>",
            # same as StraussEtAl2013
            'LISM_idx': 1,
            'lam_inner': 1.0,
            'lam_outer': 1.0,
            'R0_nf': 1.0,
            'B0_nf': 1.0,
            'kappa_ratio_inner': 1.0,
            'kappa_ratio_outer': 1.0,
        },
        'PotgieterEtAl2015': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>",
            'Fields': "Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>",
            'LISM_idx': 1,
            'kappa_inner': 1.0,
            'kappa_outer': 1.0,
            'R0_nf': 1.0,
            'B0_nf': 1.0,
            'kappa_ratio_inner': 1.0,
            'kappa_ratio_outer': 1.0,
        },
        'EmpiricalSOQLTandUNLT': {
            'Coordinates': "Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>",
            'Fields': "Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>",
            'lam_para': 1.0,
            'lam_perp': 1.0,
            'R0_nf': 1.0,
            'B0_nf': 1.0,
            'Bmix_idx': 1,
            'kappa_ratio_red': 1.0,
            'radial_limit_perp_red': 1.0,
            'solar_cycle_idx': 1,
            'solar_cycle_effect': 1.0,
        },
    },
}




