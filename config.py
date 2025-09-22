#!/usr/bin/env python3
"""
File config.py

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.

Configuration file generator for SPECTRUM particle trajectory trace solver.

This script is developed for setting compile-time constants (parameters)
while implementing a system of varying, non-independent, contingent defaults.
It is not possible to do such a thing using compile-time language features
(namely templates and the C preprocessor), as they are currently defined in the standard (C++20).
Therefore a solution is necessary. This is one possibility.

DRAFT
Lucius Schoenbaum
9/13/2025

"""

import argparse
import sys
from math import modf as math_modf


# switch for debugging purposes
# debug = False
debug = True


if debug:
    fname = "config.test.hh"
else:
    fname = "config.hh"



# The following dict of dicts define a reasonable set of defaults
# for particular trajectory/background/diffusion types.
# These defaults are intended to be more useful
# than the fallthrough defaults. This information makes this script
# the location in the code where these default values are maintained.
#
# To aid both the maintainer/developer and user,
# only the *relevant* parameters are included
# in each list.
#
# todo documentation for each parameter


reasonable_defaults = {
    'background':{
        'CylindricalObstacle': {
        },
        'Dipole': {
        },
        'Discontinuity': {
        },
        'MagnetizedCylinder': {
        },
        'MagnetizedSphere': {
        },
        'Server': {
        },
        'ServerBATL': {
        },
        'ServerCartesian': {
        },
        'Shock': {
        },
        'SmoothDiscontinuity': {
        },
        'SmoothShock': {
        },
        'SolarWind': {
            'solarwind_current_sheet': 'disabled',
            'solarwind_sectored_region': 'nowhere',
            'solarwind_polar_correction': 'none',
            'solarwind_speed_latitude_profile': 'constant',
            'solarwind_termshock_speed_exponent': 'square',
        },
        'SolarWindTermShock': {
        },
        'SphericalObstacle': {
        },
        'Uniform': {
        },
        'VLISMBochum': {
            'VLISMBochum_mod_type': 'scaled',
            'VLISMBochum_mod_rpos': 'scale_rel_zero',
        },
        'Waves': {
        },
    },
    'trajectory': {
        'Fieldline': {
            # nothing at this time.
        },
        'Focused': {
            'pperp_method': "scheme",
            'use_B_drifts': "none",
            'record_trajectory_segment_presize': 10000,
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
        },
        'Guiding': {
            'pperp_method': "scheme",
            'record_trajectory_segment_presize': 10000,
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
        },
        'GuidingDiff': {
            'pperp_method': "scheme",
            'record_trajectory_segment_presize': 10000,
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
        },
        'GuidingScatt': {
            'pperp_method': "scheme",
            'record_trajectory_segment_presize': 10000,
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
            'split_scatt': False,
            'split_scatt_fraction': 0.5,
            'const_dmumax': "constant_dthetamax",
            'stochastic_method': "Euler",
            'cfl_pitchangle': 0.5,
        },
        'GuidingDiffScatt': {
            'pperp_method': "scheme",
            'record_trajectory_segment_presize': 10000,
            'cfl_advection': 0.5,
            'drift_safety': 0.5,
            'mirror_threshold': 10,
            'split_scatt': False,
            'split_scatt_fraction': 0.5,
            'const_dmumax': "constant_dthetamax",
            'stochastic_method_mu': "Euler",
            'stochastic_method_perp': "Euler",
            'cfl_pitchangle': 0.5,
        },
        'Lorentz': {
            'pperp_method': "scheme",
            'record_trajectory_segment_presize': 100000,
            'cfl_advection': 0.1,
            'drift_safety': 0.5,
            'mirror_threshold': 300,
            'steps_per_orbit': 100,
        },
        'Parker': {
            'pperp_method': "scheme",
            'record_trajectory_segment_presize': 100000,
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
    'diffusion': {
        'IsotropicConstant': {
            'D0': 1234,
        },
        'QLTConstant': {
            'A2A': 1234,
            'l_max': 1234,
            'k_min': 1234,
            'ps_index': 1234,
            'ps_minus': 1234,
        },
        'WNLTConstant': {
            'use_qlt_scatt': False,
            'A2A': 1234,
            'l_max': 1234,
            'k_min': 1234,
            'ps_index': 1234,
            'ps_minus': 1234,
            'A2T': 1234,
            'A2L': 1234,
            'ps_plus': 1234,
        },
        'WNLTRampVLISM': {
            'A2A': 1234,
            'l_max': 1234,
            'k_min': 1234,
            'ps_index': 1234,
            'ps_minus': 1234,
            'l_max_HP': 1234,
            'dl_max': 1234,
            'nose_z_nose': 1234,
            'nose_z_sheath': 1234,
            'nose_dz': 1234,
            'kappa0': 1234,
            'kappa_ratio': 1234,
        },
        'ParaConstant': {
            'D0': 1234,
        },
        'PerpConstant': {
            'D0': 1234,
        },
        'FullConstant': {
            'Dperp': 1234,
            'Dpara': 1234,
        },
        'FlowMomentumPowerLaw': {
            'kappa0': 1234,
            'kappa_ratio': 1234,
            'U0': 1234,
            'p0': 1234,
            'pow_law_U': 1234,
            'pow_law_p': 1234,
        },
        'KineticEnergyRadialDistancePowerLaw': {
            'kappa0': 1234,
            'kappa_ratio': 1234,
            'T0': 1234,
            'r0': 1234,
            'pow_law_T': 1234,
            'pow_law_r': 1234,
        },
        'RigidityMagneticFieldPowerLaw': {
            'kappa_ratio': 1234,
            'lam0': 1234,
            'R0': 1234,
            'B0': 1234,
            'pow_law_R': 1234,
            'pow_law_B': 1234,
        },
        'StraussEtAl2013': {
            'R0': 1234,
            'B0': 1234,
            'LISM_idx': 1234,
            'LISM_ind': 1234,
            'kappa_ratio_inner': 1234,
            'kappa_ratio_outer': 1234,
            'lam_inner': 1234,
            'lam_outer': 1234,
        },
        'GuoEtAl2014': {
            'R0': 1234,
            'B0': 1234,
            'LISM_idx': 1234,
            'LISM_ind': 1234,
            'kappa_ratio_inner': 1234,
            'kappa_ratio_outer': 1234,
            'lam_inner': 1234,
            'lam_outer': 1234,
        },
        'PotgieterEtAl2015': {
            'R0': 1234,
            'B0': 1234,
            'LISM_idx': 1234,
            'LISM_ind': 1234,
            'kappa_ratio_inner': 1234,
            'kappa_ratio_outer': 1234,
            'kappa_inner': 1234,
            'kappa_outer': 1234,
        },
        'EmpiricalSOQLTandUNLT': {
            'R0': 1234,
            'B0': 1234,
            'lam_para': 1234,
            'lam_perp': 1234,
            'kappa_ratio_red': 1234,
            'radial_limit_perp_red': 1234,
            'solar_cycle_idx': 1234,
            'solar_cycle_effect': 1234,
            'Bmix_idx': 1234,
            'Bmix_ind': 1234,
        },
    },
}



# This dict defines a reasonable set of defaults that
# can be applied in general, i.e., it is a "fallthrough" set of defaults.
reasonable_defaults_baseline = {
    'local': False,
    'trajectory': "Guiding",
    'background': "Dipole",
    # ----- General ----- #
    'build_mode': "debug",
    'specie': "proton_core",
    'derivative_method': "numeric",
    'time_flow': "forward",
    'rk_integrator': "DormandPrince_54E",
    'num_trajectories': 1,
    'batch_size': 1,
    # ----- Background ----- #
    'num_numeric_grad_evals': 1,
    'incr_dmax_ratio': "{1,10000}",
    'server_interpolation_order': 1,
    'smooth_discontinuity_order': 4,
    'server_num_ghost_cells': 2,
    'solarwind_current_sheet': 'disabled',
    'solarwind_sectored_region': 'nowhere',
    'solarwind_polar_correction': 'none',
    'solarwind_speed_latitude_profile': 'constant',
    'solarwind_termshock_speed_exponent': 'square',
    'VLISMBochum_mod_type': 'scaled',
    'VLISMBochum_mod_rpos': 'scale_rel_zero',
    # ----- Trajectory ----- #
    'record_mag_extrema': False,
    'record_trajectory': False,
    'record_trajectory_segment_presize': 0,
    'trajectory_adv_safety_level': 2,
    'max_trajectory_steps': 1324,
    'max_time_adaptations': 1234,
    'n_max_calls': 1234,
    'cfl_advection': 1234,
    'cfl_diffusion': 1234,
    'cfl_acceleration': 1234,
    'cfl_pitchangle': 1234,
    'drift_safety': 1234,
    'mirror_threshold': 1234,
    'pperp_method': 1234,
    'use_B_drifts': 1234,
    'stochastic_method': 1234,
    'stochastic_method_mu': 1234,
    'stochastic_method_perp': 1234,
    'split_scatt': 1234,
    'split_scatt_fraction': 1234,
    'const_dmumax': 1234,
    'steps_per_orbit': 1234,
    'divk_method': 1234,
    'dlnp_max': 1234,
    # ----- Diffusion ----- #
    'use_qlt_scatt': False,
    'D0': 1234,
    'Dperp': 1234,
    'Dpara': 1234,
    'A2A': 1234,
    'l_max': 1234,
    'k_min': 1234,
    'ps_index': 1234,
    'ps_minus': 1234,
    'A2T': 1234,
    'A2L': 1234,
    'ps_plus': 1234,
    'l_max_HP': 1234,
    'dl_max': 1234,
    'nose_z_nose': 1234,
    'nose_z_sheath': 1234,
    'nose_z_dz': 1234,
    'kappa0': 1234,
    'kappa_ratio': 1234,
    'U0': 1234,
    'p0': 1234,
    'T0': 1234,
    'r0': 1234,
    'lam0': 1234,
    'R0': 1234,
    'B0': 1234,
    'pow_law_U': 1234,
    'pow_law_p': 1234,
    'pow_law_T': 1234,
    'pow_law_r': 1234,
    'pow_law_R': 1234,
    'pow_law_B': 1234,
    'LISM_idx': 1234,
    'LISM_ind': 1234,
    'kappa_ratio_inner': 1234,
    'kappa_ratio_outer': 1234,
    'lam_inner': 1234,
    'lam_outer': 1234,
    'kappa_inner': 1234,
    'kappa_outer': 1234,
    'lam_para': 1234,
    'lam_perp': 1234,
    'kappa_ratio_red': 1234,
    'radial_limit_perp_red': 1234,
    'solar_cycle_idx': 1234,
    'solar_cycle_effect': 1234,
    'Bmix_idx': 1234,
    'Bmix_ind': 1234,
}



def get(key, args, special_name = None, special_value = None):
    """
    This function implements resolutions of values based on the default lists.

    :param key: string
    :param args: argparse structure
    :param special_name: optional string
    :param special_value: optional string
    :return: value
    """
    if not key in args:
        raise KeyError(f"The option {key} is undefined.")
    if args.__dict__[key] is not None:
        return args.__dict__[key]
    elif special_name is not None:
        if special_name in reasonable_defaults[special_name]:
            if key in reasonable_defaults[special_name][special_value]:
                return reasonable_defaults[special_name][special_value][key]
            else:
                # the value is not in the list. Idiomatically, this means the option is not used by the special class.
                pass
        else:
            print(f"[get_value] Warning: {special_name} {special_value} not found in list. Using general default for key {key}.")
    return reasonable_defaults_baseline[key]

def check_defaults():
    """
    Check that there is always a fallthrough default value for every default used.
    This helps ensure that developers keep the script up to date.
    """
    for special_name in reasonable_defaults:
        for special_value in reasonable_defaults[special_name]:
            for key in reasonable_defaults[special_name][special_value]:
                if not key in reasonable_defaults_baseline:
                    print(f"[check_defaults] Warning: key {key} found in the list of {special_name} {special_value} defaults, but not found in the list of baseline defaults.")

def ratio(dbl):
    """
    Parse a Python immediate floating point expression
    to produce an initializer expression for a std::ratio.
    This makes setting default values in the lists (above) more convenient.
    A cost is incurred due to there occasionally being an
    unrecognizable std::ratio generated, due to roundoff error.
    For example, 2000.003 might become 20000030000000004/100000000000000
    instead of 2000003/1000.
    :param dbl: scalar
    :return: string
    """
    denom = 1
    numer = abs(dbl)
    sgn = 1 if dbl >= 0 else -1
    fp, ip = math_modf(numer)
    while fp != 0:
        denom *= 10
        numer *= 10
        fp, ip = math_modf(numer)
    return f"{{{int(sgn*numer)},{int(denom)}}}"


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Configuration file generator for SPECTRUM particle trajectory trace solver.",
    )
    for key in reasonable_defaults_baseline:
        parser.add_argument(
            f"--{key}",
            type=type(reasonable_defaults_baseline[key]),
            help="(todo)",
            required=False,
        )
    args = parser.parse_args()

    there_is_a_trajectory = '--trajectory' in sys.argv
    there_is_a_background = '--background' in sys.argv
    it_is_a_local_hconfig = '--local' in sys.argv

    if there_is_a_trajectory:
        print("[hconfig] found traj")
        # todo
    else:
        print("[hconfig] NO traj")

    if it_is_a_local_hconfig:
        print("[hconfig] local hconfig")
        # todo
    else:
        print("[hconfig] NOT local hconfig")


    # code template TODO
    content = f"""

#ifndef SPECTRUM_HCONFIG
#define SPECTRUM_HCONFIG



using HConfig = HConfig<
    trajectory = {get('trajectory', args)},
    background = {get('background', args)},
    etc.
>;


#endif

"""

    if debug:
        print(content)
    else:
        with open(fname, 'w') as f:
            f.write(content)
    check_defaults()

