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

# switch for debugging purposes
# debug = False
debug = True


if debug:
    fname = "config.test.hh"
else:
    fname = "config.hh"



# The following dicts define a reasonable set of defaults
# for particular trajectory/background types.
# These defaults are probably more useful
# than the fallthrough defaults.
# To aid both the maintainer/developer and user,
# only the *relevant* values are included
# in each list, and this information
# is propagated into the configuration
# via a comment. todo impl


reasonable_defaults_per_background = {
    'CylindricalObstacle': {
        # todo
    },
    'Dipole': {
        # todo
    },
    'MagnetizedCylinder': {
        # todo
    },
    # todo finish
}



reasonable_defaults_per_trajectory = {
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
        'split_scatt_fraction': "{1,2}",
        'const_dmumax': "constant_dthetamax",
        'stochastic_method': "Euler",
        'cfl_pitchangle': "{1,2}",
    },
    'GuidingDiffScatt': {
        'pperp_method': "scheme",
        'record_trajectory_segment_presize': 10000,
        'cfl_advection': 0.5,
        'drift_safety': 0.5,
        'mirror_threshold': 10,
        'split_scatt': False,
        'split_scatt_fraction': "{1,2}",
        'const_dmumax': "constant_dthetamax",
        'stochastic_method': "Euler",
        'cfl_pitchangle': "{1,2}",
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
}
# todo DONE (from src)




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
    'num_numeric_grad_evals': 0,
    'incr_dmax_ratio': "{1,10000}",
    'server_interpolation_order': 1,
    'smooth_discontinuity_order': 4,
    'server_num_ghost_cells': 2,
    # todo
    # ----- Trajectory ----- #
    'record_mag_extrema': False,
    'record_trajectory': False,
    'record_trajectory_segment_presize': 0,
    'trajectory_adv_safety_level': 2,
    # todo
}



def get(key, args, trajectory = None, background = None):
    """
    This function implements resolutions of values based on the default lists.

    :param key: string
    :param args: argparse structure
    :param trajectory: optional string
    :param background: optional string
    :return: value
    """
    if trajectory is not None and background is not None:
        print(f"[get_value] Warning: invalid inputs. Check implementation/inputs.")
    if not key in args:
        raise KeyError(f"The option {key} is undefined.")
    if args.__dict__[key] is not None:
        return args.__dict__[key]
    elif trajectory is not None:
        if trajectory in reasonable_defaults_per_trajectory:
            if key in reasonable_defaults_per_trajectory[trajectory]:
                return reasonable_defaults_per_trajectory[trajectory][key]
            else:
                # the value is not in the list. Idiomatically, this means the option is not used by the trajectory.
                pass
        else:
            print(f"[get_value] Warning: trajectory {trajectory} not found in list. Using general default for key {key}.")
    elif background is not None:
        if background in reasonable_defaults_per_background:
            if key in reasonable_defaults_per_background[background]:
                return reasonable_defaults_per_background[background][key]
            else:
                # the value is not in the list. Idiomatically, this means the option is not used by the background.
                pass
        else:
            print(f"[get_value] Warning: background {background} not found in list. Using general default for key {key}.")
    return reasonable_defaults_baseline[key]

def check_defaults():
    """
    Check that there is always a fallthrough default value for every default used.
    This helps ensure that developers keep the script up to date.
    """
    for trajectory in reasonable_defaults_per_trajectory:
        for key in reasonable_defaults_per_trajectory[trajectory]:
            if not key in reasonable_defaults_baseline:
                print(f"[check_defaults] Warning: key {key} found in the list of trajectory {trajectory } defaults, but not found in the list of baseline defaults.")
    for background in reasonable_defaults_per_background:
        for key in reasonable_defaults_per_background[background]:
            if not key in reasonable_defaults_baseline:
                print(f"[check_defaults] Warning: key {key} found in the list of background {background } defaults, but not found in the list of baseline defaults.")


if __name__ == "__main__":

    # argparse frontend
    import argparse
    import sys
    parser = argparse.ArgumentParser(
        description="Configuration file generator for SPECTRUM particle trajectory trace solver.",
    )
    for key in reasonable_defaults_baseline:
        parser.add_argument(
            f"--{key}",
            type=type(reasonable_defaults_baseline[key]),
            help="todo",
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




    # code template
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

