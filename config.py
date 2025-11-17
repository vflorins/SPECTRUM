#!/usr/bin/env python3
"""
File config.py Created by Lucius Schoenbaum 9/13/2025

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.

Configuration file generator for SPECTRUM particle trajectory trace solver.

This script is developed for setting compile-time constants (parameters)
while implementing a system of varying, non-independent, contingent defaults.
It is not possible to do such a thing using compile-time language features of C++
(namely templates and the C preprocessor) as they are defined in the standard used (C++20).
Therefore a custom-built solution is necessary, and this is one possibility.

"""

import argparse
import sys
from math import modf as math_modf
from importlib.machinery import SourceFileLoader


physical_defaults = SourceFileLoader("config.physical", "config.physical.py").load_module().physical_defaults
PM = SourceFileLoader("config.parameters", "config.parameters.py").load_module()
backgrounds = PM.backgrounds
servers = PM.servers
trajectories = PM.trajectories
diffusions = PM.diffusions
parameters_general = PM.parameters_general
parameters_background = PM.parameters_background
parameters_trajectory = PM.parameters_trajectory
parameters_diffusion = PM.parameters_diffusion
parameters = PM.parameters


only_generate_test_files = False

def get(key, args, special_name = None, special_value = None):
    """
    This function implements resolutions of values based on the default lists.
    This logic is the heart of the config system setup.

    Arguments:

        key (string):
        args (argparse structure): argparse structure
        special_name (optional string): optional string
        special_value (optional string): optional string

    Returns:

        value for given key

    """
    # todo better doc'n, nomenclature for special_name, special_value
    parameterinfo = parameters[key]
    preamble = parameterinfo.parameter_type + "::" if isinstance(parameterinfo.parameter_type, str) else ""
    if not key in args:
        raise KeyError(f"[get] The option {key} is undefined.")
    if args.__dict__[key] is not None:
        return args.__dict__[key]
    elif special_name is not None:
        if special_value in physical_defaults[special_name]:
            if key in physical_defaults[special_name][special_value]:
                if parameterinfo.parameter_type == float:
                    return str(physical_defaults[special_name][special_value][key])
                elif parameterinfo.parameter_type == int:
                    return str(physical_defaults[special_name][special_value][key])
                elif parameterinfo.parameter_type == bool:
                    return str(physical_defaults[special_name][special_value][key]).lower()
                else:
                    return preamble + physical_defaults[special_name][special_value][key]
            else:
                # the value is not in the list. Idiomatically, this means the option is not used by the special class.
                pass
        else:
            print(f"[get] Warning: {special_name} {special_value} not found in physical default list. Using general default for key {key}.")
    if parameterinfo.parameter_type == float:
        return str(parameterinfo.default)
    elif parameterinfo.parameter_type == int:
        # todo this branch handles some cases where the type is not int.
        #  It works but should still be revised.
        return str(parameterinfo.default)
    elif parameterinfo.parameter_type == bool:
        return str(parameterinfo.default).lower()
    return preamble + str(parameterinfo.default)


def check_defaults():
    """
    Check that there is always a fallthrough default value for every default used.
    This helps ensure that developers keep the script up to date.
    """
    for special_name in physical_defaults:
        for special_value in physical_defaults[special_name]:
            for key in physical_defaults[special_name][special_value]:
                if not key in parameters:
                    print(f"[check_defaults] Warning: key {key} found in the list of {special_name} {special_value} physical defaults, but not found in the list of baseline defaults.")


def ratio(fp):
    """
    Parse a Python immediate floating point expression
    to produce an initializer expression for a std::ratio.
    This makes setting default values in the lists (above) more convenient.
    A cost is incurred due to there occasionally being an
    unrecognizable std::ratio generated, due to base-ten roundoff error.
    (Thanks, scientists, for using base ten.)
    For example, 2000.003 might become 20000030000000004/100000000000000
    instead of 2000003/1000.
    :param fp: scalar
    :return: string
    """
    denom = 1
    numer = abs(fp)
    sgn = 1 if fp >= 0 else -1
    fp, ip = math_modf(numer)
    while fp != 0:
        denom *= 10
        numer *= 10
        fp, ip = math_modf(numer)
    return f"std::ratio<{int(sgn*numer)},{int(denom)}>"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Configuration file generator for SPECTRUM particle trajectory trace solver.",
    )
    parser.add_argument(
        "--mkdefaults",
        action=argparse.BooleanOptionalAction,
        help="Whether to operate in default mode (generate source config files) or in user mode (generate a local config file).",
    )
    parser.set_defaults(mkdefaults=True)
    for key in parameters:
        parser.add_argument(
            f"--{key}",
            type=parameters[key].argparse_parameter_type,
            # todo cf. config.py --help - these descriptions in parameters/physical data files?
            help="(todo)",
            required=False,
        )
    # todo check loading of values of enum types via command line
    args = parser.parse_args()

    check_defaults()


    ################################################################
    # argparse draft - wip

    found_traj = '--trajectory' in sys.argv
    found_background = '--background' in sys.argv
    found_local_hconfig = '--local' in sys.argv

    # if found_traj:
    #     print("[config] found traj")
    # else:
    #     print("[config] not found traj")

    ################################################################
    # generate defaults (mkdefaults)


    if args.mkdefaults:

        # todo automated config should modify/update background.hh, trajectory.hh, diffusion.hh (lists)

        spectrum_types = [
            "background",
            "trajectory",
            "diffusion",
        ]

        for spectrum_type in spectrum_types:
            Spectrum_type = spectrum_type[0].upper() + spectrum_type[1:]
            include = 'common/fields.hh' if spectrum_type != 'trajectory' else 'trajectory.config.fields.hh'
            content = f"""/*!
\\file {spectrum_type}.config.hh
\\brief (Hyper)parameters and config(uration) options for a SPECTRUM {spectrum_type} class
\\author Lucius Schoenbaum
\\date 09/29/2025

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

/*
 * This file is automatically generated by config.py. Do not edit this file, instead edit config.py. 
 */
 
#ifndef SPECTRUM_{spectrum_type.upper()}_CONFIG_HH
#define SPECTRUM_{spectrum_type.upper()}_CONFIG_HH

#include "common/compiletime_lists.hh"
#include "{include}"

namespace Spectrum {{

template<Config::{Spectrum_type} {spectrum_type}_, SpecieId specieid_>
struct {Spectrum_type}Config;
"""
            special_types = backgrounds if spectrum_type == 'background' else trajectories if spectrum_type == 'trajectory' else diffusions
            for special_type in special_types:
                    content += f"""
/*!
\\brief (Hyper)parameters and config(uration) options for a SPECTRUM {Spectrum_type} class
\\author Lucius Schoenbaum
\\date 09/29/2025
{Spectrum_type}: {special_type}
*/
template<SpecieId specieid_>
struct {Spectrum_type}Config<Config::{Spectrum_type}::{special_type}, specieid_> {{
   static constexpr Specie<specieid_> specie;"""
                    # > set fields/coords types.
                    # todo review/improve
                    if spectrum_type == 'trajectory':
                        content += f"""
   using Coordinates = TrajectoryCoordinates<Config::{Spectrum_type}::{special_type}, specieid_>;
   using Fields = TrajectoryFields<Config::{Spectrum_type}::{special_type}, specieid_>;
   using RecordCoordinates = Coordinates;
"""
                        if special_type == 'Fieldline':
                            # todo make configurable
                            content += "   using FieldlineField_t = Mag_t;\n"
                    elif spectrum_type == 'diffusion':
                        # todo Flum or AbsFlum in DiffusionFlowMomentumPowerLaw
                        # todo those that have indicator fields
                        content += f"""
   using Coordinates = Fields<FConfig<specieid_, CoordinateSystem::cartesian, CoordinateSystem::pitchangle>, Pos_t, Time_t, Rad_t, AbsVel_t, Mom_t>;
   using Fields = Fields<FConfig<>, Mag_t, AbsMag_t, DelMag_t, DelAbsMag_t, DotMag_t, DotAbsMag_t>;
"""
                    else:
                        content += "\n"
                    for key in physical_defaults[spectrum_type][special_type]:
                        content += parameters[key].str() + "\n"
                        # todo for double/int types, make the type explicit
                        content += f"   static constexpr auto {key} = {get(key, args, spectrum_type, special_type)};\n"
                    content += f"""}};

"""
            content += f"""

}}

#endif

"""
            if only_generate_test_files:
                with open(f"{spectrum_type}.config.TEST.hh", 'w') as f:
                    f.write(content)
            else:
                with open(f"src/{spectrum_type}.config.hh", 'w') as f:
                    f.write(content)


    ################################################################
    # generate a local config file


    else:
        print("[config.py] running in !mkdefaults mode")
        # todo generate a local script - can do this later
        raise NotImplementedError



