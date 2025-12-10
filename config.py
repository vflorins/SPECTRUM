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


only_generate_test_files = False


def argkey(key, spectrum_type):
    """
    This is to avoid name collisions among parameter lists for different spectrum types
    when they are dumped into the common args dict by argparse.
    """
    return f"{key}_{spectrum_type}"

def get(key, args, spectrum_type, special_type):
    """
    This function implements resolutions of values based on the default lists.
    This logic is the heart of the config system setup.

    Arguments:

        key (string):
        args (argparse structure): argparse structure
        spectrum_type (string): spectrum type, or 'general'
        special_type (string): special type

    Returns:

        value for given key

    """
    parameterinfo = PM.parameters[spectrum_type][key]
    preamble = parameterinfo.parameter_type + "::" if isinstance(parameterinfo.parameter_type, str) else ""
    if parameterinfo.parameter_type != type:
        ak = argkey(key, spectrum_type)
        if not ak in args:
            raise KeyError(f"[get] The option {ak} is undefined.")
        if args.__dict__[ak] is not None:
            return args.__dict__[ak]
    if special_type in physical_defaults[spectrum_type]:
        if key in physical_defaults[spectrum_type][special_type]:
            if parameterinfo.parameter_type == float:
                return str(physical_defaults[spectrum_type][special_type][key])
            elif parameterinfo.parameter_type == int:
                return str(physical_defaults[spectrum_type][special_type][key])
            elif parameterinfo.parameter_type == bool:
                return str(physical_defaults[spectrum_type][special_type][key]).lower()
            elif parameterinfo.parameter_type == str:
                return f'"{physical_defaults[spectrum_type][special_type][key]}"'
            elif parameterinfo.parameter_type == "GeoVector" or \
                parameterinfo.parameter_type == "GeoMatrix" or \
                parameterinfo.parameter_type == "MultiIndex":
                return physical_defaults[spectrum_type][special_type][key]
            else:
                return preamble + physical_defaults[spectrum_type][special_type][key]
        else:
            ValueError(f"[get] {spectrum_type}:{special_type}:{key} not found in physical default list.")
            pass
    else:
        ValueError(f"[get] {spectrum_type}:{special_type} not found in physical default list.")


def check_defaults():
    """
    Check that there is always a fallthrough default value for every default used.
    This helps ensure that developers keep the script up to date.
    """
    for spectrum_type in physical_defaults:
        for special_type in physical_defaults[spectrum_type]:
            for key in physical_defaults[spectrum_type][special_type]:
                if not key in PM.parameters[spectrum_type]:
                    print(f"[check_defaults] Warning: key {key} found in the list of {spectrum_type}:{special_type} physical defaults, but not found in the list of {spectrum_type} parameters.")


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
    for spectrum_type in PM.parameters:
        parameters_ = PM.parameters[spectrum_type]
        for key in parameters_:
            if parameters_[key].parameter_type != type:
                parser.add_argument(
                    f"--{argkey(key, spectrum_type)}",
                    type=parameters_[key].argparse_parameter_type,
                    help=PM.parameters[spectrum_type][key].description,
                    required=False,
                )
    args = parser.parse_args()

    check_defaults()

    ################################################################
    # argparse draft - wip

    # todo check loading of values of enum types via command line

    found_traj = '--trajectory' in sys.argv
    found_background = '--background' in sys.argv
    found_local_hconfig = '--local' in sys.argv

    # if found_traj:
    #     print("[config] found traj")
    # else:
    #     print("[config] not found traj")

    ################################################################

    local_config = ""

    PM.update_special_types_source(PM.special_types, only_generate_test_files)

    for spectrum_type in PM.spectrum_types:
        content = f"""/*!
\\file {spectrum_type.lower()}.config.hh
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
#include "common/fields.hh"

namespace Spectrum {{

// forward declaration of vector type
// class GeoVector;

template<Config::{spectrum_type} {spectrum_type.lower()}_, SpecieId specieid_>
struct {spectrum_type}Config;
"""
        special_types_ = PM.special_types[spectrum_type]
        for special_type in special_types_:
                content += f"""
/*!
\\brief (Hyper)parameters and config(uration) options for a SPECTRUM {spectrum_type} class
\\author Lucius Schoenbaum
\\date 09/29/2025
{spectrum_type}: {special_type}
*/
template<SpecieId specieid_>
struct {spectrum_type}Config<Config::{spectrum_type}::{special_type}, specieid_> {{
"""
                for key in physical_defaults[spectrum_type][special_type]:
                    try:
                        parameterinfo = PM.parameters[spectrum_type][key]
                    except:
                        raise ValueError(f"Failed to find parameter {key} in the full parameter list for type {spectrum_type}.")
                    content += parameterinfo.str() + "\n"
                    try:
                        value_definition = f"{key} = {get(key, args, spectrum_type, special_type)}"
                    except:
                        raise ValueError(f"Failed while trying to set config for parameter {key} of {spectrum_type}: {special_type}.")
                    if parameterinfo.parameter_type == int:
                        type_definition = "static constexpr int"
                    elif parameterinfo.parameter_type == bool:
                        type_definition = "static constexpr bool"
                    elif parameterinfo.parameter_type == float:
                        type_definition = "static constexpr double"
                    elif parameterinfo.parameter_type == str:
                        type_definition = "static constexpr std::string_view"
                    elif parameterinfo.parameter_type == type:
                        type_definition = "using"
                    elif parameterinfo.parameter_type == "GeoVector":
                        type_definition = "static constexpr GeoVector"
                    elif parameterinfo.parameter_type == "GeoMatrix":
                        type_definition = "static constexpr GeoMatrix"
                    else:
                        # todo deprecate in favor of string, like GeoVector
                        type_definition = "static constexpr auto"
                    content += f"   {type_definition} {value_definition};\n"
                content += f"""}};

"""
        content += f"""

}}

#endif

"""
        if only_generate_test_files:
            with open(f"CONFIG.{spectrum_type.lower()}.config.TEST.hh", 'w') as f:
                f.write(content)
        elif args.mkdefaults:
            with open(f"src/{spectrum_type.lower()}.config.hh", 'w') as f:
                f.write(content)
        else:
            print("[config.py] running in !mkdefaults mode")
            # todo operate on --trajectory, --diffusion, --background
            # todo push config to local_config
            # todo writing into the source should be non-default behavior
            raise NotImplementedError


# gen (generate) mode
# sketch:
# python $SPECTRUM/config.py --gen config.hh --trajectory Guiding --background Dipole --diffusion None --diffusion_coefficient 1.0 --etc.
# g++ compile it
# ./a.out

# re (reconfigure) mode
# sketch:
# python $SPECTRUM/config.py --re main.cc --param1 1.234 --param2 2.345



