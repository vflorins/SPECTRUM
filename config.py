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
                    return ratio(physical_defaults[special_name][special_value][key])
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
        return ratio(parameterinfo.default)
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
        help="Whether to operate in default mode (generate default headers) or in user mode (generate a local config file).",
    )
    parser.set_defaults(mkdefaults=True)
    for key in parameters:
        parser.add_argument(
            f"--{key}",
            type=type(parameters[key]),
            # todo cf. config.py --help - these descriptions in parameters/physical data files?
            help="(todo)",
            required=False,
        )
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
        # todo here, generate default files (clobber old ones -- thus, ground truth lives here)

        spectrum_types = [
            # "background",
            # "trajectory",
            "diffusion",
        ]

        for spectrum_type in spectrum_types:
            content = f"""/*!
\\file {spectrum_type}.config.default.hh
\\brief (Hyper)parameters and config(uration) options for a SPECTRUM {spectrum_type} class
\\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

/*
 * This file is automatically generated by config.py. Do not edit this file, instead edit config.py. 
 */
 
#ifndef SPECTRUM_{spectrum_type.upper()}_CONFIG_DEFAULT_HH
#define SPECTRUM_{spectrum_type.upper()}_CONFIG_DEFAULT_HH

#include "{spectrum_type}.config.hh"

namespace Spectrum {{

"""
            if spectrum_type == "background":

                for background in backgrounds:
                    if background in servers:
                        content += f"template <typename HConfig, typename ServerFront>\nclass Background{background};\n"
                    else:
                        content += f"template <typename HConfig>\nclass Background{background};\n"
                content += "\ntemplate <typename Background>\nclass BackgroundDefault;\n\n"
                for background in backgrounds:
                    if background in servers:
                        server1 = ", typename ServerFront"
                        server2 = ", ServerFront"
                    else:
                        server1 = ""
                        server2 = ""
                    content += f"""
template <typename HConfig{server1}>
class BackgroundDefault<Background{background}<HConfig{server2}>> {{
   using type = BackgroundConfig<
      HConfig::specieid,\n"""
                    for key in parameters_background:
                        content += f"{parameters[key].str()}\n    {get(key, args, spectrum_type, background)},\n"
                    # comma...
                    if len(parameters_background) > 0:
                        content = content[:-2] + "\n"
                    content += f"""   >;
}};

"""
            elif spectrum_type == "trajectory":
                content += "template <TrajectoryId trajectoryid, SpecieId specieid>\nclass TrajectoryDefault;\n\n"
                for trajectory in trajectories:
                    content += f"""
template <SpecieId specieid>
class TrajectoryDefault<TrajectoryId::{trajectory}, specieid> {{
   using type = TrajectoryConfig<
                specieid,
                TrajectoryId::{trajectory},
"""
                    for key in parameters_trajectory:
                        content += f"{parameters[key].str()}\n    {get(key, args, spectrum_type, trajectory)},\n"
                    # comma...
                    if len(parameters_trajectory) > 0:
                        content = content[:-2] + "\n"
                    content += f"""   >;
}};

"""
            elif spectrum_type == "diffusion":
                for diffusion in diffusions:
                    content += f"template <typename HConfig>\nclass Diffusion{diffusion};\n"
                content += "\ntemplate <typename Diffusion>\nclass DiffusionDefault;\n\n"
                for diffusion in diffusions:
                    content += f"""
template <typename HConfig>
class DiffusionDefault<Diffusion{diffusion}<HConfig>> {{
   using type = DiffusionConfig<
                HConfig::specieid,
"""
                    for key in parameters_diffusion:
                        content += f"{parameters[key].str()}\n    {get(key, args, spectrum_type, diffusion)},\n"
                    # comma...
                    if len(parameters_diffusion) > 0:
                        content = content[:-2] + "\n"
                    content += f"""   >;
}};

"""
            content += f"""

}}

#endif

"""
            if only_generate_test_files:
                with open(f"{spectrum_type}.config.default.TEST.hh", 'w') as f:
                    f.write(content)
            else:
                with open(f"src/{spectrum_type}.config.default.hh", 'w') as f:
                    f.write(content)


    ################################################################
    # generate a local config file


    else:
        print("[config.py] running in !mkdefaults mode")
        # todo generate a local script - can do this later
        raise NotImplementedError



