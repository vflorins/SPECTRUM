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
from importlib.machinery import SourceFileLoader
from os.path import (
    join as os_path_join,
    exists as os_path_exists,
)
from config_impl import (
    spectrum_path,
)


physical_defaults = SourceFileLoader("config.physical", os_path_join(spectrum_path, "config.physical.py")).load_module().physical_defaults
PM = SourceFileLoader("config.parameters", os_path_join(spectrum_path, "config.parameters.py")).load_module()


only_generate_test_files = False


def argkey(key, spectrum_type):
    """
    This is to avoid name collisions among parameter lists for different spectrum types
    when they are dumped into the common args dict by argparse.
    """
    if key in ["background", "trajectory", "diffusion"]:
        return key
    return f"{spectrum_type}_{key}"

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
            return preamble + args.__dict__[ak]
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


def gen_main_block(spectrum_type, special_type):
    content = ""
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
    return content


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Configuration file generator for SPECTRUM particle trajectory trace solver.",
    )
    parser.add_argument(
        "--mkdefaults",
        action=argparse.BooleanOptionalAction,
        help="Whether to operate in default mode (generate source config files) or in user mode (generate a local config file).",
    )
    parser.add_argument(
        "--out", "-o",
        help="Name of output filename (will be appended with .hyperconfig.hh)",
        required=False,
    )
    parser.add_argument(
        "--reconfig", "-re",
        action=argparse.BooleanOptionalAction,
        help="Whether to run config.py in reconfigure mode."
    )
    parser.set_defaults(mkdefaults=False)
    parser.set_defaults(reconfig=False)
    for spectrum_type in PM.parameters:
        parameters_ = PM.parameters[spectrum_type]
        for key in parameters_:
            if parameters_[key].parameter_type != type:
                k = argkey(key, spectrum_type)
                flaglist = [f"--{argkey(key, spectrum_type)}"]
                # custom shortened flags for spectrum types
                if k in ['background', 'trajectory', 'diffusion']:
                    flaglist.append(f"-{k[0]}")
                parser.add_argument(
                    *flaglist,
                    type=parameters_[key].argparse_parameter_type,
                    help=PM.parameters[spectrum_type][key].description,
                    required=False,
                )
    args = parser.parse_args()

    check_defaults()

    if args.mkdefaults:

        PM.update_special_types_source(PM.special_types, spectrum_path, only_generate_test_files)

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
                content += gen_main_block(spectrum_type, special_type)
                content += f"""}};

"""
            content += f"""

}}

#endif

"""
            if only_generate_test_files:
                fname = f"CONFIG.TEST.{spectrum_type.lower()}.config.TEST.hh"
            else:
                fname = os_path_join(spectrum_path, "src", f"{spectrum_type.lower()}.config.hh")
            with open(fname, 'w') as f:
                f.write(content)
            print(f"[config.py] file {fname} was updated")
    else:
        print("[config.py] Starting...")
        if not args.reconfig:
            # gen (generate) mode
            # sketch:
            # python $SPECTRUM/config.py -o config --trajectory Guiding --background Dipole --diffusion None --diffusion_use_QLTscatt true
            # todo - test definition of parameter values from command line (test all types incl. GeoVector)
            stem = args.out + "." if args.out else ""
            fname = f"{stem}hyperconfig.hh"
            if only_generate_test_files:
                fname = "TEST." + fname + ".TEST.hh"
            else:
                # a configuration file might have been edited,
                # so the policy is never to allow files to be overwritten.
                while os_path_exists(fname):
                    print(f"[config.py] Existing file {fname} will not be overwritten. If you wish to generate, use a different name, or first remove/rename the file. If you wish to reconfigure, use flag -re.")
                    fname += '.copy'
            specie = "proton_core" # todo ___________________
            B = args.background
            T = args.trajectory
            D = args.diffusion
            content = f"""// File {fname}
// {B} {T} {D}

#include "common/compiletime_lists.hh"
#include "common/vectors.hh"
#include "common/fields.hh"
#include "src/hyperconfigure.hh"

using namespace Spectrum;

constexpr auto specieid_ = SpecieId::{specie};
constexpr auto specie = Specie<specieid_>();


struct SimulationConfig1 {{
   static constexpr auto specieid = specieid_;
   static constexpr auto build_mode = BuildMode::debug;
// Whether to print the last trajectory
   static constexpr auto print_last_trajectory = true;
// Whether there is a supervisor process. This does not require there to be server processes.
   static constexpr auto supervisor = false;
// These can also be set from the command line (argv) at runtime.
   static constexpr auto num_trajectories = 1;
   static constexpr auto batch_size = 1;
   static constexpr auto max_trajectories_per_worker = 1;
}};

"""
            if B:
                content += f"""struct BackgroundConfig1{{\n"""
                content += gen_main_block("Background", B)
                content += "};\n\n"
                B = "BackgroundConfig1"
            else:
                B = "Default"
            if T:
                content += f"""struct TrajectoryConfig1{{\n"""
                content += gen_main_block("Trajectory", T)
                content += "};\n\n"
                T = "TrajectoryConfig1"
            else:
                T = "Default"
            if D:
                content += f"""struct DiffusionConfig1{{\n"""
                content += gen_main_block("Diffusion", D)
                content += "};\n\n"
                D = "DiffusionConfig1"
            else:
                D = "Default"
            content += f"""

using HConfig = HyperConfigure<
      SimulationConfig1,
      {B},
      {T},
      {D}
>;

"""
            with open(fname, 'w') as f:
                f.write(content)
            print(f"[config.py] Wrote to file {fname}")
        else:
            print("[config.py] reconfig")
            raise NotImplementedError
            # re (reconfigure) mode
            # sketch:
            # python $SPECTRUM/config.py --re main.cc --param1 1.234 --param2 2.345
        print("[config.py] Done.")






