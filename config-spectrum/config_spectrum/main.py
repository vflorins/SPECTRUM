"""
File config.py Created by Lucius Schoenbaum 9/13/2025

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.

"""

import argparse
from ._impl.configdata import (
    ConfigData
)


def main():
    """
    Configuration routine entry point.
    This method's only responsibility is to set up the argparse arguments.
    Data initialization, file I/O, and code-generation
    are handled by `ConfigData`.
    """
    parser = argparse.ArgumentParser(
        prog="config-spectrum",
        description="Configuration file generator for SPECTRUM particle trajectory trace solver.",
    )
    parser.add_argument(
        "relative_dir",
        help="Relative path to directory where generated files are placed.",
    )
    # parser.add_argument(
    #     "-v", "--verbose",
    #     action="count",
    #     default=0,
    #     help="Increase verbosity level (-v: info -vv: debug -vvv: trace)"
    # )
    parser.add_argument(
        "-V", "--version",
        action="version",
        version = "config-spectrum 0.2.0",
    )
    parser.add_argument(
        "--mkdefaults", "-mk",
        action=argparse.BooleanOptionalAction,
        help="Whether to operate in default mode (generate source config files) or in user mode (generate a local config file).",
    )
    parser.add_argument(
        "--mkdefaults-test", "-mktest",
        action=argparse.BooleanOptionalAction,
        help="Whether to operate in user mode (generate a local config file) or in developer mode (update defaults and special types in source code).",
    )
    parser.add_argument(
        "--test", "-test",
        action=argparse.BooleanOptionalAction,
        help="Explicitly set the operating mode to test-only. When set, clearly labeled test files will be generated in the local directory only.",
    )
    parser.add_argument(
        "--out", "-o",
        help="Name of output filename stem. This stem will be appended with `.config.hh`. If not given, output filename is `config.hh`",
        required=False,
    )
    parser.add_argument(
        "--reconfig", "-re",
        action=argparse.BooleanOptionalAction,
        help="Whether to run config.py in reconfigure mode (modify an existing config file)."
    )
    parser.set_defaults(mkdefaults=False, mkdefaults_test=False, test=False, reconfig=False)
    CD = ConfigData()
    for spectrum_type in CD.parameters:
        parameters_ = CD.parameters[spectrum_type]
        for key in parameters_:
            if parameters_[key].parameter_type != type:
                k = CD.argkey(key, spectrum_type)
                flaglist = [f"--{k}"]
                # custom shortened flags for spectrum types
                if k in ['background', 'trajectory', 'diffusion']:
                    # b, t, d
                    flaglist.append(f"-{k[0]}")
                parser.add_argument(
                    *flaglist,
                    type=parameters_[key].argparse_parameter_type,
                    help=CD.parameters[spectrum_type][key].description,
                    required=False,
                )
    args = parser.parse_args()
    CD.init(args)
    if args.mkdefaults:
        CD.mkdefaults(args)
    else:
        print("[config-spectrum] Starting...")
        if args.reconfig:
            CD.gen_reconfigure(args)
        else:
            CD.gen_configure(args)
        print("[config-spectrum] Done.")


