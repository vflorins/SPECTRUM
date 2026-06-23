


from os import (
    environ as os_environ,
)
import sys
import re
from pathlib import Path
from ..config_parameters import (
    parameters,
)
from ..config_physical import (
    physical_defaults,
)



class ConfigData:
    """
    Class to manage file I/O
    and code-generation routines
    for the configuration CLT `spectrum-config`.
    """

    def __init__(self):
        """
        Read the source/environment to obtain the basic data Config
        needs to perform virtually any operation.
        """
        spectrum_env = 'SPECTRUM'
        if spectrum_env in os_environ:
            self.spectrum_path = Path(os_environ[spectrum_env])
        else:
            print(f"Could not find spectrum path using {spectrum_env} path variable.")
            print("Could not proceed. Exiting...")
            exit()
        self.spectrum_types = ["Background", "Trajectory", "Diffusion"]
        self.spectrum_type_keys = ["background", "trajectory", "diffusion"]
        self.special_types = self.read_special_types()
        self.parameters = parameters(self.special_types)
        self.check_defaults()
        self.tgt_dir = None

    def init(self, args):
        """
        Initialize after parsing arguments
        :param args: argparse structure
        """
        self.tgt_dir = Path(args.relative_dir).resolve()
        if not self.tgt_dir.is_dir():
            self.msg(f"The directory resolved from relative path {args.relative_dir} does not exist. Tried full path: {self.tgt_dir}.")
            sys.exit()


    def read_special_types(self):
        """
        Read the source to obtain the current list of special types.
        These lists are maintained in the Config namespace.
        """
        with open(self.spectrum_path / 'common' / 'compiletime_lists.hh', 'r') as f:
            ctl = f.read()
        special_types = {}
        for st in self.spectrum_types:
            m = re.search(f"enum class {st} {{(.*?)}}", ctl, flags=re.DOTALL)
            special_types[st] = [x[:-1] for x in m.group(1).split()]
        return special_types


    def argkey(self, key, spectrum_type):
        """
        This is to avoid name collisions among parameter lists for different spectrum types
        when they are dumped into the common args dict by argparse.
        """
        if key in self.spectrum_type_keys:
            return key
        return f"{spectrum_type}_{key}"


    def value(self, key, args, spectrum_type, special_type = None):
        """
        This function implements resolutions of values based on the default lists.
        This logic is the heart of the config system setup:
        it is what enables management, lookup, and overriding of defaults.

        Arguments:

            key (string): parameter name, e.g., `dmax`
            args (argparse structure): argparse structure
            spectrum_type (string): spectrum type, or 'general'
            special_type (optional string): optional special type, if None,
                there are no specializations of the `spectrum_type` (e.g., Simulation)

        Returns:

            value for given key

        """
        parameterinfo = self.parameters[spectrum_type][key]
        preamble = parameterinfo.parameter_type + "::" if isinstance(parameterinfo.parameter_type, str) else ""
        if parameterinfo.parameter_type != type:
            ak = self.argkey(key, spectrum_type)
            if not ak in args:
                raise KeyError(f"The option {ak} is undefined.")
            if args.__dict__[ak] is not None:
                return preamble + args.__dict__[ak]
        if special_type is None:
            if key in physical_defaults[spectrum_type]:
                if parameterinfo.parameter_type == float:
                    return str(physical_defaults[spectrum_type][key])
                elif parameterinfo.parameter_type == int:
                    return str(physical_defaults[spectrum_type][key])
                elif parameterinfo.parameter_type == bool:
                    return str(physical_defaults[spectrum_type][key]).lower()
                elif parameterinfo.parameter_type == str:
                    return f'"{physical_defaults[spectrum_type][key]}"'
                elif parameterinfo.parameter_type == "GeoVector" or \
                        parameterinfo.parameter_type == "GeoMatrix" or \
                        parameterinfo.parameter_type == "MultiIndex":
                    return physical_defaults[spectrum_type][key]
                else:
                    return preamble + physical_defaults[spectrum_type][key]
            else:
                ValueError(f"{spectrum_type}:{key} not found in physical default list.")
                pass
        elif special_type in physical_defaults[spectrum_type]:
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
                ValueError(f"{spectrum_type}:{special_type}:{key} not found in physical default list.")
                pass
        else:
            ValueError(f"{spectrum_type}:{special_type} not found in physical default list.")



    def check_defaults(self):
        """
        Check that there is always a fallthrough default value for every default used.
        This helps ensure that developers keep the configuration script up to date.
        A soft warning is emitted if a missing default is found.
        """
        for spectrum_type in physical_defaults:
            if spectrum_type == 'Simulation':
                for key in physical_defaults[spectrum_type]:
                    if not key in self.parameters[spectrum_type]:
                        print(f"[check_defaults] Warning: key {key} found in the list of {spectrum_type} physical defaults, but not found in the list of {spectrum_type} parameters.")
            else:
                for special_type in physical_defaults[spectrum_type]:
                    for key in physical_defaults[spectrum_type][special_type]:
                        if not key in self.parameters[spectrum_type]:
                            print(f"[check_defaults] Warning: key {key} found in the list of {spectrum_type}:{special_type} physical defaults, but not found in the list of {spectrum_type} parameters.")



    def gen_main_block(self, args, spectrum_type, special_type = None):
        """
        Generate text of main code block of a Configuration data structure.
        This is the routine used to inject defaults in `config_physical.py`
        into the SPECTRUM source when `mkdefaults` is called.

        args (argparse structure): argparse structure
        spectrum_type (string): spectrum type, or 'general'
        special_type (optional string): optional special type, if None,
            there are no specializations of the `spectrum_type` (e.g., Simulation)
        :return: string
        """
        content = ""
        keydict = physical_defaults[spectrum_type][special_type] if special_type is not None else physical_defaults[spectrum_type]
        keydict_print = spectrum_type + ": " + special_type if special_type is not None else spectrum_type
        for key in keydict:
            try:
                parameterinfo = self.parameters[spectrum_type][key]
            except:
                raise ValueError(f"Failed to find parameter {key} in the full parameter list for type {spectrum_type}.")
            content += parameterinfo.str() + "\n"
            try:
                value_definition = f"{key} = {self.value(key, args, spectrum_type, special_type)}"
            except:
                raise ValueError(f"Failed while trying to set config for parameter {key} of {keydict_print}.")
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




    def update_special_types_source(self, args):
        """
        Modify the source to incorporate new special types that may have been added
        since the last config run. This modifies the list in {special_type}.hh.
        """
        for st in self.spectrum_types:
            srcname = self.spectrum_path / 'src' / f'{st.lower()}.hh'
            with open(srcname, 'r') as f:
                content = f.read()
            m = re.search(f"^(.*?)Fields<(.*?)>;(.*?)$", content, flags=re.DOTALL)
            newlist = ""
            for special_type in self.special_types[st]:
                print(f"update_special_types_source {special_type}")
                newlist += f"{st}{special_type}<HConfig>,\n"
            content = m.group(1) + "Fields<\nFConfig<>,\n" + newlist[:-2] + "\n>;" + m.group(3)
            if args.test:
                fname = self.tgt_dir / f"CONFIG.TEST.{st.lower()}.TEST.hh"
            else:
                fname = srcname
            with open(fname, 'w') as f:
                f.write(content)
            self.msg(f"file {fname} was updated")



    @staticmethod
    def msg(x):
        print(f"[config-spectrum] {x}")





    #### CODE GENERATORS ####

    ############################

    ############################

    ############################

    ############################











    def mkdefaults(self, args):
        """
        Update the SPECTRUM source with defaults stored in
        file `config_physical.py`.

        """
        from datetime import datetime
        self.update_special_types_source(args)

        for spectrum_type in self.spectrum_types:
            content = f"""/*!
\\file {spectrum_type.lower()}.config.hh
\\brief (Hyper)parameters and config(uration) options for a SPECTRUM {spectrum_type} class
\\date {datetime.now().strftime("%m/%d/%Y")}

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

/*
* This file is automatically generated by `config-spectrum`. 
 * Do not edit this file, instead edit the parameter and default input files 
 * in `config-spectrum`. 
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
            special_types_ = self.special_types[spectrum_type]
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
                content += self.gen_main_block(args, spectrum_type, special_type)
                content += f"""}};

    """
            content += f"""
    
}}

#endif

"""
            if args.test:
                fname = self.tgt_dir / f"CONFIG.TEST.{spectrum_type.lower()}.config.TEST.hh"
            else:
                fname = self.spectrum_path / "src" / f"{spectrum_type.lower()}.config.hh"
            with open(fname, 'w') as f:
                f.write(content)
            self.msg(f"file {fname} was updated")




    def gen_configure(self, args):
        """
        Routine to execute configuration script in "configure" mode.
        sketch:
        ```
        spectrum-config -o <fname-stem> --trajectory Guiding --background Dipole --diffusion None --diffusion_use_QLTscatt true
        ```

        :param args: argparse args structure
        """
        # todo - test definition of parameter values from command line (test all types incl. GeoVector)
        stem = args.out + "." if args.out else ""
        if args.test:
            fname = self.tgt_dir / f"TEST.{stem}config.TEST.hh"
        else:
            fname = self.tgt_dir / f"{stem}config.hh"
            # a configuration file might have been edited,
            # so the policy is never to allow files to be overwritten.
            while fname.exists():
                self.msg(f"Existing file {fname} will not be overwritten. If you wish to generate a file with this exact name, first remove/rename the file. If you wish to reconfigure, use flag -re.")
                fname = fname.with_name(fname.stem + ".copy.hh")
        specie = self.value("specieid", args, 'Simulation')
        B = args.background
        T = args.trajectory
        D = args.diffusion
        content = f"""// File {fname.name}
// Automatically generated by config-spectrum
// Background/Trajectory/Diffusion:
// {B}/{T}/{D}

#include "common/compiletime_lists.hh"
#include "common/vectors.hh"
#include "common/fields.hh"
#include "src/hyperconfigure.hh"

using namespace Spectrum;

constexpr auto specieid_ = {specie};
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
            content += self.gen_main_block(args, "Background", B)
            content += "};\n\n"
            B = "BackgroundConfig1"
        else:
            B = "Default"
        if T:
            content += f"""struct TrajectoryConfig1{{\n"""
            content += self.gen_main_block(args, "Trajectory", T)
            content += "};\n\n"
            T = "TrajectoryConfig1"
        else:
            T = "Default"
        if D:
            content += f"""struct DiffusionConfig1{{\n"""
            content += self.gen_main_block(args, "Diffusion", D)
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
        self.msg(f"Wrote to file {fname}")





    def gen_reconfigure(self, args):
        """
        Routine to execute configuration script in "reconfigure" mode.
        This mode is not implemented, it may be implemented in the future,
        but this is subject to further review.
        It achieves what can also be achieved by manually editing the config file,
        so it is not needed provided there is a pause for human activity
        between runs. However, this may not be the case if the
        runs are automated by scripting, or handled by an automated optimizer.

        sketch:
        ```
        spectrum-config --re main.cc --param1 1.234 --param2 2.345
        ```
        :param args: argparse args structure
        """
        self.msg("Reconfigure mode")
        raise NotImplementedError


