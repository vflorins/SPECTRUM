# File config.parameters.py Created by Lucius Schoenbaum October 20, 2025

"""
This file contains a specification of all (hyper)parameters
for the SPECTRUM test particle trajectory solver.
It provides type information, possible ranges in some cases,
a docstring, and a globally-defined fallthrough default value.

***Notes for Developers***

New parameters can be added to this file. In order to do so,
be sure to make note of the interoperability requirement between
the lists here, and the template arguments defined in config.hh files
(background.config.hh, trajectory.config.hh, diffusion.config.hh)
in the C++ source. This should become clear after a moment.
Once these files are all in agreement, the main
`config.py` script can be run to generate new default values,
incorporating physical defaults defined in `config.physical.py`.
These "physical" defaults are case-specific, and do not
need to be defined until a specific case demands/suggests doing so.

In order to "install" any updates to definitions in this file, or to
physical defaults (in `config.physical.py`), you only need to run
the main `config.py` configuration script.

"""

import re


spectrum_types = ["Background", "Trajectory", "Diffusion"]


def get_special_types(spectrum_types):
    with open('common/compiletime_lists.hh', 'r') as f:
        ctl = f.read()
    special_types = {}
    for st in spectrum_types:
        m = re.search(f"enum class {st} {{(.*?)}}", ctl, flags=re.DOTALL)
        special_types[st] = [x[:-1] for x in m.group(1).split()]
    return special_types


def update_special_types_source(special_types, test_only):
    for st in spectrum_types:
        srcname = f'src/{st.lower()}.hh'
        with open(srcname, 'r') as f:
            content = f.read()
        m = re.search(f"^(.*?)Fields<(.*?)>;(.*?)$", content, flags=re.DOTALL)
        newlist = ""
        for special_type in special_types[st]:
            newlist += f"{st}{special_type}<HConfig>,\n"
        content = m.group(1) + "Fields<\nFConfig<>,\n" + newlist[:-2] + "\n>;" + m.group(3)
        if test_only:
            with open(f"CONFIG.{st.lower()}.TEST.hh", 'w') as f:
                f.write(content)
        else:
            with open(srcname, 'w') as f:
                f.write(content)


special_types = get_special_types(spectrum_types)


class ParameterInfo:
    """
    A data structure setting up the essential
    information needed to define a default
    config class. Provides documentation
    that can be placed where would be appropriate.

    Arguments:

        name (string):
            parameter name
        description (string):
            parameter description
        parameter_type (type or string):
            The parameter type. If string, then the
            type refers to an enum, class, or namespace.
        possible_values (optional list of string or int or float):
            A list of possible values, if discrete
        secular: If True, the parameter not configurable,
            but instead is hard-coded in the Config data structure.
    """

    def __init__(self, name, description, parameter_type, possible_values = None, secular = False):
        self.name = name
        self.description = description
        self.possible_values = possible_values
        self.parameter_type = parameter_type
        # kludge:
        self.argparse_parameter_type = str if isinstance(parameter_type, str) else parameter_type
        self.secular = secular

    def str(self, cpp_comment = True):
        """
        A plain text string summarizing the parameter.
        Formatted as C++ comment if `cpp_comment`.
        """
        out = "//! " if cpp_comment else ""
        nl = "\n// " if cpp_comment else "\n"
        out += f"Name: {self.name}"
        out += f"{nl}Description: {self.description}"
        if self.possible_values:
            out += f"{nl}Options: " + " | ".join(self.possible_values)
        if cpp_comment and self.secular:
            out += f"{nl}Secular (do not modify)"
        elif not cpp_comment:
            out += f"{nl}Secular: True"
        return out




parameters_general = {
    'specieid': ParameterInfo(
        name = 'specieid',
        description = "The specie (particle species) used for simulation",
        possible_values=["proton_core", "electron_core"], # todo
        parameter_type = "SpecieId",
    ),
    'build_mode': ParameterInfo(
        name = 'build_mode',
        description = "The build mode used for simulation",
        possible_values=["release", "debug"],
        parameter_type="BuildMode",
    ),
    'num_trajectories': ParameterInfo(
        name = 'num_trajectories',
        description = "The number of trajectories solved during simulation",
        parameter_type = int,
    ),
    'batch_size': ParameterInfo(
        name = 'batch_size',
        description = "The number of trajectories in a single batch during simulation",
        parameter_type = int,
    ),
    'max_trajectories_per_worker': ParameterInfo(
        name = 'max_trajectories_per_worker',
        description = "The number of trajectories ordinarily assigned to a worker",
        parameter_type = int,
    ),
}


parameters_background = {
    'background': ParameterInfo(
        name = 'background',
        description = "The background type used for simulation",
        possible_values=special_types['Background'],
        parameter_type="Config::Background",
    ),
    'derivative_method': ParameterInfo(
        name = 'derivative_method',
        description = "The method used to evaluate derivatives of spatially located field quantities. "
                      "Usually, the only information checked is whether or not this value is numeric.",
        possible_values=["numeric", "nonnumeric", "analytic", "datadefined"],
        parameter_type="BackgroundOptions::DerivativeMethod",
    ),
    'num_numeric_grad_evals': ParameterInfo(
        name = 'num_numeric_grad_evals',
        description = "The number of derivative evaluations applied in the derivative-averaging method",
        parameter_type = int,
    ),
    ########################################
    # DataBackground Only
    'n_servers_per_node': ParameterInfo(
        name = 'n_servers_per_node',
        description = "The number servers in each node",
        parameter_type = int,
    ),
    'allow_server_worker': ParameterInfo(
        name = 'servers_are_workers',
        description = "Whether servers are workers. Servers only exist if the background is a data-serving background.",
        parameter_type = bool,
    ),
    'file_name_pattern': ParameterInfo(
        name = 'file_name_pattern',
        description = "File name pattern used for stored data.",
        parameter_type = str,
    ),
    #
    ########################################
    'dmax0': ParameterInfo(
        name = 'dmax0',
        description = "baseline simulation-wide dmax value",
        parameter_type = float,
    ),
    'dmax_fraction': ParameterInfo(
        name = 'dmax_fraction',
        # todo background-specific documentation?
        description = "Fraction of dmax in a reduced region (background dependent)",
        parameter_type = float,
    ),
    'incr_dmax_ratio': ParameterInfo(
        name = 'incr_dmax_ratio',
        description = "What fraction of _dmax to use to calculate the field increment",
        parameter_type = float,
    ),
    'r0': ParameterInfo(
        name = 'r0',
        description = "value (todo: description)",
        parameter_type = "GeoVector",
    ),
    'u0': ParameterInfo(
        name = 'u0',
        description = "u0",
        parameter_type = "GeoVector",
    ),
    'B0': ParameterInfo(
        name = 'B0',
        description = "B0",
        parameter_type = "GeoVector",
    ),
    # 'M': ParameterInfo(
    #     name = 'M',
    #     description = "Moment",
    #     parameter_type = "GeoVector",
    # ),
    'axis': ParameterInfo(
        name = 'axis',
        description = "special field axis",
        parameter_type = "GeoVector",
    ),
    'radius': ParameterInfo(
        name = 'radius',
        description = "radial parameter",
        parameter_type = float,
    ),
    # todo name is obscure
    'n_discont': ParameterInfo(
        name = 'n_discont',
        description = "discontinuity normal",
        parameter_type = "GeoVector",
    ),
    # todo name is obscure
    'v_discont': ParameterInfo(
        name = 'v_discont',
        description = "discontinuity velocity",
        parameter_type = float,
    ),
    'compression': ParameterInfo(
        name = 'compression ratio',
        description = "discontinuity velocity",
        parameter_type = float,
    ),
    'u_down': ParameterInfo(
        name = 'u_down',
        description = "downstream velocity parameter",
        parameter_type = "GeoVector",
    ),
    'B_down': ParameterInfo(
        name = 'B_down',
        description = "downstream magnetic field parameter",
        parameter_type = "GeoVector",
    ),
    'width_discont': ParameterInfo(
        name = 'width_discont',
        description = "discontinuity width",
        parameter_type = "GeoVector",
    ),
    'width_shock': ParameterInfo(
        name = 'width_shock',
        description = "shock width",
        parameter_type = "GeoVector",
    ),
    'eprime': ParameterInfo(
        name = 'eprime',
        description = "todo",
        parameter_type = "GeoVector",
    ),
    'w0': ParameterInfo(
        name = 'w0',
        description = "todo",
        parameter_type = "GeoVector",
    ),
    'ur0': ParameterInfo(
        name = 'ur0',
        description = "todo",
        parameter_type = "GeoVector",
    ),
    'r_ref': ParameterInfo(
        name = 'r_ref',
        description = "Reference equatorial distance",
        parameter_type = float,
    ),
    'z_nose': ParameterInfo(
        name = 'z_nose',
        description = "todo",
        parameter_type = float,
    ),
    'fsl_mns': ParameterInfo(
        name = 'fsl_mns',
        description = "todo",
        parameter_type = "GeoVector",
    ),
    'fsl_pls': ParameterInfo(
        name = 'fsl_pls',
        description = "todo",
        parameter_type = "GeoVector",
    ),
    'Omega': ParameterInfo(
        name = 'Omega',
        description = "todo",
        parameter_type = "GeoVector",
    ),
    'r_TS': ParameterInfo(
        name = 'r_TS',
        description = "todo",
        parameter_type = "GeoVector",
    ),
    'w_TS': ParameterInfo(
        name = 'w_TS',
        description = "todo",
        parameter_type = "GeoVector",
    ),
    's_TS_inv': ParameterInfo(
        name = 's_TS_inv',
        description = "todo",
        parameter_type = "GeoVector",
    ),
    'dmax_TS': ParameterInfo(
        name = 'dmax_TS',
        description = "todo",
        parameter_type = "GeoVector",
    ),
    'server_interpolation_order': ParameterInfo(
        name = 'server_interpolation_order',
        description = "server_interpolation_order",
        parameter_type = int,
    ),
    'smooth_discontinuity_order': ParameterInfo(
        name = 'smooth_discontinuity_order',
        description = "Parameter controlling smoothness of discontinuity/shock",
        possible_values=["0: not continuous", "1: differentiable", "2: twice differentiable", "3: thrice differentiable", "4, 5, >5: smooth"],
        parameter_type = int,
    ),
    'server_num_ghost_cells': ParameterInfo(
        name = 'server_num_ghost_cells',
        description = "number of ghost cells (server parameter)",
        parameter_type = int,
    ),
    # [test_parker_spiral.cc] Make sure that "SOLARWIND_CURRENT_SHEET" and SOLARWIND_POLAR_CORRECTION are (#)defined as 0 in src/background_solarwind.hh.
    'solarwind_current_sheet': ParameterInfo(
        name = 'solarwind_current_sheet',
        description = "Heliospheric current sheet",
        possible_values=["disabled", "flat", "wavy_static: wavy (Jokipii-Thomas 1981) and static", "wavy_time_dependent: wavy and time-dependent"],
        parameter_type="BackgroundOptions::CurrentSheet",
    ),
    'solarwind_sectored_region': ParameterInfo(
        name = 'solarwind_sectored_region',
        description = "Magnetic topology region",
        possible_values=["nowhere", "HCS"],
        parameter_type="BackgroundOptions::SectoredRegion",
    ),
    'solarwind_polar_correction': ParameterInfo(
        name = 'solarwind_polar_correction',
        description = "Correction to Parker Spiral, mainly for polar regions",
        possible_values=["none", "Smith_Bieber: Smith-Bieber 1991", "Zurbuchen_etal: Zurbuchen et al. 1997", "Schwadron_McComas: Schwadron-McComas 2003"],
        parameter_type="BackgroundOptions::PolarCorrection",
    ),
    'solarwind_speed_latitude_profile': ParameterInfo(
        name = 'solarwind_speed_latitude_profile',
        description = "Latitudinal profile for bulk speed",
        possible_values=["constant", "linear_step", "smooth_step"],
        parameter_type="BackgroundOptions::SpeedLatitudeProfile",
    ),
    'with_termination_shock': ParameterInfo(
        name = 'with_termination_shock',
        description = "Whether the model has a spherical termination shock feature (requires extra setup, see source/documentation)",
        parameter_type=bool,
    ),
    'termshock_speed_exponent': ParameterInfo(
        name = 'termshock_speed_exponent',
        description = "Integer exponent of decrease of solar wind speed beyond the termination shock",
        possible_values=["zero", "one", "square", "cube"],
        parameter_type="BackgroundOptions::TermShockSpeedExponent",
    ),
    'mod_type': ParameterInfo(
        name = 'mod_type',
        description = "What function to use within 'get_ampfactor'",
        possible_values=["none", "zero", "constant", "scaled"],
        parameter_type="BackgroundOptions::ModType",
    ),
    'mod_rpos': ParameterInfo(
        name = 'mod_rpos',
        description = "Whether to scale relative to s=0 or s=+inf",
        possible_values=["scale_rel_zero", "scale_rel_inf"],
        parameter_type="BackgroundOptions::ModRPos",
    ),
    'tanh_width_factor': ParameterInfo(
        name = 'tanh_width_factor',
        description = "Scaling factor to better match discontinuity width when using smooth discontinuity (tanh)",
        parameter_type=float,
    ),
}


parameters_trajectory = {
    'trajectory': ParameterInfo(
        name = 'trajectory',
        description = "The trajectory type used for simulation",
        possible_values=special_types['Trajectory'],
        parameter_type="Config::Trajectory",
    ),
    'Coordinates': ParameterInfo(
        name = 'Coordinates',
        description = "The coordinates of the trajectory during the simulation.",
        parameter_type = type,
    ),
    'Fields': ParameterInfo(
        name = 'Fields',
        description = "The fields computed in the local environment of the trajectory during the simulation.",
        parameter_type = type,
    ),
    'FieldlineField_t': ParameterInfo(
        name = 'FieldlineField_t',
        description = "The tracked field for a Fieldline trajectory class.",
        parameter_type = type,
    ),
    'RecordCoordinates': ParameterInfo(
        name = 'RecordCoordinates',
        description = "The coordinates (and format spec) used to record the trajectory progress.",
        parameter_type = type,
    ),
    'time_flow': ParameterInfo(
        name = 'time_flow',
        description = "The time flow direction used for simulation",
        possible_values=["forward", "backward"],
        parameter_type="TrajectoryOptions::TimeFlow",
    ),
    'rk_integrator': ParameterInfo(
        name = 'rk_integrator',
        description = "The discretization scheme used for simulation, a Runge-Kutta scheme",
        possible_values=["Euler1E", "Euler1I", "Midpoint_2I", "Midpoint_2E", "RungeKutta_4E", "Kutta38_4E", "RungeKuttaFehlberg_54E", "CashKarp_54E", "DormandPrince_54E", "RungeKuttaFehlberg65E", "RungeKuttaFehlberg76E", "RungeKuttaFehlberg_87E", "...others (see source)"],
        parameter_type="RKIntegrator",
    ),
    'record_mag_extrema': ParameterInfo(
        name = 'record_mag_extrema',
        description = "Whether to record magnetic field extrema",
        parameter_type = bool,
    ),
    'record_trajectory': ParameterInfo(
        name = 'record_trajectory',
        description = "Whether to record segments",
        parameter_type = bool,
    ),
    'record_trajectory_segment_presize': ParameterInfo(
        name = 'record_trajectory_segment_presize',
        description = "Default initial size of trajectory segment record",
        parameter_type = int,
    ),
    'advance_safety_level': ParameterInfo(
        name = 'advance_safety_level',
        description = "Trajectory advance routine safety level",
        possible_values=["low: no checks", "medium: check dt only", "high: check dt, number of segments, and time adaptations per step"],
        parameter_type="TrajectoryOptions::SafetyLevel"
    ),
    'max_trajectory_steps': ParameterInfo(
        name = 'max_trajectory_steps',
        description = "Largest length for single trajectory",
        parameter_type = int,
    ),
    'max_time_adaptations': ParameterInfo(
        name = 'max_time_adaptations',
        description = "Largest number of time step adaptations for a single time step",
        parameter_type = int,
    ),
    'n_max_calls': ParameterInfo(
        name = 'n_max_calls',
        description = "Upper limit on the number of steps in debug mode",
        possible_values=["-1: unlimited", "n ≥ 0: limited, upper bound n"],
        parameter_type = int,
    ),
    'cfl_advection': ParameterInfo(
        name = 'cfl_advection',
        description = "CFL condition for advection",
        parameter_type=float,
    ),
    'cfl_diffusion': ParameterInfo(
        name = 'cfl_diffusion',
        description = "CFL condition for diffusion",
        parameter_type=float,
    ),
    'cfl_acceleration': ParameterInfo(
        name = 'cfl_acceleration',
        description = "CFL condition for acceleration",
        parameter_type=float,
    ),
    'cfl_pitchangle': ParameterInfo(
        name = 'cfl_pitchangle',
        description = "CFL condition for pitch angle scattering",
        parameter_type=float,
    ),
    'drift_safety': ParameterInfo(
        name = 'drift_safety',
        description = "Safety factor for drift-based time step (to modify \"drift_vel\" with a small fraction of the particle's velocity)",
        parameter_type=float,
    ),
    'mirror_threshold': ParameterInfo(
        name = 'mirror_threshold',
        description = "How many time steps to allow before recording a mirror event",
        parameter_type = int,
    ),
    'pperp_method': ParameterInfo(
        name = 'pperp_method',
        description = "Switch controlling how to calculate mu. updating mu according to the scheme does not guarantee conservation of magnetic moment, but can be used with non-adiabatic terms.",
        possible_values=["moment_cons: compute mu from magnetic moment conservation", "scheme: update according to scheme"],
        parameter_type="TrajectoryOptions::PPerpMethod",
    ),
    'use_B_drifts': ParameterInfo(
        name = 'use_B_drifts',
        description = "Flag to use gradient and curvature drifts in drift velocity calculation",
        parameter_type="TrajectoryOptions::UseBDrifts",
    ),
    'stochastic_method': ParameterInfo(
        name = 'stochastic_method',
        description = "Which stochastic method to use for scattering",
        possible_values=["Euler", "Milstein", "RK2"],
        parameter_type="TrajectoryOptions::StochasticMethod",
    ),
    'stochastic_method_mu': ParameterInfo(
        name = 'stochastic_method_mu',
        description = "which stochastic method to use for pitch angle scattering",
        possible_values=["Euler", "Milstein", "RK2"],
        parameter_type="TrajectoryOptions::StochasticMethod",
    ),
    'stochastic_method_perp': ParameterInfo(
        name = 'stochastic_method_perp',
        description = "Which stochastic method to use for perpendicular diffusion",
        possible_values=["Euler", "Milstein", "RK2"],
        parameter_type="TrajectoryOptions::StochasticMethod",
    ),
    'split_scatt_fraction': ParameterInfo(
        name = 'split_scatt_fraction',
        description = "Whether to split the diffusive advance into two (one before and one after the advection).",
        possible_values=["0.0: do not split", ">0.0: fraction of stochastic step to take before deterministic step"],
        parameter_type=float,
    ),
    'const_dmumax': ParameterInfo(
        name = 'const_dmumax',
        description = "Desired accuracy in pitch angle cosine or in pitch angle mu",
        possible_values=["constant_dtheta_max: dtheta_max = 2π/180 (deg to rad conversion factor)", "constant_dmumax: dmumax = 0.02 (desired accuracy in pitch angle cosine)"],
        parameter_type="TrajectoryOptions::ConstDmumax",
    ),
    'steps_per_orbit': ParameterInfo(
        name = 'steps_per_orbit',
        description = "Number of time steps per one orbit",
        parameter_type = int,
    ),
    'divk_method': ParameterInfo(
        name = 'divk_method',
        description = "Which method of computation to use for divK",
        possible_values=["direct: using direct central finite differences", "gradients: using background-computed gradient quantities"],
        parameter_type="TrajectoryOptions::DivkMethod",
    ),
    'dlnp_max': ParameterInfo(
        name = 'dlnp_max',
        description = "Maximum allowed fraction of momentum change per step",
        parameter_type=float,
    ),
}

parameters_diffusion = {
    'diffusion': ParameterInfo(
        name = 'diffusion',
        description = "The diffusion type used for simulation",
        possible_values=special_types['Diffusion'],
        parameter_type="Config::Diffusion",
    ),
    'Coordinates': ParameterInfo(
        name = 'Coordinates',
        description = "The coordinates where the diffusion is computed during the simulation.",
        parameter_type = type,
    ),
    'Fields': ParameterInfo(
        name = 'Fields',
        description = "The fields computed in the local environment for the diffusion computation during the simulation.",
        parameter_type = type,
    ),
    'use_qlt_scatt': ParameterInfo(
        name = 'use_qlt_scatt',
        description = "Whether to use QLT pitch angle scattering with WLNT perpendicular diffusion",
        parameter_type = bool,
    ),
    'D0': ParameterInfo(
        name = 'D0',
        description = "diffusion coefficient (if constant)",
        parameter_type=float,
    ),
    'Dperp': ParameterInfo(
        name = 'Dperp',
        description = "diffusion coefficient (if constant)",
        parameter_type=float,
    ),
    'Dpara': ParameterInfo(
        name = 'Dpara',
        description = "diffusion coefficient (if constant)",
        parameter_type=float,
    ),
    'A2A': ParameterInfo(
        name = 'A2A',
        description = "Alfven turbulence relative variance",
        parameter_type=float,
    ),
    'A2A_ref': ParameterInfo(
        name = 'A2A_ref',
        description = "Reference Alfven turbulence relative variance",
        parameter_type=float,
    ),
    'l_max': ParameterInfo(
        name = 'l_max',
        description = "Maximum turbulent lengthscale",
        parameter_type=float,
    ),
    'k_min': ParameterInfo(
        name = 'k_min',
        description = "Characteristic wavenumber",
        parameter_type=float,
    ),
    'k_min_ref': ParameterInfo(
        name = 'k_min_ref',
        description = "Reference characteristic wavenumber",
        parameter_type=float,
    ),
    'ps_index': ParameterInfo(
        name = 'ps_index',
        description = "Power spectral index",
        parameter_type=float,
    ),
    'ps_minus': ParameterInfo(
        name = 'ps_minus',
        description = "Power spectral index minus one",
        parameter_type=float,
    ),
    'A2T': ParameterInfo(
        name = 'A2T',
        description = "Transverse turbulence relative variance",
        parameter_type=float,
    ),
    'A2T_ref': ParameterInfo(
        name = 'A2T_ref',
        description = "Reference transverse turbulence relative variance",
        parameter_type=float,
    ),
    'A2L': ParameterInfo(
        name = 'A2L',
        description = "Longitudinal turbulence relative variance",
        parameter_type=float,
    ),
    'A2L_ref': ParameterInfo(
        name = 'A2L_ref',
        description = "Reference longitudinal turbulence relative variance",
        parameter_type=float,
    ),
    'ps_plus': ParameterInfo(
        name = 'ps_plus',
        description = "Power spectral index plus one",
        parameter_type=float,
    ),
    'l_max_HP': ParameterInfo(
        name = 'l_max_HP',
        description = "Largest scale of turbulence at HP (Heliopause)",
        parameter_type=float,
    ),
    'dl_max': ParameterInfo(
        name = 'dl_max',
        description = "Difference between largest scales",
        parameter_type=float,
    ),
    'z_nose': ParameterInfo(
        name = 'z_nose',
        description = "Extent of the HP in the nose direction",
        parameter_type=float,
    ),
    'fzoom': ParameterInfo(
        name = 'fzoom',
        description = "Amplification factor at the nose",
        parameter_type=float,
    ),
    'z_sheath': ParameterInfo(
        name = 'z_sheath',
        description = "Extent of the sheath in the nose direction",
        parameter_type=float,
    ),
    'nose_dz': ParameterInfo(
        name = 'nose_dz',
        description = "Difference between nose distance",
        parameter_type=float,
    ),
    'kappa0': ParameterInfo(
        name = 'kappa0',
        description = "Reference diffusion coefficient",
        parameter_type=float,
    ),
    'kappa_ratio': ParameterInfo(
        name = 'kappa_ratio',
        description = "Ratio of perpendicular to parallel diffusion",
        parameter_type=float,
    ),
    'U0': ParameterInfo(
        name = 'U0',
        description = "Flow velocity normalization factor",
        parameter_type=float,
    ),
    'p0': ParameterInfo(
        name = 'p0',
        description = "Momentum normalization factor",
        parameter_type=float,
    ),
    'T0': ParameterInfo(
        name = 'T0',
        description = "Kinetic Energy normalization factor",
        parameter_type=float,
    ),
    'r0_nf': ParameterInfo(
        name = 'r0_nf',
        description = "Radial distance normalization factor",
        parameter_type=float,
    ),
    'stream_dep_idx': ParameterInfo(
        name = 'stream_dep_idx',
        description = "Downstream dependance index",
        parameter_type=int,
    ),
    'u_upstream': ParameterInfo(
        name = 'u_upstream',
        description = "Upstream flow",
        parameter_type=float,
    ),
    'w_sh': ParameterInfo(
        name = 'w_sh',
        description = "Width of shock",
        parameter_type=float,
    ),
    's_sh': ParameterInfo(
        name = 's_sh',
        description = "Shock strength",
        parameter_type=float,
    ),
    'dn_up_ratio': ParameterInfo(
        name = 'dn_up_ratio',
        description = "Ratio of downstream to upstream value",
        parameter_type=float,
    ),
    'lam0': ParameterInfo(
        name = 'lam0',
        description = "Parallel mean free path",
        parameter_type=float,
    ),
    'R0_nf': ParameterInfo(
        name = 'R0_nf',
        description = "Rigidity normalization factor",
        parameter_type=float,
    ),
    'B0_nf': ParameterInfo(
        name = 'B0_nf',
        description = "Magnetic field normalization factor",
        parameter_type=float,
    ),
    'pow_law_U': ParameterInfo(
        name = 'pow_law_U',
        description = "Power law slope for flow velocity",
        parameter_type=float,
    ),
    'pow_law_p': ParameterInfo(
        name = 'pow_law_p',
        description = "Power law slope for momentum",
        parameter_type=float,
    ),
    'pow_law_T': ParameterInfo(
        name = 'pow_law_T',
        description = "Power law slope for kinetic energy",
        parameter_type=float,
    ),
    'pow_law_r': ParameterInfo(
        name = 'pow_law_r',
        description = "Power law slope for radial distance",
        parameter_type=float,
    ),
    'pow_law_R': ParameterInfo(
        name = 'pow_law_R',
        description = "Power law slope for rigidity",
        parameter_type=float,
    ),
    'pow_law_B': ParameterInfo(
        name = 'pow_law_B',
        description = "Power law slope for magnetic field",
        parameter_type=float,
    ),
    'LISM_idx': ParameterInfo(
        name = 'LISM_idx',
        description = "Index for LISM indicator variable",
        parameter_type=int,
    ),
    'kappa_ratio_inner': ParameterInfo(
        name = 'kappa_ratio_inner',
        description = "Ratio of perpendicular to parallel diffusion inner heliosphere",
        parameter_type=float,
    ),
    'kappa_ratio_outer': ParameterInfo(
        name = 'kappa_ratio_outer',
        description = "Ratio of perpendicular to parallel diffusion outer heliosphere",
        parameter_type=float,
    ),
    'lam_inner': ParameterInfo(
        name = 'lam_inner',
        description = "Parallel inner heliosphere mean free path",
        parameter_type=float,
    ),
    'lam_outer': ParameterInfo(
        name = 'lam_outer',
        description = "Parallel outer heliosphere mean free path",
        parameter_type=float,
    ),
    'kappa_inner': ParameterInfo(
        name = 'kappa_inner',
        description = "Parallel inner heliosphere diffusion coefficient",
        parameter_type=float,
    ),
    'kappa_outer': ParameterInfo(
        name = 'kappa_outer',
        description = "Parallel outer heliosphere diffusion coefficient",
        parameter_type=float,
    ),
    'lam_para': ParameterInfo(
        name = 'lam_para',
        description = "Parallel mean free path",
        parameter_type=float,
    ),
    'lam_perp': ParameterInfo(
        name = 'lam_perp',
        description = "Perpendicular mean free path mean free path",
        parameter_type=float,
    ),
    'kappa_ratio_red': ParameterInfo(
        name = 'kappa_ratio_red',
        description = "Reduction factor for kappa in unipolar regions",
        parameter_type=float,
    ),
    'radial_limit_perp_red': ParameterInfo(
        name = 'radial_limit_perp_red',
        description = "Radial limit to apply unipolar reduction factor",
        parameter_type=float,
    ),
    'solar_cycle_idx': ParameterInfo(
        name = 'solar_cycle_idx',
        description = "Index for solar cycle indicator variable",
        parameter_type=int,
    ),
    'solar_cycle_effect': ParameterInfo(
        name = 'solar_cycle_effect',
        description = "Solar cycle effect constant",
        parameter_type=float,
    ),
    'Bmix_idx': ParameterInfo(
        name = 'Bmix_idx',
        description = "Index for magnetic mixing indicator variable",
        parameter_type=int,
    ),
}


parameters = {
    'General': parameters_general,
    'Background': parameters_background,
    'Trajectory': parameters_trajectory,
    'Diffusion': parameters_diffusion,
}

