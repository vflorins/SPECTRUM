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



backgrounds = [
    "CylindricalObstacle",
    "Dipole",
    "Discontinuity",
    "MagnetizedCylinder",
    "MagnetizedSphere",
    "Server",
    "ServerBATL",
    "ServerCartesian",
    "Shock",
    "SmoothDiscontinuity",
    "SmoothShock",
    "SolarWind",
    "SolarWindTermShock",
    "SphericalObstacle",
    "Uniform",
    "VLISMBochum",
    "Waves",
]

# todo documentation
servers = {
    "Server": 0,
    "ServerBATL": 0,
    "ServerCartesian": 0,
}

trajectories = [
    "Fieldline",
    "Focused",
    "Guiding",
    "Lorentz",
    "Parker",
]

diffusions = [
    "None",
    "IsotropicConstant",
    "QLTConstant",
    "WNLTConstant",
    "WNLTRampVLISM",
    "ParaConstant",
    "PerpConstant",
    "FullConstant",
    "FlowMomentumPowerLaw",
    "KineticEnergyRadialDistancePowerLaw",
    "RigidityMagneticFieldPowerLaw",
    "StraussEtAl2013",
    "GuoEtAl2014",
    "PotgieterEtAl2015",
    "EmpiricalSOQLTandUNLT",
]



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
        default (string or int or float):
            a default value
        possible_values (list of string or int or float):
            A list of possible values, if discrete
        parameter_type (type or string):
            The parameter type. If string, then the
            type refers to an enum, class, or namespace.

    """

    def __init__(self, name, description, default, possible_values = None, parameter_type = int):
        self.name = name
        self.description = description
        self.possible_values = possible_values
        self.default = default
        self.parameter_type = parameter_type

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
        out += f"{nl}Default: {self.default}"
        return out




parameters_general = {
    'background': ParameterInfo(
        name = 'background',
        description = "The background type used for simulation",
        default = "Dipole",
        possible_values=backgrounds,
    ),
    'trajectory': ParameterInfo(
        name = 'trajectory',
        description = "The trajectory type used for simulation",
        default = "Lorentz",
        possible_values=trajectories,
    ),
    'diffusion': ParameterInfo(
        name = 'diffusion',
        description = "The diffusion type used for simulation",
        default = "Uniform",
        possible_values=diffusions,
    ),
    'specieid': ParameterInfo(
        name = 'specieid',
        description = "The specie (particle species) used for simulation",
        default = "proton_core",
        possible_values=["proton_core", "electron_core"], # todo
    ),
    'trajectoryid': ParameterInfo(
        name = 'trajectoryid',
        description = "The trajectory family used for simulation",
        default = "Lorentz",
        possible_values=["Guiding", "Lorentz"], # todo
    ),
    'build_mode': ParameterInfo(
        name = 'build_mode',
        description = "The build mode used for simulation",
        default = "release",
        possible_values=["release", "debug"],
        parameter_type="BuildMode",
    ),
    'num_trajectories': ParameterInfo(
        name = 'num_trajectories',
        description = "The number of trajectories solved during simulation",
        default = 1,
    ),
    'batch_size': ParameterInfo(
        name = 'batch_size',
        description = "The number of trajectories in a single batch during simulation",
        default = 1,
    ),
}


# entries must match as a list with lists in background.config.hh
parameters_background = {
    'derivative_method': ParameterInfo(
        name = 'derivative_method',
        description = "The method used to evaluate derivatives of spatially located field quantities",
        default = "numeric",
        possible_values=["analytic", "numeric"],
        parameter_type="BackgroundOptions::DerivativeMethod",
    ),
    'num_numeric_grad_evals': ParameterInfo(
        name = 'num_numeric_grad_evals',
        description = "The number of derivative evaluations applied in the derivative-averaging method",
        default = 1,
    ),
    'incr_dmax_ratio': ParameterInfo(
        name = 'incr_dmax_ratio',
        description = "What fraction of _dmax to use to calculate the field increment",
        default = 0.0001,
        parameter_type = float,
    ),
    'server_interpolation_order': ParameterInfo(
        name = 'server_interpolation_order',
        description = "server_interpolation_order",
        default = 1,
    ),
    'smooth_discontinuity_order': ParameterInfo(
        name = 'smooth_discontinuity_order',
        description = "Parameter controlling smoothness of discontinuity/shock",
        default = 4,
        possible_values=["0: not continuous", "1: differentiable", "2: twice differentiable", "3: thrice differentiable", "4, 5, >5: smooth"],
    ),
    'server_num_ghost_cells': ParameterInfo(
        name = 'server_num_ghost_cells',
        description = "number of ghost cells (server parameter)",
        default = 2,
    ),
    # [test_parker_spiral.cc] Make sure that "SOLARWIND_CURRENT_SHEET" and SOLARWIND_POLAR_CORRECTION are (#)defined as 0 in src/background_solarwind.hh.
    'solarwind_current_sheet': ParameterInfo(
        name = 'solarwind_current_sheet',
        description = "Heliospheric current sheet",
        default = 'disabled',
        possible_values=["disabled", "flat", "wavy_static: wavy (Jokipii-Thomas 1981) and static", "wavy_time_dependent: wavy and time-dependent"],
        parameter_type="BackgroundOptions::CurrentSheet",
    ),
    'solarwind_sectored_region': ParameterInfo(
        name = 'solarwind_sectored_region',
        description = "Magnetic topology region",
        default = 'nowhere',
        possible_values=["nowhere", "HCS"],
        parameter_type="BackgroundOptions::SectoredRegion",
    ),
    'solarwind_polar_correction': ParameterInfo(
        name = 'solarwind_polar_correction',
        description = "Correction to Parker Spiral, mainly for polar regions",
        default = 'none',
        possible_values=["none", "Smith_Bieber: Smith-Bieber 1991", "Zurbuchen_etal: Zurbuchen et al. 1997", "Schwadron_McComas: Schwadron-McComas 2003"],
        parameter_type="BackgroundOptions::PolarCorrection",
    ),
    'solarwind_speed_latitude_profile': ParameterInfo(
        name = 'solarwind_speed_latitude_profile',
        description = "Latitudinal profile for bulk speed",
        default = 'constant',
        possible_values=["constant", "linear_step", "smooth_step"],
        parameter_type="BackgroundOptions::SpeedLatitudeProfile",
    ),
    'solarwind_termshock_speed_exponent': ParameterInfo(
        name = 'solarwind_termshock_speed_exponent',
        description = "Integer exponent of decrease of solar wind speed beyond the termination shock",
        default = 'square',
        possible_values=["zero", "one", "square", "cube"],
        parameter_type="BackgroundOptions::TermShockSpeedExponent",
    ),
    'mod_type': ParameterInfo(
        name = 'mod_type',
        description = "What function to use within 'get_ampfactor'",
        default = 'scaled',
        possible_values=["none", "zero", "constant", "scaled"],
        parameter_type="BackgroundOptions::ModType",
    ),
    'mod_rpos': ParameterInfo(
        name = 'mod_rpos',
        description = "Whether to scale relative to s=0 or s=+inf",
        default = 'scale_rel_zero',
        possible_values=["scale_rel_zero", "scale_rel_inf"],
        parameter_type="BackgroundOptions::ModRPos",
    ),
}


# entries must match as a list with lists in trajectory.config.hh
parameters_trajectory = {
    'TrajectoryFields': ParameterInfo(
        name = 'TrajectoryFields',
        description = "The fields needed to advance the trajectory during simulation.",
        default = "Default",
    ),
    'RecordCoordinates': ParameterInfo(
        name = 'RecordCoordinates',
        description = "The coordinates (and format spec) used to record the trajectory progress.",
        default = "Fields<FConfig<>, Pos_t, Time_t>",
    ),
    'time_flow': ParameterInfo(
        name = 'time_flow',
        description = "The time flow direction used for simulation",
        default = "forward",
        possible_values=["forward", "backward"],
        parameter_type="TrajectoryOptions::TimeFlow",
    ),
    'rk_integrator': ParameterInfo(
        name = 'rk_integrator',
        description = "The discretization scheme used for simulation, a Runge-Kutta scheme",
        default = "DormandPrince_54E",
        possible_values=["Euler1E", "Euler1I", "Midpoint_2I", "Midpoint_2E", "RungeKutta_4E", "Kutta38_4E", "RungeKuttaFehlberg_54E", "CashKarp_54E", "DormandPrince_54E", "RungeKuttaFehlberg65E", "RungeKuttaFehlberg76E", "RungeKuttaFehlberg_87E", "...others (see source)"],
        parameter_type="RKIntegrator",
    ),
    'record_mag_extrema': ParameterInfo(
        name = 'record_mag_extrema',
        description = "Whether to record magnetic field extrema",
        default = False,
        parameter_type = bool,
    ),
    'record_trajectory': ParameterInfo(
        name = 'record_trajectory',
        description = "Whether to record segments",
        default = False,
        parameter_type = bool,
    ),
    'record_trajectory_segment_presize': ParameterInfo(
        name = 'record_trajectory_segment_presize',
        description = "Default initial size of trajectory segment record",
        default = 10000,
    ),
    'trajectory_adv_safety_level': ParameterInfo(
        name = 'trajectory_adv_safety_level',
        description = "Trajectory advance safety level",
        possible_values=["low: no checks", "medium: check dt only", "high: check dt, number of segments, and time adaptations per step"],
        default = 'low',
        parameter_type="TrajectoryOptions::SafetyLevel"
    ),
    'max_trajectory_steps': ParameterInfo(
        name = 'max_trajectory_steps',
        description = "Largest length for single trajectory",
        default = 100000,
    ),
    'max_time_adaptations': ParameterInfo(
        name = 'max_time_adaptations',
        description = "Largest number of time step adaptations for a single time step",
        default = 1,
    ),
    'n_max_calls': ParameterInfo(
        name = 'n_max_calls',
        description = "Upper limit on the number of steps in debug mode",
        possible_values=["-1: unlimited", "n ≥ 0: limited, upper bound n"],
        default = -1,
    ),
    'cfl_advection': ParameterInfo(
        name = 'cfl_advection',
        description = "CFL condition for advection",
        default = 0.5,
        parameter_type=float,
    ),
    'cfl_diffusion': ParameterInfo(
        name = 'cfl_diffusion',
        description = "CFL condition for diffusion",
        default = 0.5,
        parameter_type=float,
    ),
    'cfl_acceleration': ParameterInfo(
        name = 'cfl_acceleration',
        description = "CFL condition for acceleration",
        default = 0.5,
        parameter_type=float,
    ),
    'cfl_pitchangle': ParameterInfo(
        name = 'cfl_pitchangle',
        description = "CFL condition for pitch angle scattering",
        default = 0.5,
        parameter_type=float,
    ),
    'drift_safety': ParameterInfo(
        name = 'drift_safety',
        description = "Safety factor for drift-based time step (to modify \"drift_vel\" with a small fraction of the particle's velocity)",
        default = 0.5,
        parameter_type=float,
    ),
    'mirror_threshold': ParameterInfo(
        name = 'mirror_threshold',
        description = "How many time steps to allow before recording a mirror event",
        default = 1,
    ),
    'pperp_method': ParameterInfo(
        name = 'pperp_method',
        description = "Switch controlling how to calculate mu. updating mu according to the scheme does not guarantee conservation of magnetic moment, but can be used with non-adiabatic terms.",
        possible_values=["moment_cons: compute mu from magnetic moment conservation", "scheme: update according to scheme"],
        default = "mag_moment_conservation",
        parameter_type="TrajectoryOptions::PPerpMethod",
    ),
    'use_B_drifts': ParameterInfo(
        name = 'use_B_drifts',
        description = "Flag to use gradient and curvature drifts in drift velocity calculation",
        default = "none",
        parameter_type="TrajectoryOptions::UseBDrifts",
    ),
    'stochastic_method': ParameterInfo(
        name = 'stochastic_method',
        description = "Which stochastic method to use for scattering",
        possible_values=["Euler", "Milstein", "RK2"],
        default = "Euler",
        parameter_type="TrajectoryOptions::StochasticMethod",
    ),
    'stochastic_method_mu': ParameterInfo(
        name = 'stochastic_method_mu',
        description = "which stochastic method to use for pitch angle scattering",
        possible_values=["Euler", "Milstein", "RK2"],
        default = "Euler",
        parameter_type="TrajectoryOptions::StochasticMethod",
    ),
    'stochastic_method_perp': ParameterInfo(
        name = 'stochastic_method_perp',
        description = "Which stochastic method to use for perpendicular diffusion",
        possible_values=["Euler", "Milstein", "RK2"],
        default = "Euler",
        parameter_type="TrajectoryOptions::StochasticMethod",
    ),
    'split_scatt_fraction': ParameterInfo(
        name = 'split_scatt_fraction',
        description = "Whether to split the diffusive advance into two (one before and one after the advection).",
        possible_values=["0.0: do not split", ">0.0: fraction of stochastic step to take before deterministic step"],
        default = 0.0,
        parameter_type=float,
    ),
    'const_dmumax': ParameterInfo(
        name = 'const_dmumax',
        description = "Desired accuracy in pitch angle cosine or in pitch angle mu",
        possible_values=["constant_dtheta_max: dtheta_max = 2π/180 (deg to rad conversion factor)", "constant_dmumax: dmumax = 0.02 (desired accuracy in pitch angle cosine)"],
        default = "constant_dmu_max",
        parameter_type="TrajectoryOptions::ConstDmumax",
    ),
    'steps_per_orbit': ParameterInfo(
        name = 'steps_per_orbit',
        description = "Number of time steps per one orbit",
        default = 1,
    ),
    'divk_method': ParameterInfo(
        name = 'divk_method',
        description = "Which method of computation to use for divK",
        possible_values=["direct: using direct central finite differences", "gradients: using background-computed gradient quantities"],
        default = "gradients",
        parameter_type="TrajectoryOptions::DivkMethod",
    ),
    'dlnp_max': ParameterInfo(
        name = 'dlnp_max',
        description = "Maximum allowed fraction of momentum change per step",
        default = 1.0,
        parameter_type=float,
    ),
}

# entries must match as a list with lists in diffusion.config.hh
parameters_diffusion = {
    'use_qlt_scatt': ParameterInfo(
        name = 'use_qlt_scatt',
        description = "Whether to use QLT pitch angle scattering with WLNT perpendicular diffusion",
        default = False,
        parameter_type = bool,
    ),
    'D0': ParameterInfo(
        name = 'D0',
        description = "diffusion coefficient (if constant)",
        default = 1.0,
        parameter_type=float,
    ),
    'Dperp': ParameterInfo(
        name = 'Dperp',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'Dpara': ParameterInfo(
        name = 'Dpara',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'A2A': ParameterInfo(
        name = 'A2A',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'l_max': ParameterInfo(
        name = 'l_max',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'k_min': ParameterInfo(
        name = 'k_min',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'ps_index': ParameterInfo(
        name = 'ps_index',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'ps_minus': ParameterInfo(
        name = 'ps_minus',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'A2T': ParameterInfo(
        name = 'A2T',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'A2L': ParameterInfo(
        name = 'A2L',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'ps_plus': ParameterInfo(
        name = 'ps_plus',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'l_max_HP': ParameterInfo(
        name = 'l_max_HP',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'dl_max': ParameterInfo(
        name = 'dl_max',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'nose_z_nose': ParameterInfo(
        name = 'nose_z_nose',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'nose_z_sheath': ParameterInfo(
        name = 'nose_z_sheath',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'nose_z_dz': ParameterInfo(
        name = 'nose_z_dz',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'kappa0': ParameterInfo(
        name = 'kappa0',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'kappa_ratio': ParameterInfo(
        name = 'kappa_ratio',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'U0': ParameterInfo(
        name = 'U0',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'p0': ParameterInfo(
        name = 'p0',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'T0': ParameterInfo(
        name = 'T0',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'r0': ParameterInfo(
        name = 'r0',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'lam0': ParameterInfo(
        name = 'lam0',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'R0': ParameterInfo(
        name = 'R0',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'B0': ParameterInfo(
        name = 'B0',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'pow_law_U': ParameterInfo(
        name = 'pow_law_U',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'pow_law_p': ParameterInfo(
        name = 'pow_law_p',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'pow_law_T': ParameterInfo(
        name = 'pow_law_T',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'pow_law_r': ParameterInfo(
        name = 'pow_law_r',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'pow_law_R': ParameterInfo(
        name = 'pow_law_R',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'pow_law_B': ParameterInfo(
        name = 'pow_law_B',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'LISM_idx': ParameterInfo(
        name = 'LISM_idx',
        description = "todo_add_description",
        default = 1,
        parameter_type=int,
    ),
    'LISM_ind': ParameterInfo(
        name = 'LISM_ind',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'kappa_ratio_inner': ParameterInfo(
        name = 'kappa_ratio_inner',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'kappa_ratio_outer': ParameterInfo(
        name = 'kappa_ratio_outer',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'lam_inner': ParameterInfo(
        name = 'lam_inner',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'lam_outer': ParameterInfo(
        name = 'lam_outer',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'kappa_inner': ParameterInfo(
        name = 'kappa_inner',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'kappa_outer': ParameterInfo(
        name = 'kappa_outer',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'lam_para': ParameterInfo(
        name = 'lam_para',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'lam_perp': ParameterInfo(
        name = 'lam_perp',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'kappa_ratio_red': ParameterInfo(
        name = 'kappa_ratio_red',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'radial_limit_perp_red': ParameterInfo(
        name = 'radial_limit_perp_red',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'solar_cycle_idx': ParameterInfo(
        name = 'solar_cycle_idx',
        description = "todo_add_description",
        default = 1,
        parameter_type=int,
    ),
    'solar_cycle_effect': ParameterInfo(
        name = 'solar_cycle_effect',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
    'Bmix_idx': ParameterInfo(
        name = 'Bmix_idx',
        description = "todo_add_description",
        default = 1,
        parameter_type=int,
    ),
    'Bmix_ind': ParameterInfo(
        name = 'Bmix_ind',
        description = "todo_add_description",
        default = 1.0,
        parameter_type=float,
    ),
}


parameters = parameters_general | parameters_background | parameters_trajectory | parameters_diffusion


