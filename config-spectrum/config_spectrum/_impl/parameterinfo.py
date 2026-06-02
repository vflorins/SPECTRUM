




class ParameterInfo:
    """
    A data structure setting up the essential
    information needed to define a default
    config class. Provides documentation string,
    for placement where users will be
    able to easily have access.

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
    # todo the documentation should be class-specific in a few cases.
    #  get the documentation via a method, and redirect in the method optionally.

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




