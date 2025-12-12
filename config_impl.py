


from os import (
    environ as os_environ,
)
import re
from os.path import (
    join as os_path_join,
)

spectrum_types = ["Background", "Trajectory", "Diffusion"]
spectrum_path = os_environ['SPECTRUM']

def get_special_types(spectrum_types):
    """
    Read the source to obtain the current list of special types.
    These lists are maintained in the Config namespace.
    """
    with open(os_path_join(spectrum_path, 'common', 'compiletime_lists.hh'), 'r') as f:
        ctl = f.read()
    special_types = {}
    for st in spectrum_types:
        m = re.search(f"enum class {st} {{(.*?)}}", ctl, flags=re.DOTALL)
        special_types[st] = [x[:-1] for x in m.group(1).split()]
    return special_types


def update_special_types_source(special_types, test_only):
    """
    Modify the source to incorporate new special types that may have been added
    since the last config run. This modifies the list in {special_type}.hh.
    """
    for st in spectrum_types:
        srcname = os_path_join(spectrum_path, 'src', f'{st.lower()}.hh')
        with open(srcname, 'r') as f:
            content = f.read()
        m = re.search(f"^(.*?)Fields<(.*?)>;(.*?)$", content, flags=re.DOTALL)
        newlist = ""
        for special_type in special_types[st]:
            newlist += f"{st}{special_type}<HConfig>,\n"
        content = m.group(1) + "Fields<\nFConfig<>,\n" + newlist[:-2] + "\n>;" + m.group(3)
        if test_only:
            fname = f"CONFIG.TEST.{st.lower()}.TEST.hh"
        else:
            fname = srcname
        with open(fname, 'w') as f:
            f.write(content)
        print(f"[config.py] file {fname} was updated")


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




from math import modf as math_modf


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



