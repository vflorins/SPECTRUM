

# Config-SPECTRUM

![Development Status](https://img.shields.io/badge/status-experimental-orange)

> ⚠️ **DISCLAIMER:** This project is currently under active pre-alpha development.

Config-SPECTRUM is a configuration assistant for [SPECTRUM](https://github.com/vflorins/SPECTRUM), facilitating C++-native compile-time configuration and parameter management. 

## Quick Start 🥚

### Installation

After obtaining SPECTRUM and adding the variable `SPECTRUM` to the user environment, 
defining the source code location, 
install ``config-spectrum``:
```bash
pip install config-spectrum
```
This defines a command line entry point ``config-spectrum``. 
The only dependency of ``config-spectrum`` is Tomli, 
so the effect on system Python should be negligible, but to be safe, 
you may do this in a Python virtual environment if you prefer 
(this step is outside the scope of this guide). 
To remove ``config-spectrum`` from your system:
```bash
pip uninstall config-spectrum
```
in the environment where it was installed.

### Usage

Compile-time constants (parameters) may be set at the command line. 
The tool will also resolve a set of contingent defaults.
For example:
```bash
config-spectrum . -b Dipole -t Lorentz -d None
```
This will generate a config file with a default name (``config.hh``) in the current directory (``.``). 
This file will inject the default values for a Lorentz particle trajectory in a Dipole background, 
with no Diffusion model applied. The path argument is usually (``.``), but it can be modified. 
To change the name of the file, use the ``-o`` flag. As in:
```bash
config-spectrum ./my_path -b Dipole -t Lorentz -d None -o my_first_test
```
Now the generated file is in a local subdirectory, and is named ``my_first_test.config.hh``. 
Config files are C++ header files. By convention, suffixed with ``config.hh`` instead of the 
usual ``.hh`` for a C++ header file. 

For more information, use the help flag ``config-spectrum -h`` or ``config-spectrum --help``.


## Guide for Developers 💾

Except for the frontend in ``main.py``, the implementation is in the 
submodule ``_impl``, and should not need to be frequently modified, 
in the event that the SPECTRUM source has been updated.

The files ``config_parameters.py`` and ``config_physical.py`` in ``config-spectrum``, 
and the file ``common/compiletime_lists.hh`` are 
the "source of truth" for SPECTRUM configuration. 

* ``compiletime_lists.hh``. This file and included header files determine 
fundamental structure is exposed at compile-time. In particular, it determines what Backgrounds, 
Trajectories, and Diffusion classes are recognized by ``config-spectrum``. 
The enum classes there may be edited when the source code is changed, 
to reflect changes to the SPECTRUM object model. For example, 
a new Background implemented in ``background_silly.hh`` / ``background_silly.cc``
can be registered in ``compiletime_lists.hh`` by adding ``Silly`` to the 
list of backgrounds.

* ``config_parameters.py``. Contains a specification of all parameters
for the SPECTRUM test particle trajectory solver. 
It provides type information, possible ranges,
a docstring, and a globally-defined fallthrough default value.

* ``config_physical.py``. Contains the default values for individual 
cases, so that these can be tuned individually. 
See below for information about setting/modifying these. 

This leads to the following general heuristic for developers:
* **Modifications of ``compiletime_lists.hh`` should be made in tandem with review and
  possible modifications of ``config_parameters.py``, and vice versa.**
* **Modifications of ``config_physical.py`` can be freely made.**

This distinction explains why ``config_physical.py`` is isolated 
in its own file. 

Whenever changes to either of the ``.py`` files are made, 
the ``.config.hh`` files that appear in ``/src/`` need to be updated. 
This process has been automated. You only need to run:
```bash
config-spectrum . --mkdefaults
```
To proceed cautiously, it is recommended to run this first in testing mode:
```bash
config-spectrum . --mkdefaults --test
```

## Setting default values in ``config_physical.py``

### The ``physical_defaults`` dict

The dict of dicts in ``config_physical.py``, called ``physical_defaults``, 
performs two functions: 

* It defines what parameters are used by each trajectory/background/diffusion type. 
For example, if ``r_ref`` does not appear in the dict for the ``Waves`` background, 
then there will be no ``r_ref`` parameter in the ``Waves`` background ``Configure`` class.
* It defines a class-specific "reasonable" default value for each particular trajectory/background/diffusion 
type. Its purpose is to allow default values to depend
on the choice of Background, Trajectory, and Diffusion.

This is a direct evolution of the approach used in previous versions of 
SPECTRUM, where these values were defined in the ``.hh`` file. 
The value may be adjusted:

* after the ``config.hh`` file is generated, by manual editing of that file.  
The ``config.hh`` file is "owned" and managed by the user, so it may be 
changed at any time. 
* during the generation of ``config.hh`` at the command line, by passing 
an optional argument (run ``config-spectrum -h`` to see options). 
* by modifying ``config_physical.py`` (this edit makes the value 
serve as the new default), and re-running ``config-spectrum``. 

This file is the unique location in the code where
these default values are defined. All default values are class-specific. 
There is no inheritance mechanism in the lists of defaults, so
each class repeats parameters defined in base classes.
In some cases, this is indicated by commented-out dividers.

The dict of dicts ``physical_defaults`` has the following structure:
```
physical_defaults = {
    'Background': ..., # dict of backgrounds
    'Trajectory': ..., # dict of trajectories
    'Diffusion': ..., # dict of diffusions
}
```
This defines a simple but efficient "tree" of default values.

### Formatting

The formatting of values of key-value pairs 
is performed by the ``value`` method in ``_impl/configdata.py``. 
Each parameter (say, a floating point parameter) has its type resolved 
via the dictionaries in ``config.parameters.py``. If this type is ``float``, 
then the value is set as a float. If this type is ``string`` (for example,
a file name), then this value is set as a string. Default values of strings can 
be set using an ordinary double-quote "..." string.
If the type is type (a C++ type), then the type is set here
inside of a double-quote string. Brace-initializer expressions can be 
used for default values.

### Last Updated

Lucius Schoenbaum June 3, 2026
