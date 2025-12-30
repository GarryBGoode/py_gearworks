from importlib.metadata import version, PackageNotFoundError
from py_gearworks.curve import *
from py_gearworks.defs import *
from py_gearworks.function_generators import *
from py_gearworks.conv_spline import *
from py_gearworks.core import *
from py_gearworks.conv_build123d import *
from py_gearworks.wrapper import *
from py_gearworks.gearmath import *


try:
    __version__ = version("py_gearworks")
except PackageNotFoundError:
    __version__ = "unknown version"
