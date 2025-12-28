from importlib.metadata import version, PackageNotFoundError
from gggears.curve import *
from gggears.defs import *
from gggears.function_generators import *
from src.gggears.conv_spline import *
from src.gggears.core import *
from src.gggears.conv_build123d import *
from src.gggears.wrapper import *
from gggears.gearmath import *


try:
    __version__ = version("gggears")
except PackageNotFoundError:
    __version__ = "unknown version"
