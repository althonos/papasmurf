__version__ = "0.1.0"

from . import lib
from .lib import (
    Builder,
    Database,
    Mapper,
    MapperResult,
)

__all__ = ["Builder", "Database", "Mapper", "MapperResult"]
__author__ = lib.__author__
__license__ = "GPLv3"
__doc__ = lib.__doc__
