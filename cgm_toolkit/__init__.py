try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    from importlib_metadata import version, PackageNotFoundError

try:
    __version__ = version("cgm_toolkit")
except PackageNotFoundError:
    __version__ = "unknown"

from .cgm_profiles import HaloProfile
from .xray_emissivity import XrayEmissivity