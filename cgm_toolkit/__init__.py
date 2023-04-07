from pkg_resources import DistributionNotFound
from pkg_resources import get_distribution

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass

from .cgm_profiles import HaloProfile
from .xray_emissivity import XrayEmissivity 
