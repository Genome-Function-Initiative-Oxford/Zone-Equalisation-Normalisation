from . import norm_compare
from . import zone_norm
from . import reverse_norm

from importlib.metadata import version, PackageNotFoundError

# Package modules
__all__ = ["norm_compare", "zone_norm", "reverse_norm"]

try:
    # Set attribute version based on toml file
    __version__ = version("your-package-name")
except PackageNotFoundError:
    __version__ = "0.0.0"