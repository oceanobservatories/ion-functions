import warnings
import functools

from importlib.metadata import PackageNotFoundError, version as _package_version

try:
    __version__ = _package_version("ion-functions")  # set from version.py at install time
except PackageNotFoundError:  # running from a source tree that was never installed
    __version__ = "unknown"


def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emmitted
    when the function is used."""

    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.simplefilter('always', DeprecationWarning)  # turn off filter
        warnings.warn("Call to deprecated function {}.".format(func.__name__),
                      category=DeprecationWarning, stacklevel=2)
        warnings.simplefilter('default', DeprecationWarning)  # reset filter
        return func(*args, **kwargs)

    return new_func
