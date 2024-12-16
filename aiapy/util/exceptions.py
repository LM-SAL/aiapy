"""
Errors/exceptions and warnings to be used throughout aiapy.
"""

from astropy.utils.exceptions import AstropyWarning

__all__ = ["AIApyUserWarning", "AIApyWarning"]


class AIApyWarning(AstropyWarning):
    """
    The base warning class from which all aiapy warnings should inherit.

    This warning should not be issued in normal code. Use
    "AIApyUserWarning" instead or a specific sub-class.
    """


class AIApyUserWarning(UserWarning, AIApyWarning):
    """
    The primary warning class for aiapy.

    Use this if you do not need a specific type of warning.
    """
