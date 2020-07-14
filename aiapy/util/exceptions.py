"""
Errors/exceptions and warnings to be used throughout aiapy
"""

from astropy.utils.exceptions import AstropyWarning

__all__ = ['AiapyWarning', 'AiapyUserWarning']


class AiapyWarning(AstropyWarning):
    """
    The base warning class from which all aiapy warnings should inherit.
    This warning should not be issued in normal code. Use
    "AiapyUserWarning" instead or a specific sub-class.
    """


class AiapyUserWarning(UserWarning, AiapyWarning):
    """
    The primary warning class for aiapy
    Use this if you do not need a specific type of warning.
    """
