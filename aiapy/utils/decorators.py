"""
This module contains decorators for validating arguments related to AIA
channels.
"""

import inspect
import functools

import astropy.units as u

_all_channels = [
    94 * u.angstrom,
    131 * u.angstrom,
    171 * u.angstrom,
    193 * u.angstrom,
    211 * u.angstrom,
    304 * u.angstrom,
    335 * u.angstrom,
    1600 * u.angstrom,
    1700 * u.angstrom,
    4500 * u.angstrom,
]


def validate_channel(argument, *, valid_channels="all"):
    """
    Parameters
    ----------
    argument : str
        Argument name to validate.
    valid_channels : list or str, optional
        List of valid channels.
        Defaults to "all", which will validate against the list of all AIA channels.
    """
    if valid_channels == "all":
        valid_channels = _all_channels

    def outer(function):
        sig = inspect.signature(function)
        if argument not in sig.parameters:
            msg = f"Did not find {argument} in function signature ({sig})."
            raise ValueError(msg)

        @functools.wraps(function)
        def inner(*args, **kwargs):
            all_args = sig.bind(*args, **kwargs)
            channel = all_args.arguments[argument]
            if channel not in valid_channels:
                msg = f'channel "{channel}" not in list of valid channels: {valid_channels}.'
                raise ValueError(msg)
            return function(*args, **kwargs)

        return inner

    return outer
