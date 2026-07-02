import pytest

import astropy.units as u

from aiapy.utils.decorators import _all_channels, validate_channel


@validate_channel("channel")
def identity(channel):
    return channel


@pytest.mark.parametrize("channel", _all_channels)
def test_valid_channel(channel) -> None:
    assert identity(channel) == channel


@pytest.mark.parametrize("channel", [c.to(u.nm) for c in _all_channels])
def test_valid_channel_equivalent_unit(channel) -> None:
    assert identity(channel) == channel


@pytest.mark.parametrize("channel", [1 * u.angstrom, 94, "94", None])
def test_invalid_channel_raises_error(channel) -> None:
    with pytest.raises(ValueError, match="not in list of valid channels"):
        identity(channel)


def test_custom_valid_channels() -> None:
    @validate_channel("channel", valid_channels=[94 * u.angstrom])
    def f(channel):
        return channel

    assert f(94 * u.angstrom) == 94 * u.angstrom
    with pytest.raises(ValueError, match="not in list of valid channels"):
        f(171 * u.angstrom)


def test_missing_argument_raises_error() -> None:
    with pytest.raises(ValueError, match="Did not find channel in function signature"):

        @validate_channel("channel")
        def f(not_channel):
            return not_channel
