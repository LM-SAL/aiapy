import numpy as np
import pytest
from sunkit_instruments.response import SourceSpectra

import astropy.time
import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

from aiapy.response import Channel, get_source_spectra, temperature_response

pytestmark = pytest.mark.remote_data

OBSTIME = astropy.time.Time("2015-01-01T00:00:00", scale="utc")


@pytest.fixture(scope="module")
def spectra():
    return get_source_spectra()


@pytest.mark.parametrize(
    ("channel", "corrected", "expected"),
    # SSW aia_get_response(/temperature,/dn) at log T = 5.6, 6.0, 6.4 and 7.0 with evenorm=0,
    # or corrected with /evenorm,/chiantifix,timedepend_date='2015-01-01T00:00:00'
    [
        (94, False, [2.38798598e-28, 1.62852642e-27, 3.79002406e-28, 2.03470708e-27]),
        (335, False, [2.36254330e-27, 3.63886561e-27, 4.91369007e-27, 4.63133491e-28]),
        (94, True, [1.86155718e-28, 1.69469000e-27, 2.98143060e-28, 1.34597032e-27]),
        (335, True, [9.82866390e-28, 1.51384261e-27, 2.04419570e-27, 1.92673017e-28]),
    ],
)
def test_temperature_response(spectra, channel, corrected, expected) -> None:
    channel = Channel(channel * u.angstrom, include_eve_correction=corrected)
    obstime = OBSTIME if corrected else None
    temperature, response = temperature_response(channel, spectra, obstime, chiantifix=corrected)
    index = [32, 40, 48, 60]
    assert_quantity_allclose(temperature[index], 10 ** np.array([5.6, 6.0, 6.4, 7.0]) * u.K)
    # TODO: rtol=1e-5 once sunkit-instruments folds on the spectra wavelengths (sunpy/sunkit-instruments#231)
    assert_quantity_allclose(response[index], expected * u.Unit("cm5 DN / (pix s)"), rtol=2e-4)


def test_chiantifix_requires_eve(spectra) -> None:
    with pytest.raises(ValueError, match="include_eve_correction=True"):
        temperature_response(Channel(94 * u.angstrom), spectra, chiantifix=True)


def test_chiantifix_requires_ssw_spectra(spectra) -> None:
    other = SourceSpectra(spectra.temperature, spectra.wavelength, spectra.data)
    with pytest.raises(ValueError, match="get_source_spectra"):
        temperature_response(Channel(94 * u.angstrom, include_eve_correction=True), other, chiantifix=True)


def test_temperature_response_obstime_array(spectra) -> None:
    channel = Channel(94 * u.angstrom)
    obstime = astropy.time.Time(["2015-01-01T00:00:00", "2020-01-01T00:00:00"], scale="utc")
    _, response = temperature_response(channel, spectra, obstime)
    assert response.shape == (2, *spectra.temperature.shape)
    for row, t in zip(response, obstime, strict=True):
        assert_quantity_allclose(row, temperature_response(channel, spectra, t)[1], rtol=1e-12)
