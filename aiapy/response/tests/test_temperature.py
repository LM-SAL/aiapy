import numpy as np
import pytest

import astropy.time
import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

from aiapy.response import Channel, get_source_spectra, temperature_response

pytestmark = pytest.mark.remote_data

OBSTIME = astropy.time.Time("2015-01-01T00:00:00", scale="utc")


@pytest.fixture(scope="module")
def spectra():
    return get_source_spectra()


@pytest.fixture(scope="module")
def spectra_chianti11():
    return get_source_spectra("chianti11")


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


@pytest.mark.parametrize(
    ("channel", "index", "expected", "rtol"),
    # sdo.aia.temperature_response's aia_chianti11_temperature_response.ecsv (sun-data/solar-dynamics-observatory#28,
    # the same plasma, fiasco 0.8.2 and CHIANTI 11.0.2), raw, at log T 6.0, 6.4 and 7.0. Two known differences:
    # sdo boosts only one member of the He II doublet (see get_source_spectra); He II reaches 335 Å through the
    # short-wavelength wing of its filter and 94 Å through the 304 Å crosstalk, worth 1.8e-3 of 335 Å at log T 6.0,
    # 3.9e-4 at 6.4 and 7.9e-4 at 7.0 (and 1.1e-4 of 94 Å at 5.6), so 335 Å is compared at 6.4 and 7.0 only. And sdo
    # computed the continua on a grid starting at 25 Å, where CHIANTI's free-bound routine (and fiasco's port of it)
    # leaves out recombination onto the excited levels of the H-like and He-like ions (mostly O VIII, O VII and
    # N VII here), worth 2.4e-3 of 94 Å at log T 6.4 and under 1e-4 elsewhere. The rest is the fold, about 4e-5.
    [
        (94, [40, 60], [1.60677102e-27, 1.91762518e-27], 2e-4),
        (94, [48], [1.34141879e-28], 3e-3),
        (335, [48, 60], [4.95033912e-27, 4.59150350e-28], 1.5e-3),
    ],
)
def test_temperature_response_chianti11(spectra_chianti11, channel, index, expected, rtol) -> None:
    assert spectra_chianti11.meta["chianti_version"] == "11.0.2"
    assert "chiantifix" not in spectra_chianti11.meta
    temperature, response = temperature_response(Channel(channel * u.angstrom), spectra_chianti11)
    assert_quantity_allclose(np.log10(temperature[[40, 48, 60]].to_value(u.K)), [6.0, 6.4, 7.0])
    assert_quantity_allclose(response[index], expected * u.Unit("cm5 DN / (pix s)"), rtol=rtol)


def test_chiantifix_requires_ssw_spectra(spectra_chianti11) -> None:
    channel = Channel(94 * u.angstrom, include_eve_correction=True)
    with pytest.raises(ValueError, match="get_source_spectra"):
        temperature_response(channel, spectra_chianti11, chiantifix=True)


def test_temperature_response_obstime_array(spectra) -> None:
    channel = Channel(94 * u.angstrom)
    obstime = astropy.time.Time(["2015-01-01T00:00:00", "2020-01-01T00:00:00"], scale="utc")
    _, response = temperature_response(channel, spectra, obstime)
    assert response.shape == (2, *spectra.temperature.shape)
    for row, t in zip(response, obstime, strict=True):
        assert_quantity_allclose(row, temperature_response(channel, spectra, t)[1], rtol=1e-12)
