import astropy.units as u
from astropy.table import Table
from astropy.tests.helper import assert_quantity_allclose

from aiapy.calibrate import correct_aec


def test_correct_aec() -> None:
    table = Table(
        {
            "T_OBS": ["2014-09-09T01:00:00Z", "2014-09-09T01:00:10Z", "2014-09-09T01:00:20Z", "2014-09-09T01:00:30Z"],
            "AECTYPE": [0, 2, 0, 2],
            "DATAMEAN": [10.0, 1.0, 10.0, 1.0],
            "EXPTIME": [2.0, 0.1, 2.0, 0.1],
            "NSATPIX": [100, 50, 300, 50],
        }
    )
    assert_quantity_allclose(correct_aec(table), [5, 10, 5, 10] * u.DN / u.s)
    # The exposure time bias only applies to AEC-on frames.
    rate = correct_aec(table, dn_bias=-0.1 * u.DN, exptime_bias=-0.01 * u.s)
    assert_quantity_allclose(rate, [5.05, 10, 5.05, 10] * u.DN / u.s)
    # Each frame loses the NSATPIX of the other AEC state, interpolated to its time.
    rate = correct_aec(table, saturation_scale=100, saturation_shift=0 * u.s)
    assert_quantity_allclose(rate, [4.5, 8, 4.5, 7] * u.DN / u.s)
