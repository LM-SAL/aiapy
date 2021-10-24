import os
from contextlib import nullcontext

import numpy as np
import pytest

import astropy.units as u

from aiapy.calibrate import estimate_error
from aiapy.calibrate.util import get_error_table
from aiapy.tests.data import get_test_filepath

# These are not fixtures so that they can be easily used in the parametrize mark
CHANNELS = [94, 131, 171, 193, 211, 304, 335, 1600, 1700, 4500] * u.angstrom
table_local = get_error_table(get_test_filepath('aia_V3_error_table.txt'))


@pytest.mark.parametrize('channel', CHANNELS)
def test_error_all_channels(channel):
    intensity = 10.0 * u.ct / u.pix
    error = estimate_error(intensity, channel, error_table=table_local)
    assert error.unit == intensity.unit


@pytest.mark.parametrize('counts', [
    1,
    np.random.rand(1),
    np.random.rand(10),
    np.random.rand(100, 200),
    np.random.rand(10, 10, 5),
])
def test_counts_shapes(counts):
    counts = counts * 1000 * u.ct / u.pix
    errors = estimate_error(counts, 171*u.angstrom, error_table=table_local)
    if counts.shape == ():
        assert errors.shape == (1,)
    else:
        assert counts.shape == errors.shape


@pytest.mark.parametrize('include_preflight,include_eve,include_chianti,expectation', [
    (False, True, False, nullcontext()),
    (True, False, False, nullcontext()),
    (False, False, False, nullcontext()),
    (False, False, True, nullcontext()),
    (True, True, False, pytest.raises(ValueError, match='Cannot include both EVE and pre-flight correction.'))
])
def test_flags(include_preflight, include_eve, include_chianti, expectation):
    with expectation:
        errors = estimate_error(1*u.ct / u.pix,
                                94 * u.angstrom,
                                error_table=table_local,
                                include_eve=include_eve,
                                include_preflight=include_preflight,
                                include_chianti=include_chianti)
        assert isinstance(errors, u.Quantity)


@pytest.mark.parametrize(
    'channel,counts,include_eve,include_preflight,include_chianti',
    [[c, 10*u.ct / u.pixel] + 3*[False] for c in CHANNELS] +
    [
        [171*u.angstrom, 1000*u.ct/u.pix, True, False, False],
        [171*u.angstrom, 1000*u.ct/u.pix, False, False, True],
    ]
)
def test_error_consistent(idl_environment, channel, counts, include_eve, include_preflight, include_chianti):
    idl = """
    common aia_bp_error_common,common_errtable
    common_errtable=aia_bp_read_error_table('{{ error_table }}')
    data = {{ data }}
    channel = {{ channel }}
    error=aia_bp_estimate_error(data,channel,n_sample=1{{ include_eve }}{{ include_preflight }}{{ include_chianti }})
    """
    error_table = os.path.join(
        idl_environment.ssw_home,
        'sdo',
        'aia',
        'response',
        'aia_V3_error_table.txt'
    )
    ssw = idl_environment.run(
        idl,
        save_vars=['error'],
        args={
            'channel': channel.to('angstrom').value,
            'data': counts.to('ct pixel-1').value,
            'error_table': error_table,
            'include_eve': ',/evenorm' if include_eve else '',
            # NOTE: use of this keyword is actually broken in SSW so these
            # tests only set it to False until it works, but consistency with
            # these results has been verified
            'include_preflight': ',/cal' if include_preflight else '',
            'include_chianti': ',/temperature' if include_chianti else '',
        },
        verbose=False
    )
    error_ssw = ssw['error'] * counts.unit
    error = estimate_error(counts,
                           channel,
                           include_eve=include_eve,
                           include_preflight=include_preflight,
                           include_chianti=include_chianti,
                           error_table=error_table,
                           compare_idl=True)
    assert u.allclose(error, error_ssw, rtol=1e-4)
