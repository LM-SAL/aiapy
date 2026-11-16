"""
Contains all the the IDL specific tests for aiapy.
"""

from pathlib import Path

import numpy as np
import pytest

import astropy.time
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.tests.helper import assert_quantity_allclose

import aiapy.psf
from aiapy.calibrate import correct_aec, estimate_error, fetch_spikes, register, respike, update_pointing
from aiapy.calibrate.utils import get_correction_table, get_error_table
from aiapy.response import Channel, get_source_spectra, temperature_response
from aiapy.tests.data import get_test_filepath
from aiapy.tests.helpers import CHANNELS

FUV_CHANNELS = [1600, 1700, 4500] * u.angstrom
# An epoch inside the correction table and the latest epoch
OBSTIMES = astropy.time.Time(["2015-01-01T00:00:00", "2025-01-01T00:00:00"], scale="utc")


@pytest.mark.parametrize(
    ("channel", "n_sample", "include_eve", "include_chianti"),
    [(c, n, False, False) for c in [*CHANNELS, *FUV_CHANNELS] for n in (1, 4)]
    + [
        (171 * u.angstrom, 1, True, False),
        (171 * u.angstrom, 1, False, True),
    ],
)
def test_error_consistent(idl_environment, channel, n_sample, include_eve, include_chianti) -> None:
    # Spans the 25 DN switch in the compression noise
    counts = [0.5, 10, 24.9, 25, 1000, 1e4] * u.DN / u.pix
    idl = """
    common aia_bp_error_common,common_errtable
    common_errtable=aia_bp_read_error_table('{{ error_table }}')
    data = [{{ data }}]
    channel = replicate({{ channel }}, n_elements(data))
    error=aia_bp_estimate_error(data,channel,n_sample={{ n_sample }}{{ include_eve }}{{ include_chianti }})
    """
    error_table_path = Path(idl_environment.ssw_home) / "sdo" / "aia" / "response" / "aia_V3_error_table.txt"
    error_table = get_error_table(error_table_path)
    ssw = idl_environment.run(
        idl,
        save_vars=["error"],
        args={
            "channel": channel.to("angstrom").value,
            "data": ",".join(map(str, counts.to_value("DN pix-1"))),
            "n_sample": n_sample,
            "error_table": error_table_path.as_posix(),
            "include_eve": ",/evenorm" if include_eve else "",
            # NOTE: The /cal keyword (include_preflight) is broken in SSW so it is not tested here
            "include_chianti": ",/temperature" if include_chianti else "",
        },
        verbose=False,
    )
    error_ssw = ssw["error"] * counts.unit
    error = estimate_error(
        counts,
        channel,
        n_sample=n_sample,
        include_eve=include_eve,
        include_chianti=include_chianti,
        error_table=error_table,
        compare_idl=True,
    )
    assert_quantity_allclose(error, error_ssw, rtol=1e-4)


@pytest.mark.remote_data
@pytest.mark.parametrize("channel", [*CHANNELS, *FUV_CHANNELS])
@pytest.mark.parametrize(
    ("kwargs", "keywords"),
    [
        ({}, ",evenorm=0"),
        ({"include_crosstalk": False}, ",/noblend,evenorm=0"),
        *[
            (
                {"obstime": t, "include_eve_correction": eve},
                f",evenorm={int(eve)},timedepend_date='{t.tai.isot}',version=8,respversion='8'",
            )
            for t in OBSTIMES
            for eve in (False, True)
        ],
        ({"include_eve_correction": True}, ",/evenorm,version=8,respversion='8'"),
    ],
)
def test_wavelength_response_consistent(idl_environment, channel, kwargs, keywords) -> None:
    if channel == 4500 * u.angstrom and ("obstime" in kwargs or "include_eve_correction" in kwargs):
        pytest.skip("The version 8 correction table has no 4500 Å calibration")
    fuv = channel in FUV_CHANNELS
    kwargs = dict(kwargs)
    obstime = kwargs.pop("obstime", None)
    c = Channel(
        channel,
        correction_table=get_correction_table(get_test_filepath("aia_V8_20171210_050627_response_table.txt")),
        calibration_version=8,
        **kwargs,
    )
    r = c.wavelength_response(obstime=obstime) / c.pixel_solid_angle
    ssw = idl_environment.run(
        "r = aia_get_response(/area,/dn{{ uv }}{{ keywords }})",
        args={"uv": ",/uv" if fuv else "", "keywords": keywords},
        save_vars=["r"],
        verbose=False,
    )
    r_ssw = ssw["r"][f"A{channel.to_value('angstrom'):.0f}"][0]["ea"][0] * u.cm**2 * u.DN / u.ph
    assert_quantity_allclose(r, r_ssw, rtol=1e-4, atol=0.0 * u.cm**2 * u.DN / u.ph)


@pytest.mark.remote_data
@pytest.mark.parametrize("channel", CHANNELS)
@pytest.mark.parametrize("obstime", [None, OBSTIMES[0]])
@pytest.mark.parametrize(("include_eve_correction", "chiantifix"), [(False, False), (True, False), (True, True)])
def test_temperature_response_consistent(idl_environment, channel, obstime, include_eve_correction, chiantifix) -> None:
    c = Channel(
        channel,
        include_eve_correction=include_eve_correction,
        correction_table=get_correction_table(get_test_filepath("aia_V10_20201119_190000_response_table.txt")),
        calibration_version=10,
    )
    _, response = temperature_response(c, get_source_spectra(), obstime, chiantifix=chiantifix)
    # Always set evenorm: with timedepend_date and no evenorm, IDL waits for input
    keywords = f",evenorm={int(include_eve_correction)}"
    if chiantifix:
        keywords += ",/chiantifix"
    if obstime is not None:
        keywords += f",timedepend_date='{obstime.tai.isot}'"
    ssw = idl_environment.run(
        "r = aia_get_response(/temperature,/dn,version=10,respversion='20201119_190000'{{ keywords }})",
        args={"keywords": keywords},
        save_vars=["r"],
        verbose=False,
    )
    r_ssw = ssw["r"][f"A{channel.to_value('angstrom'):.0f}"][0]["tresp"][0] * u.Unit("cm5 DN / (pix s)")
    # TODO: rtol=1e-5 once sunkit-instruments folds on the spectra wavelengths (sunpy/sunkit-instruments#231)
    assert_quantity_allclose(response, r_ssw, rtol=2e-4)


@pytest.mark.parametrize("channel", CHANNELS)
@pytest.mark.parametrize("use_preflightcore", [False, True])
def test_filter_mesh_parameters_consistent(idl_environment, channel, use_preflightcore) -> None:
    ssw = idl_environment.run(
        """
        resolve_routine, 'aia_calc_psf', /is_function
        m = call_function('aia_getmeshinfo', {{ channel }}, use_preflightcore={{ preflight }})
        """,
        args={"channel": f"{channel.to_value('angstrom'):.0f}", "preflight": int(use_preflightcore)},
        save_vars=["m"],
        verbose=False,
    )
    m = ssw["m"][0]
    params = aiapy.psf.filter_mesh_parameters(use_preflightcore=use_preflightcore)[channel]
    for key, expected in [
        ("angle_arm", [m[f"angle{i}"] for i in range(1, 5)] * u.deg),
        ("error_angle_arm", [m[f"delta{i}"] for i in range(1, 5)] * u.deg),
        ("spacing_e", m["spacing"] * u.pix),
        ("spacing_fp", m["fp_spacing"] * u.pix),
        ("mesh_pitch", m["meshpitch"] * u.um),
        ("mesh_width", m["meshwidth"] * u.um),
        ("width", m["gs_width"] * u.pix),
    ]:
        assert_quantity_allclose(params[key], expected, rtol=1e-6)


@pytest.mark.remote_data
def test_register_consistent(idl_environment, aia_193_level1_file, aia_193_level1_map) -> None:
    ssw = idl_environment.run(
        """
        read_sdo, '{{ file }}', index, data, /uncomp_delete
        aia_prep, index, data, oindex, odata, /no_mpo_update
        header = [oindex.crpix1, oindex.crpix2, oindex.cdelt1, oindex.cdelt2, oindex.crota2]
        """,
        args={"file": aia_193_level1_file},
        save_vars=["header", "odata"],
        verbose=False,
    )
    registered = register(aia_193_level1_map)
    meta = registered.meta
    np.testing.assert_allclose([meta["crpix1"], meta["crpix2"], meta["cdelt1"], meta["cdelt2"]], ssw["header"][:4])
    np.testing.assert_allclose(registered.rotation_matrix, np.identity(2), atol=1e-8)
    assert ssw["header"][4] == 0
    # SSW interpolates with cubic convolution and aiapy with a cubic spline,
    # so only compare the bulk intensities
    bright = ssw["odata"] > 50
    np.testing.assert_allclose(np.median(registered.data[bright] / ssw["odata"][bright]), 1, atol=5e-3)
    np.testing.assert_allclose(registered.data.sum() / ssw["odata"].sum(), 1, atol=5e-3)


@pytest.mark.remote_data
def test_update_pointing_consistent(idl_environment, aia_193_level1_file, aia_193_level1_map, pointing_table) -> None:
    # aia2wcsmin queries the JSOC for the master pointing table
    ssw = idl_environment.run(
        """
        read_sdo, '{{ file }}', index
        nindex = aia2wcsmin(index)
        header = [nindex.crpix1, nindex.crpix2, nindex.cdelt1, nindex.cdelt2, nindex.crota2]
        """,
        args={"file": aia_193_level1_file},
        save_vars=["header"],
        verbose=False,
    )
    meta = update_pointing(aia_193_level1_map, pointing_table=pointing_table).meta
    np.testing.assert_allclose(
        [meta["crpix1"], meta["crpix2"], meta["cdelt1"], meta["cdelt2"], meta["crota2"]],
        ssw["header"],
        rtol=1e-6,
    )


@pytest.mark.remote_data
def test_respike_consistent(idl_environment, aia_193_level1_map, tmp_path) -> None:
    # SSW's own spike download needs a registered JSOC email, so pass it the full-frame spikes
    coords, values = fetch_spikes(aia_193_level1_map)
    index = coords.y.to_value("pix").round() * aia_193_level1_map.data.shape[1] + coords.x.to_value("pix").round()
    fits.writeto(tmp_path / "spikes.fits", np.array([index, values]).astype(np.int32))
    cutout = aia_193_level1_map.submap([500, 1000] * u.pix, top_right=[2499, 2999] * u.pix)
    # BLANK only applies to integer data
    cutout.meta.pop("blank")
    cutout.save(tmp_path / "cutout.fits")
    ssw = idl_environment.run(
        """
        read_sdo, '{{ path }}/cutout.fits', index, data
        spikes = readfits('{{ path }}/spikes.fits')
        aia_respike, index, data, oindex, odata, ispikedd=spikes
        """,
        args={"path": tmp_path.as_posix()},
        save_vars=["odata"],
        verbose=False,
    )
    respiked = respike(cutout)
    assert (respiked.data != cutout.data).any()
    np.testing.assert_array_equal(respiked.data, ssw["odata"])


# For now, we only need the 94 Angstrom PSF
# In future, we will want to test all channels
@pytest.fixture(params=[CHANNELS[0]])
def psf_idl(idl_environment, request):
    """
    The point spread function as calculated by aia_calc_psf.pro.
    """
    r = idl_environment.run(
        "psf = aia_calc_psf({{channel}},/use_preflightcore)",
        args={"channel": f"{request.param.value:.0f}"},
        save_vars=["psf"],
        verbose=False,
    )
    return r["psf"]


def test_psf_consistent(psf_94, psf_idl) -> None:
    """
    Check whether PSF is consistent with IDL calculation.

    **WARNING** This test will take a very long time to run.
    """
    # NOTE: Only compare values above some threshold as the
    # rest of the PSF is essentially noise
    i_valid = np.where(psf_idl > 1e-10)
    # NOTE: The IDL and Python PSF functions have been found to
    # agree within 0.2% for all points along the PSF arms for
    # both the preflight and non-preflight cases.
    np.testing.assert_allclose(psf_94[i_valid], psf_idl[i_valid], rtol=2e-3)


def test_deconvolve_consistent(idl_environment, aia_171_positive_map, psf_128, tmp_path) -> None:
    fits.writeto(tmp_path / "image.fits", aia_171_positive_map.data)
    fits.writeto(tmp_path / "psf.fits", psf_128)
    ssw = idl_environment.run(
        """
        image = readfits('{{ path }}/image.fits')
        psf = readfits('{{ path }}/psf.fits')
        out = aia_deconvolve_richardsonlucy(image, psf, niter=25)
        """,
        args={"path": tmp_path.as_posix()},
        save_vars=["out"],
        verbose=False,
    )
    map_decon = aiapy.psf.deconvolve(aia_171_positive_map, psf=psf_128, iterations=25)
    np.testing.assert_allclose(map_decon.data, ssw["out"], rtol=1e-10)


def test_aec_consistent(idl_environment, tmp_path) -> None:
    # Interleaved AEC-off/on frames; NSATPIX is nonzero at both ends, where IDL INTERPOL extrapolates.
    n = 40
    flare = np.sin(np.linspace(0, np.pi, n)) ** 2
    aec_on = np.arange(n) % 2 == 1
    table = Table(
        {
            "T_OBS": (astropy.time.Time("2014-09-09T00:00:00") + np.arange(n) * 12 * u.s).isot,
            "AECTYPE": np.where(aec_on, 2, 0),
            "DATAMEAN": 10 + 5 * flare,
            "EXPTIME": np.where(aec_on, 2.9 - 2.7 * flare, 2.9),
            "NSATPIX": np.round((0.1 + flare) * np.where(aec_on, 600, 6000)).astype(int),
        }
    )
    table.write(tmp_path / "hdrs.fits")
    ssw = idl_environment.run(
        """
        hdrs = mrdfits('{{ path }}', 1, /silent)
        aia_aec_correct, hdrs, dnbias = -0.084d, exptimebias = -0.0032d, /shortonly, satscale = 18000, $
            satshift = -12, fixcounts_on = fixcounts_on, fixcounts_off = fixcounts_off, /noplot
        """,
        args={"path": (tmp_path / "hdrs.fits").as_posix()},
        save_vars=["fixcounts_on", "fixcounts_off"],
        verbose=False,
    )
    rate = correct_aec(
        table,
        dn_bias=-0.084 * u.DN,
        exptime_bias=-3.2 * u.ms,
        saturation_scale=18000,
        saturation_shift=-12 * u.s,
    )
    assert_quantity_allclose(rate[aec_on], ssw["fixcounts_on"] * u.DN / u.s, rtol=1e-6)
    assert_quantity_allclose(rate[~aec_on], ssw["fixcounts_off"] * u.DN / u.s, rtol=1e-6)
