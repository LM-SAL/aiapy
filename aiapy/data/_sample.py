from pathlib import Path
from urllib.parse import urljoin

from sunpy.data._sample import _download_sample_data
from sunpy.util.config import get_and_create_sample_dir
from sunpy.util.parfive_helpers import Downloader

_BASE_URLS = (
    'https://github.com/sunpy/sample-data/raw/master/aiapy/',
    'http://data.sunpy.org/aiapy/',
)  # NOQA
# Shortcut requirements:
# start with the instrument name then
# the wavelength or energy if needed then
# an optional description if needed then
# a reference name for the class into which the file will be opened
# (e.g. IMAGE for Maps, TIMESERIES for TimeSeries, SPECTRUM for Spectrum)
# All separated by underscores
# the files should include necessary extensions
_SAMPLE_DATA = {
    "AIA_094_IMAGE": "aia_lev1_94a_2019_01_01t00_00_11_12z_image_lev1.fits",
    "AIA_193_IMAGE": "aia_lev1_193a_2013_03_15t12_01_06_84z_image_lev1.fits",
    "AIA_171_IMAGE": "aia_lev1_171a_2019_01_01t00_00_09_35z_image_lev1.fits",
}
# Reverse the dict because we want to use it backwards, but it is nicer to
# write the other way around
_SAMPLE_FILES = {v: k for k, v in _SAMPLE_DATA.items()}  # NOQA


def _retry_sample_data(results):
    # In case we have a broken file on disk, overwrite it.
    dl = Downloader(overwrite=True, progress=True, headers={'Accept-Encoding': 'identity'})
    for err in results.errors:
        file_name = err.filepath_partial().name
        # Update the url to a mirror and requeue the file.
        new_url = urljoin(_BASE_URLS[1], file_name)
        dl.enqueue_file(new_url, filename=err.filepath_partial)
    extra_results = dl.download()
    for err in extra_results.errors:
        file_name = err.filepath_partial().name
    return results + extra_results


def download_sample_data(overwrite=False):
    """
    Download all sample data at once. This will overwrite any existing files.

    Parameters
    ----------
    overwrite: `bool`
        Overwrite existing sample data.
    """
    sampledata_dir = Path(get_and_create_sample_dir()) / Path("aiapy")
    already_downloaded = []
    to_download = []
    for url_file_name in _SAMPLE_FILES.keys():
        fname = sampledata_dir/url_file_name
        if fname.exists() and not overwrite:
            already_downloaded.append(fname)
        else:
            to_download.append((url_file_name, fname))
    if to_download:
        results = _download_sample_data(_BASE_URLS[0], to_download, overwrite=overwrite)
    else:
        return already_downloaded
    return results + already_downloaded
