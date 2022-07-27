from pathlib import Path
from urllib.parse import urljoin

from sunpy import log
from sunpy.util.config import get_and_create_sample_dir
from sunpy.util.parfive_helpers import Downloader

_BASE_URLS = (
    "https://github.com/sunpy/data/raw/main/aiapy/",
    "http://data.sunpy.org/aiapy/",
)
_SAMPLE_DATA = {
    "AIA_094_IMAGE": "aia_lev1_94a_2019_01_01t00_00_11_12z_image_lev1.fits",
    "AIA_193_IMAGE": "aia_lev1_193a_2013_03_15t12_01_06_84z_image_lev1.fits",
    "AIA_171_IMAGE": "aia_lev1_171a_2019_01_01t00_00_09_35z_image_lev1.fits",
}
_SAMPLE_FILES = {v: k for k, v in _SAMPLE_DATA.items()}


def _download_sample_data(base_url, sample_files, overwrite):
    """
    Downloads a list of files.

    Parameters
    ----------
    base_url : str
        Base URL for each file.
    sample_files : list of tuples
        List of tuples that are (URL_NAME, SAVE_NAME).
    overwrite : bool
        Will overwrite a file on disk if True.

    Returns
    -------
    `parfive.Results`
        Download results. Will behave like a list of files.
    """
    dl = Downloader(overwrite=overwrite, progress=True, headers={"Accept-Encoding": "identity"})
    for url_file_name, fname in sample_files:
        url = urljoin(base_url, url_file_name)
        dl.enqueue_file(url, filename=fname)
    results = dl.download()
    return results


def _retry_sample_data(results, new_url_base):
    dl = Downloader(overwrite=True, progress=True, headers={"Accept-Encoding": "identity"})
    for err in results.errors:
        file_name = err.url.split("/")[-1]
        log.debug(f"Failed to download {_SAMPLE_FILES[file_name]} from {err.url}: {err.exception}")
        new_url = urljoin(new_url_base, file_name)
        log.debug(f"Attempting redownload of {_SAMPLE_FILES[file_name]} using {new_url}")
        dl.enqueue_file(new_url, filename=err.filepath_partial)
    extra_results = dl.download()
    new_results = results + extra_results
    new_results._errors = extra_results._errors
    return new_results


def _handle_final_errors(results):
    for err in results.errors:
        file_name = err.url.split("/")[-1]
        log.debug(f"Failed to download {_SAMPLE_FILES[file_name]} from {err.url}: {err.exception}")
        log.error(f"Failed to download {_SAMPLE_FILES[file_name]} from all mirrors," "the file will not be available.")


def download_sample_data(overwrite=False):
    """
    Download all sample data at once. This will overwrite any existing files.

    Parameters
    ----------
    overwrite : `bool`
        Overwrite existing sample data.
    """
    sampledata_dir = Path(get_and_create_sample_dir()).parent / Path("aiapy")
    already_downloaded = []
    to_download = []
    for url_file_name in _SAMPLE_FILES.keys():
        fname = sampledata_dir / url_file_name
        # We want to avoid calling download if we already have all the files.
        if fname.exists() and not overwrite:
            already_downloaded.append(fname)
        else:
            to_download.append((url_file_name, fname))
    if to_download:
        results = _download_sample_data(_BASE_URLS[0], to_download, overwrite=overwrite)
    else:
        return already_downloaded
    if results.errors:
        for next_url in _BASE_URLS[1:]:
            results = _retry_sample_data(results, next_url)
            if not results.errors:
                break
        else:
            _handle_final_errors(results)
    return results + already_downloaded
