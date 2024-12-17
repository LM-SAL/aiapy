"""
This only exists to bypass the gzipped issue with parfive.

See https://github.com/Cadair/parfive/issues/121
"""

from pathlib import Path

from sunpy import config
from sunpy.data.data_manager.cache import Cache
from sunpy.data.data_manager.downloader import DownloaderBase, DownloaderError
from sunpy.data.data_manager.manager import DataManager
from sunpy.data.data_manager.storage import SqliteStorage
from sunpy.util.parfive_helpers import Downloader

__all__ = ["manager"]


class AIAParfiveDownloader(DownloaderBase):
    def download(self, url, path) -> None:
        downloader = Downloader()
        path = Path(path)
        filename = path.name
        directory = path.parent
        downloader.enqueue_file(url, directory, filename, max_splits=1)
        try:
            output = downloader.download()
        except Exception as e:
            raise DownloaderError from e
        if output.errors:
            raise DownloaderError(output.errors[0].exception)


_download_dir = config.get("downloads", "remote_data_manager_dir")
manager = DataManager(
    Cache(
        AIAParfiveDownloader(),
        SqliteStorage(f"{_download_dir}/data_manager.db"),
        _download_dir,
    ),
)
