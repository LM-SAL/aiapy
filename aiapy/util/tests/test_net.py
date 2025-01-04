import pytest

from aiapy.util.net import get_data_from_jsoc


@pytest.mark.remote_data
def test_get_data_from_jsoc():
    assert get_data_from_jsoc("aia.master_pointing3h[2010-05-13T00:00:00Z]", key="**ALL**") is not None


def test_get_data_from_jsoc_error():
    with pytest.raises(OSError, match="Unable to query the JSOC"):
        get_data_from_jsoc("abc", key="def")
