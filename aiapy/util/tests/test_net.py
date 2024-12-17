import pytest

from aiapy.util.net import get_data_from_jsoc


@pytest.mark.remote_data
@pytest.mark.xfail(reason="JSOC is currently down")
def test_get_data_from_jsoc():
    assert get_data_from_jsoc("aia.lev1[2001-01-01]", key="T_OBS") is not None


def test_get_data_from_jsoc_error():
    with pytest.raises(OSError, match="Unable to query the JSOC."):
        get_data_from_jsoc("aia.lev1[2001-01-01]", key="T_OBS")
