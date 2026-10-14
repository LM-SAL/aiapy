import re

import pandas as pd
import pytest

import drms

from aiapy.utils.net import _get_data_from_jsoc


def _fake_client(result):
    class FakeClient:
        def query(self, *args, **kwargs):  # NOQA: ARG002
            if isinstance(result, Exception):
                raise result
            return result

    return FakeClient


def test__get_data_from_jsoc_error(monkeypatch):
    error = drms.DrmsQueryError("Bad query")
    monkeypatch.setattr(drms, "Client", _fake_client(error))
    with pytest.raises(OSError, match=re.escape("Unable to query the JSOC.\n Error message: Bad query")) as excinfo:
        _get_data_from_jsoc("abc", key="def")
    assert excinfo.value.__cause__ is error


def test__get_data_from_jsoc_no_data(monkeypatch):
    monkeypatch.setattr(drms, "Client", _fake_client(pd.DataFrame()))
    with pytest.raises(RuntimeError, match="No data found for this query: abc, key: def, seg: None"):
        _get_data_from_jsoc("abc", key="def")
