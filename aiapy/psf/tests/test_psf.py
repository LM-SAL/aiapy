import numpy as np
import pytest

import aiapy.psf
from aiapy.conftest import CHANNELS

MESH_PROPERTIES = [
    "angle_arm",
    "error_angle_arm",
    "spacing_e",
    "spacing_fp",
    "mesh_pitch",
    "mesh_width",
    "width",
]


@pytest.mark.parametrize("use_preflightcore", [True, False])
def test_filter_mesh_parameters(use_preflightcore) -> None:
    params = aiapy.psf.filter_mesh_parameters(use_preflightcore=use_preflightcore)
    assert isinstance(params, dict)
    assert all(c in params for c in CHANNELS)
    assert all(all(p in params[c] for p in MESH_PROPERTIES) for c in CHANNELS)


def test_psf(psf) -> None:
    assert isinstance(psf, np.ndarray)
    assert psf.shape == (4096, 4096)
