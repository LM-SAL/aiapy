from sunpy.image.transform import _rotation_registry, add_rotation_function

__all__ = ["_rotation_cupy", "_rotation_function_names"]


@add_rotation_function(
    "cupy", allowed_orders=range(6), handles_clipping=False, handles_image_nans=False, handles_nan_missing=True
)
def _rotation_cupy(image, matrix, shift, order, missing, clip):
    """
    * Rotates using `cupyx.scipy.ndimage.affine_transform` from `cupy <https://docs.cupy.dev/en/stable/index.html>`__
    * Coverts from a numpy array to a cupy array and then back again.
    * The ``order`` parameter is the order of the spline interpolation, and ranges
      from 0 to 5.
    * The ``mode`` parameter for :func:`~cupyx.scipy.ndimage.affine_transform` is fixed to
      be ``'constant'``
    """
    try:
        import cupy
        import cupyx.scipy.ndimage
    except ImportError:
        raise ImportError(
            "cupy or cupy-cuda* (pre-compiled for each cuda version) is required to use this rotation method."
        )
    rotated_image = cupyx.scipy.ndimage.affine_transform(
        cupy.array(image).T, cupy.array(matrix), offset=cupy.array(shift), order=order, mode="constant", cval=missing
    ).T
    return cupy.asnumpy(rotated_image)


# Generate the string with allowable rotation-function names for use in docstrings
_rotation_function_names = ", ".join([f"``'{name}'``" for name in _rotation_registry])
