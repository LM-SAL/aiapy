from sunpy.image.transform import add_rotation_function

__all__ = []


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
        raise ImportError("The cupy package is required to use this rotation method.")
    rotated_image = cupyx.scipy.ndimage.affine_transform(
        cupy.array(image.T), cupy.array(matrix), offset=cupy.array(shift), order=order, mode="constant", cval=missing
    ).T
    return cupy.asnumpy(rotated_image.T)
