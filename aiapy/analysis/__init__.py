"""
Subpackage for AIA analysis tools.
"""


def isolate_fe18_intensity(m_94, m_171, m_193):
    """
    Isolate the hot Fe XVIII emission from the 94 Å channel observations.

    Parameters
    ----------
    """
    # NOTE: Should we check that they have the same WCS
    # NOTE: Document these numbers and where they come from
    c_fit = [-7.19e-2, 9.75e-1, 9.79e-2, -2.81e-3]
    A_fit = 0.39
    B_fit = 116.32
    f_fit = 0.31
    I_max_fit = 27.5
    I_94_warm = 0.0
    x = (f_fit * m_171.data + (1 - f_fit) * m_193.data) / B_fit
    x[x >= I_max_fit] = I_max_fit
    for i, cf in enumerate(c_fit):
        I_94_warm += cf * x**i
    # Make metadata modifications to indicate this has been filtered?
    return m_94 - A_fit * I_94_warm * m_193.unit
