"""
Make the figure for the JOSS paper
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
import astropy.time
from astropy.visualization import time_support, ImageNormalize, AsinhStretch
import sunpy.map
from sunpy.net import Fido, attrs

from aiapy.response import Channel
from aiapy.calibrate import degradation
from aiapy.psf import deconvolve


channels = [94, 131, 171, 193, 211, 304, 335] * u.angstrom

# Query data, deconvolve, and crop
q = Fido.search(
    attrs.Time('2011-06-07T06:52:00', '2011-06-07T06:52:10'),
    attrs.Instrument('AIA'),
    attrs.Wavelength(wavemin=171*u.angstrom, wavemax=171*u.angstrom),
)
f = Fido.fetch(q)
m = sunpy.map.Map(f)
# NOTE: this is slow without a GPU
m_decon = deconvolve(m)
left_corner = 500*u.arcsec, -600*u.arcsec
right_corner = 1000*u.arcsec, -100*u.arcsec
m_sub = m.submap(
    SkyCoord(*left_corner, frame=m.coordinate_frame),
    top_right=SkyCoord(*right_corner, frame=m.coordinate_frame),
)
m_decon_sub = m_sub.submap(
    SkyCoord(*left_corner, frame=m_decon.coordinate_frame),
    top_right=SkyCoord(*right_corner, frame=m_decon.coordinate_frame),
)

# Wavelength response
wave_response = {}
for c in channels:
    chan = Channel(c)
    wave_response[c] = (chan.wavelength, chan.wavelength_response())

# Degradation
time_0 = astropy.time.Time('2010-03-25T00:00:00', scale='utc')
now = astropy.time.Time.now()
time = time_0 + np.arange(0, (now - time_0).to(u.day).value, 7) * u.day
degrad = {c: degradation(c,
                         time,
                         calibration_version=10,
                         correction_table='./aia_V10_20201023_180000_response_table.txt')
          for c in channels}

fig = plt.figure(figsize=(8, 8))
norm = ImageNormalize(vmin=0, vmax=1.5e4, stretch=AsinhStretch(0.01))
# Image before deconvolution
ax = fig.add_subplot(221, projection=m_sub)
m_sub.plot(axes=ax, title=r'171 $\mathrm{\AA}$ Level 1', norm=norm)
ax.grid(alpha=0)
# Image after deconvolution
ax = fig.add_subplot(222, projection=m_decon_sub)
m_decon_sub.plot(axes=ax,
                 title=r'171 $\mathrm{\AA}$ Level 1, Deconvolved',
                 norm=norm)
ax.grid(alpha=0)
lon, lat = ax.coords
lat.set_axislabel(' ')
lat.set_ticklabel_visible(False)
# Degradation Curve
ax = fig.add_subplot(223)
with time_support(format='jyear'):
    for c in channels:
        ax.plot(time, degrad[c],)
    ax.set_xlim(time[[0, -1]])
ax.set_xlabel('Date')
ax.set_ylabel('Degradation')
# Wavelength response
ax = fig.add_subplot(224)
for i, c in enumerate(channels):
    w, r = wave_response[c]
    ax.plot(w, r, label=f'{c.value:.0f} {c.unit.to_string(format="latex")}')
ax.set_yscale('log')
ax.set_xlim(80, 360)
ax.set_ylim(1e-3, 5)
ax.legend(frameon=False)
ax.set_xlabel(f'Wavelength [{w.unit.to_string(format="latex")}]')
ax.set_ylabel(f'Response [{r.unit.to_string(format="latex")}]')
# Save figure
plt.subplots_adjust(wspace=0.22)
plt.savefig('figure-1.pdf')
