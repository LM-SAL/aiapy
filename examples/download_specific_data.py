"""
============================================
Requesting specific AIA images from the JSOC
============================================

This example shows how to request a specific series of AIA images from the JSOC.

We will be filtering the data we require by keywords and requesting short exposure images from a recent flare.

Unfortunately, this can not be done using the sunpy downloader,
`Fido <https://docs.sunpy.org/en/stable/api/sunpy.net.Fido.html>`__
and instead we will use the `DRMS <https://docs.sunpy.org/projects/drms/en/stable/>`__ Python library.
"""
import os
from pathlib import Path

import drms
import matplotlib.pyplot as plt

import sunpy.map as smap

from aiapy.calibrate import normalize_exposure, register, update_pointing

#####################################################
# Exporting data from the JSOC requires registering your
# email first. Please replace this with your email
# address once you have registered.
# See `this page <http://jsoc.stanford.edu/ajax/register_email.html>`__
# for more details.

jsoc_email = os.environ.get("JSOC_EMAIL", "nabil.freij@gmail.com")

#####################################################
# Our goal is to request data of a recent (of time of writing)
# X-class flare. However, we will request the explanation of
# the keywords we want from the JSOC.

client = drms.Client(email=jsoc_email)
keys = [
    "EXPTIME",
    "QUALITY",
    "T_OBS",
    "T_REC",
    "WAVELNTH"
]

print("Querying series info")
# We plan to only use the EUV 12s data for this example.
series_info = client.info("aia.lev1_euv_12s")
for key in keys:
    linkinfo = series_info.keywords.loc[key].linkinfo
    note_str = series_info.keywords.loc[key].note
    print(f"{key:>10} : {note_str}")

#####################################################
# We will construct the query. The X-class flare occurred
# on the 2021/07/03 at 14:30:00 UTC. We will focus on the 5 minutes
# before and after this time.

qstr = "aia.lev1_euv_12s[2021-07-03T14:25:00Z-2021-07-03T14:35:00Z]"
print(f"Querying data -> {qstr}")
results = client.query(qstr, key=keys)
print(f"{len(results)} records retrieved.")

#####################################################
# As you can see from the output, we have received a
# a list of AIA images that were taken during the flare.
# What we want to do now is to filter the list of images
# to only include shorter expsoures.
# However, before we do this, let us check what the exposure times are.

# Filter out entries with EXPTIME > 2 seconds
results = results[results.EXPTIME < 2]
print(results)

#####################################################
# This style of filtering can be done to any column
# in the results. For example, we can filter the WAVELNTH
# column to only include 171 data with short expsoures.

# Only use entries with WAVELNTH == 211
results = results[results.WAVELNTH == 211]
print(results)

#####################################################
# .. note::
#   **Only complete searches can be downloaded from JSOC**,
#   this means that no slicing operations performed on the results object
#   will affect the number of files downloaded.
#
# We can filter and do analysis on the metadata that was returned.
# The issue is is that if we only want this data, you can not use
# this "filtered results" to download only the data we want.
# To do this, we will have to do a second query to the JSOC.
# This time using the query string syntax the
# `lookdata <http://jsoc.stanford.edu/ajax/lookdata.html>`__ web page.
# You can use the website to validate the string before you export the query.

updated_qstr = "aia.lev1_euv_12s[2021-07-03T14:25:00Z-2021-07-03T14:35:00Z][? EXPTIME<2.0 AND WAVELNTH=211 ?]"
print(f"Querying data -> {updated_qstr}")
records = client.query(updated_qstr, key=keys)
print(f"{len(records)} records retrieved. \n")

# We do a quick comparision to ensure the final results are the same.
# For this to work, we just need to deal with the different indexes.
print("Quick Comparison")
print(results.reset_index(drop=True) == records.reset_index(drop=True))

#####################################################
# From here you can now request (export) the data.
# This will download this specific subset of data to your
# local machine when the export request has been completed.
# Depending on the status of the JSOC, this might take a while.
#
# Please be aware the script will hold until the export is complete.

# TODO: Avoid downloading the spike data
export = client.export(updated_qstr, method="url", protocol="fits")
files = export.download(Path("~/sunpy/").expanduser().as_posix())

#####################################################
# Now we will "prep" the data to level 1.5 and plot the
# data sequence using sunpy.

img_files = sorted([f for f in files if f.endswith("image_lev1.fits")])
maps = []
for img_file in img_files:
    basic_map = smap.Map(img_file)
    map_updated_pointing = update_pointing(basic_map)
    map_registered = register(map_updated_pointing)
    map_normalized = normalize_exposure(map_registered)
    maps.append(map_normalized)
sequence = smap.Map(maps, sequence=True)
sequence.peek()

plt.show()
