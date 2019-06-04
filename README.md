# aiapy

Python tools for AIA data analysis

Functionality needed
aia_prep
   - update master pointing keywords by fetching from JSOC
   - image interpolation onto a desired grid (equivalent to AIA level 1.5)
   - optionally normalize to 1 AU
aia_get_response
   - effective areas as function of time
   - temperature response functions
aia_deconv
   - point-spread function deconvolution and desaturation
DEMs
   - Differential Emission Measure inversions a la Cheung et al. (2015), and 
   - by DNNs (http://helioml.org/04/Differential_Emission_Measurements.html)
aia_intscale
   - Byte-scale AIA images for display with "standard" color tables. 
