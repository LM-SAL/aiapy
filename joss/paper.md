---
title: "aiapy: A Python Package for Analyzing Solar EUV Image Data from AIA"
tags:
- Python
- Astronomy
- Solar Physics
authors:
- name: Will T. Barnes^[Other than the first two authors, the author list in this paper is sorted alphabetically.]
  orcid: 0000-0001-9642-6089
  affiliation: "1, 2"
- name: Mark C. M. Cheung
  orcid: 0000-0003-2110-9753
  affiliation: 2
- name: Monica G. Bobra
  orcid: 0000-0002-5662-9604
  affiliation: 3
- name: Paul F. Boerner
  affiliation: 2
  orcid: 0000-0002-4490-9860
- name: Georgios Chintzoglou
  orcid: 0000-0002-1253-8882
  affiliation: "2, 4"
- name: Drew Leonard
  orcid: 0000-0001-5270-7487
  affiliation: 5
- name: Stuart J. Mumford
  orcid: 0000-0003-4217-4642
  affiliation: 6
- name: Nicholas Padmanabhan
  affiliation: "2, 7"
  orcid: 0000-0001-8067-6788
- name: Albert Y. Shih
  orcid: 0000-0001-6874-2594
  affiliation: 8
- name: Nina Shirman
  affiliation: 2
  orcid: 0000-0001-7136-8628
- name: David Stansby
  affiliation: 9
  orcid: 0000-0002-1365-1908
- name: Paul J. Wright
  orcid: 0000-0001-9021-611X
  affiliation: 3

affiliations:
- name: National Research Council Postdoctoral Research Associate residing at the Naval Research Laboratory, Washington, D.C. 20375, USA
  index: 1
- name: Lockheed Martin Solar and Astrophysics Laboratory, Palo Alto, CA 94304, USA
  index: 2
- name: W. W. Hansen Experimental Physics Laboratory, Stanford University, Stanford, CA 94305, USA
  index: 3
- name: University Corporation for Atmospheric Research, Boulder, CO 80301, USA
  index: 4
- name: Aperio Software Ltd, Leeds LS6 3HN, UK
  index: 5
- name: School of Mathematics and Statistics, The University of Sheffield, Sheffield S3 7RH, UK
  index: 6
- name: Princeton University, Princeton, NJ 08544, USA
  index: 7
- name: NASA Goddard Space Flight Center, Greenbelt, MD 20771, USA
  index: 8
- name: Mullard Space Science Laboratory, University College London, Holmbury St. Mary, Surrey RH5 6NT, UK
  index: 9

date: 7 October 2020
bibliography: paper.bib
---

# Summary

The Atmospheric Imaging Assembly [AIA; @lemen_atmospheric_2012] instrument onboard the NASA *Solar Dynamics Observatory* [SDO; @pesnell_solar_2012] spacecraft has observed the full-disk of the Sun, nearly continuously, for the last ten years.
It is one of three instruments on SDO, along with the Helioseismic and Magnetic Imager [HMI; @schou_design_2012] and the Extreme Ultraviolet Variability Experiment [EVE; @woods_extreme_2012].
With its high spatial ($0.6''$ per pixel) and temporal (up to 12 seconds for most channels) resolution, AIA has greatly enhanced our understanding of our closest star in a number of different areas, including the initiation of flares and coronal mass ejections as well as the quiescent heating of the corona, the outermost layer of the Sun's atmosphere.
AIA is a narrowband imaging instrument comprised of four separate telescopes that collectively observe the full-disk of the Sun at ten different wavelengths: seven extreme ultraviolet (EUV) wavelengths, two far UV wavelengths, and one visible wavelength.
It produces nearly 60,000 images per day with a new 4K-by-4K image produced by each EUV channel every 12 seconds.
The image data are provided to the community by the [Joint Science Operations Center (JSOC)](http://jsoc.stanford.edu/) at Stanford University in the Flexible Image Transport System [FITS; @wells_fits_1981] format.

# Statement of Need

`aiapy` is a Python package for analyzing calibrated (level 1) EUV imaging data from AIA.
It includes capabilities for aligning images between channels, deconvolving images with the instrument point-spread function (PSF), computing channel sensitivity as a function of wavelength, and correcting images for telescope degradation, among others.
Historically, most data analysis in solar physics has been done using the Interactive Data Language (IDL), a proprietary interpreted language commonly used throughout astronomy.
The AIA instrument team has developed and continues to maintain a comprehensive set of IDL tools for querying, calibrating, and analyzing AIA data that has been widely used for nearly a decade.
A full description of this software can be found in the [SDO Data Analysis Guide](https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/index.html).

As the solar physics community, and the astronomy community as a whole, began transitioning to Python, some of this functionality was absorbed into `sunpy`, the core Python package for solar data analysis [@sunpy_project_sunpy_2020].
As the `sunpy` package grew and Python adoption amongst solar physicists increased [see @bobra_survey_2020], the need for instrument-specific packages maintained by instrument teams has become apparent.
As such, the aim of `aiapy` is to provide instrument-specific functionality while being fully interoperable with both the SunPy and broader Astropy [@the_astropy_collaboration_astropy_2018] ecosystems.

![Some examples of the capabilities of the `aiapy` package. The top row shows a cutout of an AIA 171 Å image before (top left) and after (top right) being convolved with the instrument point spread function (PSF). The bottom left panel shows the degradation as a function of time for all seven EUV channels since the launch of SDO. The bottom right panel shows the wavelength response functions for all seven EUV channels.\label{fig:figure1}](figure-1.pdf)

# Package Structure

`aiapy` has three primary subpackages.
\autoref{fig:figure1} shows several examples of the functionality included in these subpackages.
The `aiapy.calibrate` subpackage contains functions for correcting images for telescope degradation, normalizing to the exposure time, and replacing "hot" pixels removed during earlier calibration procedures.
Additionally, this subpackage also contains the `register` function for removing the roll angle, aligning the center of the image with the center of the Sun, and scaling the image to a common resolution across channels.
An image that has been rotated, aligned, and rescaled in this manner is a level 1.5 image.
This function, combined with the ability to update the image metadata that describes the satellite pointing provided by the `update_pointing` function, replicates the commonly-used `aiaprep.pro` IDL procedure, provided in the SolarSoftware ecosystem [@freeland_data_1998], used to align images across channels.

The `aiapy.psf` subpackage includes functions for calculating the instrument point spread function (PSF) for each EUV channel and deconvolving images with the PSF via Richardson-Lucy deconvolution.
Because computing the PSF and performing the deconvolution is computationally expensive for a full-frame AIA image, `aiapy.psf` also includes optional GPU acceleration via `cupy` [@cupy_learningsys2017].
The top row of \autoref{fig:figure1} shows an example of a 171 Å image that has been deconvolved with the PSF.

Finally, the `aiapy.response` subpackage provides the `Channel` class for computing the wavelength response functions for each channel as a function of wavelength.
Optional corrections to the response functions for instrument degradation and channel crosstalk are also provided.
Additionally, individual components of the wavelength response, including the primary and secondary reflectance of each telescope as well as the efficiency of the secondary and focal plane filters, are accessible via the `Channel` class.
The AIA wavelength response functions are described in detail in @boerner_initial_2012.
Combined with atomic data from the CHIANTI atomic database [@dere_chianti_1997], these wavelength response functions can be used to compute the temperature sensitivity of each EUV channel.
The lower right panel of \autoref{fig:figure1} shows the seven EUV wavelength response functions as a function of wavelength.

# Development and Infrastructure

Version `0.3.0` of `aiapy` was released on 6 October 2020 and is available through the [Python Package Index](https://pypi.org/project/aiapy/) via `pip`.
`aiapy` is compatible with Python 3.6+ and is built on top of `sunpy` and `astropy` and utilizes the `drms` package [@glogowski_drms_2019] for retrieving metadata information from the JSOC.
`aiapy` is also part of the SunPy ecosystem and is [a SunPy-affiliated package](https://sunpy.org/project/affiliated#sunpy-affiliated-packages) [See Section 6 of @sunpy_project_sunpy_2020].
As such, it is reviewed on an annual basis to ensure it meets certain established criteria in the categories of functionality, documentation, testing, and community engagement, among others.
Additionally, the code is developed openly [on GitLab](https://gitlab.com/LMSAL_HUB/aia_hub/aiapy) and the documentation is hosted online on [Read the Docs](https://aiapy.readthedocs.io/en/stable/).
`aiapy` includes a comprehensive test suite built on top of the [`pytest` testing framework](https://pytest.org).
The full test suite, including tests for all supported versions of Python, online tests, documentation builds, and code style checks, is run on every single code contribution using the built-in continuous integration pipelines in GitLab and test coverage is monitored and reported using [Codecov](https://codecov.io/).
This test suite is also run weekly to monitor any failures that may occur due to upstream changes in package dependencies.

# Acknowledgements

The authors acknowledge support from NASA's SDO/AIA contract (NNG04EA00C) to the Lockheed Martin Solar and Astrophysics Laboratory.
AIA is an instrument onboard the *Solar Dynamics Observatory*, a mission for NASA's Living With a Star program.
WTB was supported by NASA’s *Hinode* program.
*Hinode* is a Japanese mission developed and launched by ISAS/JAXA with NAOJ as a domestic partner and NASA and STFC (UK) as international partners.
It is operated by these agencies in cooperation with ESA and NSC (Norway).
GC acknowledges support by NASA HSR grant 80NSSC19K0855.
DS is supported by STFC grant ST/S000240/1.
PJW acknowledges support from NASA Contract NAS5-02139 (HMI) to Stanford University.

# References
